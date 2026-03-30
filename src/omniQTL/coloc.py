from .utils import *

class Coloc:
    def __init__(self):
        self.coloc_R_script = BASE / 'scripts/coloc_SuSiE.R'
        if not os.path.exists(self.coloc_R_script):
            raise FileNotFoundError(f"Coloc R script not found at {self.coloc_R_script}.")

    def download_gwas_harmoniser_reference(self, out_dir='harmoniser_reference', url='https://ftp.ebi.ac.uk/pub/databases/gwas/harmonisation_resources'):
        if os.path.exists(out_dir):
            print(f"Directory {out_dir} already exists. Skipping download.")
        else:
            cwd = os.getcwd()
            os.makedirs(out_dir, exist_ok=True)
            os.chdir(out_dir)
            chrs = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY', 'chrMT']
            for ch in chrs:
                for suffix in ['.parquet', '.vcf.gz', '.vcf.gz.tbi']:
                    print(f'Downloading {ch}{suffix}...')
                    cmd = f'wget {url}/homo_sapiens-{ch}{suffix}'
                    subprocess.run(cmd, shell=True)
            for x in ['rsID.sql', 'md5sums.txt']:
                print(f'Downloading {x}...')
                cmd = f'wget {url}/{x}'
                subprocess.run(cmd, shell=True)
            os.chdir(cwd)

    def prepare_gwas_harmoniser_input(self, in_file, params={'chromosome':'Chromsome', 'base_pair_location':'Position', 'effect_allele':'EffectAllele', 'other_allele':'NonEffectAllele', 'p_value':'Pval', 'beta':'Beta', 'standard_error':'SE'}, genome_assembly='37', coordinate_system='1-based', file_type='GWAS-SSF v0.1', is_harmonised='false', is_sorted='false'):
        df = pd.read_table(in_file, header=0, sep=r'\s+', low_memory=False)
        cols = df.columns
        df_out = pd.DataFrame()
        for key, value in params.items():
            if value in cols:
                df_out[key] = df[value]
            else:
                df_out[key] = 'NA'
        for col in cols:
            if col not in params.values():
                df_out[col] = df[col]
        df_out[f'chr_pos_GRCh{genome_assembly}'] = df_out[['chromosome', 'base_pair_location']].astype(str).agg('_'.join, axis=1)
        output_file = in_file.split('.tsv')[0].split('.txt')[0] + '_harmoniser.tsv'
        df_out.to_csv(output_file, sep='\t', index=False)

        output_config = output_file + '-meta.yaml'
        today = datetime.datetime.today().strftime('%Y-%m-%d')
        md5 = subprocess.check_output(f'md5sum {in_file}', shell=True, text=True).strip().split()[0]
        with open(output_config, 'w') as f:
            f.write(f'date_metadata_last_modified: {today}\n')
            f.write(f'genome_assembly: {genome_assembly}\n')
            f.write(f'coordinate_system: {coordinate_system}\n')
            f.write(f'data_file_name: {output_file}\n')
            f.write(f'file_type: {file_type}\n')
            f.write(f'data_file_md5sum: {md5}\n')
            f.write(f'is_harmonised: {is_harmonised}\n')
            f.write(f'is_sorted: {is_sorted}\n')

    def harmonise_gwas_sumstats(self, in_file, reference_dir='harmonisation_reference', from_build='37', to_build='38', cooridnate='1-based', chroms=None, version='v1.1.11', profile='standard,singularity', resume=True):
        if chroms is None:
            df = pd.read_table(in_file, header=0, sep='\t', low_memory=False)
            chroms = list(sorted(df['chromosome'].astype(str).unique()))
        cmd = f'nextflow run  EBISPOT/gwas-sumstats-harmoniser -profile {profile} -r {version} --harm --file {in_file} --ref {reference_dir} --to_build {to_build} --from_build {from_build} --coordinate {cooridnate} --chromlist {",".join(chroms)}'
        if resume:
            cmd += ' -resume'
        print('Ruuning the harmoniser requires >28Gb memory. Please ensure you have sufficient resources before running this command.')
        print(cmd)
        subprocess.run(cmd, shell=True)

    def prepare_coloc_input(self, sumstats1='sumstats_from_QTLtools.txt', sumstats2='sumstats_gwas_harmonised.txt', sumstats1_type='qtl', sumstats2_type='gwas', sumstats1_sample_size=100, sumstats2_sample_size=1000000, sumstats1_study_type='quant', sumstats2_study_type='cc', sumstats1_sig_file='sumstats_from_QTLtools_permute_sig.txt', out_dir=None, bfile_for_ld=None, external_ld=None, sumstats_suffixes=['_ss1', '_ss2'], params1={'var_key':['var_id'], 'feature_id':'phe_id', 'pos_col':'var_from', 'beta_col':'slope', 'se_col':'slope_se', 'maf_col':'MAF'}, params2={'var_key':['rsid'], 'pos_col':'base_pair_location', 'beta_col':'beta', 'se_col':'standard_error', 'maf_col':'MAF'}):
        if out_dir is None:
            out_dir = sumstats1.split('.txt')[0].split('.tsv')[0] + '_' + sumstats2.split('.txt')[0].split('.tsv')[0] + '_coloc'
        os.makedirs(out_dir, exist_ok=True)
        if bfile_for_ld is not None and external_ld is not None:
            raise ValueError("Please provide either bfile_for_ld or external_ld, not both.")
        if bfile_for_ld is not None:
            if not os.path.exists(bfile_for_ld + '.bed'):
                raise FileNotFoundError(f"Genotype data not found at {bfile_for_ld}.")
        elif external_ld is not None:
            if not os.path.exists(external_ld):
                raise FileNotFoundError(f"External LD file not found at {external_ld}.")

        df1 = pd.read_table(sumstats1, sep='\t', header=0)
        if sumstats1_sig_file is not None and os.path.exists(sumstats1_sig_file):
            df_sig = pd.read_table(sumstats1_sig_file, header=0, sep=r'\s+')
            wh = df1.iloc[:, 0].isin(df_sig.iloc[:, 0])
            df1 = df1[wh].copy()

        df1['var_key'] = df1[params1['var_key']].apply(lambda x: '_'.join(x.astype(str)), axis=1)

        df2 = pd.read_table(sumstats2, sep='\t', header=0)
        df2['var_key'] = df2[params2['var_key']].apply(lambda x: '_'.join(x.astype(str)), axis=1)
        df1.columns = [x + sumstats_suffixes[0] for x in df1.columns]
        df2.columns = [x + sumstats_suffixes[1] for x in df2.columns]

        features = sorted(df1[params1['feature_id'] + sumstats_suffixes[0]].unique())
        for feature in features:
            df1_sub = df1[df1[params1['feature_id'] + sumstats_suffixes[0]] == feature]
            if sumstats2_type == 'gwas':
                df2_sub = df2
            elif sumstats2_type == 'qtl':
                df2_sub = df2[df2[params2['feature_id'] + sumstats_suffixes[1]] == feature]
            df_merged = pd.merge(df1_sub, df2_sub, left_on='var_key' + sumstats_suffixes[0], right_on='var_key' + sumstats_suffixes[1]).sort_values(by=params1['pos_col'] + sumstats_suffixes[0])
            if df_merged.shape[0]:
                print(df_merged)
                df1x = pd.DataFrame() 
                df2x = pd.DataFrame() 
                df1x['snp'] = df_merged['var_key' + sumstats_suffixes[0]]
                df1x['pos'] = df_merged[params1['pos_col'] + sumstats_suffixes[0]]
                df1x['beta'] = df_merged[params1['beta_col'] + sumstats_suffixes[0]]
                df1x['varbeta'] = df_merged[params1['se_col'] + sumstats_suffixes[0]] ** 2
                df1x['type'] = sumstats1_study_type
                df1x['N'] = sumstats1_sample_size
                if sumstats1_study_type == 'quant':
                    df1x['MAF'] = df_merged[params1['maf_col'] + sumstats_suffixes[0]]

                df2x['snp'] = df_merged['var_key' + sumstats_suffixes[1]]
                df2x['pos'] = df_merged[params2['pos_col'] + sumstats_suffixes[1]]
                df2x['beta'] = df_merged[params2['beta_col'] + sumstats_suffixes[1]]
                df2x['varbeta'] = df_merged[params2['se_col'] + sumstats_suffixes[1]] ** 2
                df2x['type'] = sumstats2_study_type
                df2x['N'] = sumstats2_sample_size
                if sumstats2_study_type == 'quant':
                    df2x['MAF'] = df_merged[params2['maf_col'] + sumstats_suffixes[1]]

                output_file1 = os.path.join(out_dir, f'{feature}{sumstats_suffixes[0]}.txt')
                output_file2 = os.path.join(out_dir, f'{feature}{sumstats_suffixes[1]}.txt')
                output_snps = os.path.join(out_dir, f'{feature}_snps.txt')
                output_ld = os.path.join(out_dir, f'{feature}')

                df1x.to_csv(output_file1, sep='\t', index=False)
                df2x.to_csv(output_file2, sep='\t', index=False)
                df1x['snp'].to_csv(output_snps, sep='\t', index=False, header=False)

                if bfile_for_ld is not None:
                    cmd = f'plink --bfile {bfile_for_ld} --r square --extract {output_snps} --out {output_ld}'
                    print(cmd)
                    subprocess.run(cmd, shell=True, check=True)
                elif external_ld is not None:
                    print('Using external LD file, to be implemented.')

    def get_coloc_script(self, in_dir, out_script='run_coloc.sh', R_env='QTLtools'):
        ld_files = sorted([f for f in os.listdir(in_dir) if f.endswith('.ld')])
        with open(out_script, 'w') as f:
            for ld in ld_files:
                ss1 = ld.split('.ld')[0] + '_ss1.txt'
                ss2 = ld.split('.ld')[0] + '_ss2.txt'
                cmd = f'Rscript {self.coloc_R_script} {in_dir}/{ss1} {in_dir}/{ss2} {in_dir}/{ld}'
                if R_env is not None:
                    cmd = f'conda run -n {R_env} ' + cmd
                f.write(cmd + '\n')

    def merge_coloc_results(self, in_dir, parmas={'PP.H4.abf':0.8}):
        out_file = in_dir + '_results.txt'
        out_file_sub = in_dir + '_results_H4.txt'
        coloc_files = sorted([f for f in os.listdir(in_dir) if f.endswith('_coloc.txt')])
        df_list = []
        for coloc in coloc_files:
            df = pd.read_table(os.path.join(in_dir, coloc), sep='\t')
            df['feature'] = coloc.split('_coloc')[0]
            df_list.append(df)
        if not df_list:
            print("No coloc result files found in the directory.")
            return
        df_merged = pd.concat(df_list, axis=0)
        k = 'PP.H4.abf'
        v = parmas['PP.H4.abf']
        df_sub = df_merged[df_merged[k] >= v]
        df_merged.to_csv(out_file, sep='\t', index=False)
        df_sub.to_csv(out_file_sub, sep='\t', index=False)
