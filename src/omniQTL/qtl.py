from .utils import *

class QTL:
    def __init__(self):
        self.window_size_name = {1000:'1k', 5000:'5k', 10000:'10k', 50000:'50k', 100000:'100k', 500000:'500k', 1000000:'1M'}
        self.QTLtools_env = 'QTLtools'

    def subset_rename_samples(self, in_file='ATACseq_peakCounts_closestGene_TPM.txt', sample_file='sample_subset.txt'):
        '''
        subset samples to include those shared between genotype and phenotype data, and rename samples to match those in genotype data if needed. The sample file should be a two-column tab-delimited file, with the first column being the sample names in the phenotype data, and the second column being the corresponding sample names in the genotype data. If the second column is not provided, the sample names will not be renamed. sample order should be the same as in the genotype data, which is important for downstream QTL analysis.
        '''
        df = pd.read_table(in_file, header=0, sep='\t')
        df_sample = pd.read_table(sample_file, header=None, sep='\t')
        flag = all(s in df.columns[2:] for s in df_sample.iloc[:, 0])
        if not flag:
            raise ValueError(f'Some samples in {sample_file} not found in {in_file}')
        df_subset = df[list(df.columns[0:2]) + list(df_sample.iloc[:, 0])]
        if df_sample.shape[1] > 1:
            df_subset.columns = list(df.columns[0:2]) + list(df_sample.iloc[:, 1])
        out_file = in_file.replace('.txt', '_subsetRenamed.txt')
        df_subset.to_csv(out_file, header=True, index=False, sep='\t')

    def filter_phenotype_features(self, tpm_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed.txt', counts_file='ATACseq_peakCounts_closestGene_subsetRenamed.txt', params={'tpm':0.1, 'counts':6, 'sample_percent':0.2}):
        if os.path.exists(tpm_file):
            df_tpm = pd.read_table(tpm_file, header=0, sep='\t')
        else:
            raise FileNotFoundError(f'{tpm_file} not found.')

        if counts_file:
            if os.path.exists(counts_file):
                df_counts = pd.read_table(counts_file, header=0, sep='\t')
            else:
                raise FileNotFoundError(f'{counts_file} not found, using subset_rename_samples to generate it')

        wh = []
        for n in range(df_tpm.shape[0]):
            L = df_tpm.iloc[n, 2:] >= params['tpm']
            flag = False
            if sum(L)/len(L) >= params['sample_percent']:
                flag = True
            wh.append(flag)
        wh = np.array(wh)

        if counts_file:
            wh2 = []
            for n in range(df_counts.shape[0]):
                L = df_counts.iloc[n, 2:] >= params['counts']
                flag = False
                if sum(L)/len(L) >= params['sample_percent']:
                    flag = True
                wh2.append(flag)
            wh2 = np.array(wh2)
            wh = wh & wh2

        feature = 'peak' if tpm_file.find('peak') != -1 else 'gene'
        tpm_file_out = tpm_file.replace('.txt', f'_{feature}Filtered.txt')
        df_tpm = df_tpm.loc[wh, ]
        df_tpm.to_csv(tpm_file_out, header=True, index=False, sep='\t')

    def make_bed_for_QTLtools(self, in_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed_peakFiltered.txt', qtl_type='caQTL', gtf_file='Homo_sapiens.GRCh38.115.gtf'):
        out_file = in_file.replace('.txt', '.bed')
        df = pd.read_table(in_file, header=0, sep='\t')
        df_sample = df.iloc[:, 2:]
        if qtl_type in ['caQTL']:
            df_feat = df.iloc[:, 0].str.split('_', expand=True)
            df_feat.columns = ['#Chr', 'start', 'end']
            df_feat['start'] = df_feat['start'].astype(int)
            df_feat['end'] = df_feat['end'].astype(int)
            df_feat['pid'] = df['peakID'] + '_' + df['closestGene']
            df_feat['gid'] = df['closestGene']
            df_feat['strand'] = '+'
        elif qtl_type in ['eQTL', 'pQTL']:
            D = {}
            gene_table = gtf_file.replace('.gtf', '_GenePosType.txt')
            if os.path.exists(gene_table):
                with open(gene_table) as f:
                    for line in f:
                        line = line.strip()
                        fields = line.split('\t')
                        D[fields[0]] = fields
            else:
                raise ValueError(f'gene table {gene_table} is not found, run gtf_to_GenePosType in utils.py on the gtf file first')

            L = []
            for n in range(df.shape[0]):
                gene_id = df.iloc[n, 0]
                if gene_id in D:
                    gene_name = D[gene_id][1]
                    chrom = D[gene_id][2]
                    if chrom.find('chr') == -1:
                        chrom = 'chr' + chrom
                    strand = D[gene_id][5]
                    if strand == '+':
                        start = int(D[gene_id][3]) - 1
                        end = int(D[gene_id][3])
                    else:
                        start = int(D[gene_id][4]) - 1
                        end = int(D[gene_id][4])
                    gid = gene_id
                    pid = gene_id + '_' + gene_name
                else:
                    chrom, start, end, pid, gid, strand = ['.', '.', '.', '.', '.', '.']
                L.append([chrom, start, end, pid, gid, strand])
            df_feat = pd.DataFrame(L, columns=['#Chr', 'start', 'end', 'pid', 'gid', 'strand'])

        df = pd.concat([df_feat, df_sample], axis=1)
        wh = df['pid'] != '.'
        df = df.loc[wh, ]
        df.sort_values(by=['#Chr', 'start', 'end'], inplace=True)
        df.to_csv(out_file, header=True, index=False, sep='\t')
        cmd = f'bgzip -f {out_file}; tabix -p bed {out_file}.gz'
        subprocess.run(cmd, shell=True)
        print('bed file for QTLtools generated and indexed.')

    def run_PCA_on_bed(self, in_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed_peakFiltered.bed.gz', n_pcs=25, scale=True, center=True):
        out_file = in_file.split('.bed.gz')[0]
        out_file_head = in_file.replace('.bed.gz', f'_PC{n_pcs}.txt')
        cmd = f'QTLtools pca --bed {in_file} --out {out_file}'
        if scale:
            cmd += ' --scale'
        if center:
            cmd += ' --center'
        if self.QTLtools_env is not None:
            cmd = f'conda run -n {self.QTLtools_env} ' + cmd
        subprocess.run(cmd, shell=True)
        cmd = f'head -n {n_pcs + 1} {out_file}.pca > {out_file_head}'
        subprocess.run(cmd, shell=True)
        print('PCA on bed file completed.')

    def get_QTLtools_script(self, pheno_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed_peakFiltered.bed.gz', geno_file='genotype_imputed.vcf.gz', cov_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed_peakFiltered_PC25.txt', out_suffix='PC25', qtl_type='caQTL', qtl_pass=['nominal', 'permute', 'conditional'], n_chunks=30, with_normal=True, with_std_err=True, with_cov=True, window_size=None, fdr_script=None, params={'nominal':1.0, 'permute':1000, 'conditional':0.05, 'seed':42}):
        if n_chunks < 22:
            print('WARNING: QTLtools may not run properly with if chunks fewer than the number of chromosomes')

        pheno_samples = pd.read_table(pheno_file, header=0, sep='\t', nrows=1).columns[6:].values
        geno_samples = []
        with gzip.open(geno_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#CHROM'):
                    geno_samples = np.array(line.split('\t')[9:])
                    break
        cov_samples = pd.read_table(cov_file, header=0, sep=r'\s+', nrows=1).columns[1:].values
        wh1 = sum(geno_samples == pheno_samples) == len(geno_samples)
        wh2 = sum(geno_samples == cov_samples) == len(geno_samples)
        if not wh1 or not wh2:
            raise ValueError(f'Sample names in genotype, phenotype and covariate files do not match')

        for qtl in qtl_pass:
            run_script = f'run_{qtl_type}_{qtl}_{out_suffix}.sh'
            merge_script = f'merge_{qtl_type}_{qtl}_{out_suffix}.sh'
            if window_size is None:
                if qtl_type in ['caQTL']:
                    window_size = 1e3
                elif qtl_type in ['eQTL', 'pQTL']:
                    window_size = 1e6
            window_size = int(window_size)
            window_size_name_str = self.window_size_name.get(window_size, window_size)

            out_dir = f'{qtl_type}_{qtl}-{params[qtl]}_w{window_size_name_str}_{out_suffix}'
            os.makedirs(out_dir, exist_ok=True)

            with open(run_script, 'w') as out:
                for n in range(n_chunks + 1):
                    out_file = f'{out_dir}/chunk_{n}.txt'
                    if qtl in ['nominal', 'permute']:
                        cmd = f'QTLtools cis --chunk {n} {n_chunks} --vcf {geno_file} --bed {pheno_file} --{qtl} {params[qtl]} --window {window_size}'
                    elif qtl == 'conditional':
                        if n == 0:
                            print('Please make sure to run permutation analysis first to get the corresponding mapping file before the conditional analysis.')
                        mapping_file = f'{qtl_type}_permute-{params["permute"]}_w{window_size_name_str}_{out_suffix}.thresholds.txt'
                        cmd = f'QTLtools cis --chunk {n} {n_chunks} --vcf {geno_file} --bed {pheno_file} --mapping {mapping_file} --window {window_size}'

                    if with_cov and os.path.exists(cov_file):
                        cmd += f' --cov {cov_file}'
                    if with_normal:
                        cmd += ' --normal'
                    if with_std_err:
                        cmd += ' --std-err'
                    if qtl in ['permute']:
                        cmd += f' --seed {params["seed"]}'
                    if self.QTLtools_env is not None: 
                        cmd = f'conda run -n {self.QTLtools_env} ' + cmd
                    cmd += f' --out {out_file}'
                    out.write(cmd + '\n')

            with open(merge_script, 'w') as out:
                out_file_merged = f'{out_dir}.txt'
                cmd = f"cat {out_dir}/chunk_*.txt > {out_file_merged}; sed -i 's/ /\\t/g' {out_file_merged}"
                out.write(cmd + '\n')
                cmd = f'gzip -f {out_file_merged}'
                out.write(cmd + '\n')
                if qtl == 'permute':
                    if fdr_script is None:
                        if  self.QTLtools_env is not None:
                            cmd = f'conda run -n {self.QTLtools_env} which QTLtools'
                            sp = subprocess.check_output(cmd, shell=True, text=True)
                            fdr_script = os.path.dirname(sp) + '/scripts/qtltools_runFDR_cis.R'
                        else:
                            raise FileNotFoundError(f'FDR script is needed')
                    if os.path.exists(fdr_script):
                        cmd = f'Rscript {fdr_script} {out_file_merged}.gz {params['conditional']} {out_file_merged.split(".txt")[0]}'
                        if self.QTLtools_env is not None: 
                            cmd = f'conda run -n {self.QTLtools_env} ' + cmd
                    else:
                        raise FileNotFoundError(f'FDR script {fdr_script} not found')
                    out.write(cmd + '\n')

    def get_QTLtools_sig_table(self, in_file, qtl_pass=None, thresholds_file=None, significant_file=None, params={'p_col':'nom_pval', 'padj_col':'adj_beta_pval'}):
        if qtl_pass is None:
            if in_file.find('nominal') != -1:
                qtl_pass = 'nominal'
                if thresholds_file is None:
                    raise ValueError(f'For nominal pass, thresholds_file is needed to filter significant QTLs')
            elif in_file.find('permute') != -1:
                qtl_pass = 'permute'
                if significant_file is None:
                    raise ValueError(f'For permute pass, significant_file is needed to filter significant QTLs')
            elif in_file.find('conditional') != -1:
                qtl_pass = 'conditional'
            else:
                raise ValueError(f'Cannot determine the qtl_pass from the input file name, please specify it in the function argument')
        elif qtl_pass not in ['nominal', 'permute', 'conditional']:
            raise ValueError(f'Invalid qtl_pass value, should be one of nominal, permute or conditional')

        out_file = in_file.split('.txt')[0] + f'_sig.txt.gz'
        if qtl_pass == 'nominal':
            if not os.path.exists(thresholds_file):
                raise FileNotFoundError(f'{thresholds_file} not found')
            df_thresholds = pd.read_table(thresholds_file, header=None, sep=r'\s+')
            df_thresholds.dropna(inplace=True)
            D = dict(zip(df_thresholds.iloc[:, 0], df_thresholds.iloc[:, 1]))

            with gzip.open(in_file, 'rt') as fin, gzip.open(out_file, 'wt') as fout:
                header = fin.readline().strip().split('\t')
                fout.write('\t'.join(header) + '\n')
                p_idx = header.index(params['p_col'])
                for line in fin:
                    fields = line.strip().split('\t')
                    try:
                        p_value = float(fields[p_idx])
                        p_threshold = D.get(fields[0], 0)
                        if p_value < p_threshold:
                            fout.write(line)
                    except:
                        pass
        elif qtl_pass == 'permute':
            if not os.path.exists(significant_file):
                raise FileNotFoundError(f'{significant_file} not found')
            out_file = in_file.split('.txt')[0] + f'_sig.txt.gz'
            df = pd.read_table(in_file, header=0, sep='\t')
            df_sig = pd.read_table(significant_file, header=None, sep=r'\s+')
            wh = df.iloc[:, 0].isin(df_sig.iloc[:, 0])
            df_sub = df.loc[wh, ]
            df_sub.sort_values(by=params['padj_col'], inplace=True)
            df_sub.to_csv(out_file, header=True, index=False, sep='\t')
        elif qtl_pass == 'conditional':
            out_file = in_file.split('.txt')[0] + f'_sig.txt.gz'
            out_file_ranks = in_file.split('.txt')[0] + '_ranks.txt'
            df = pd.read_table(in_file, header=0, sep='\t')
            wh = df['bwd_sig'] == 1
            df_sig = df[wh]
            wh = df_sig['bwd_best_hit'] == 1
            df_best = df_sig.loc[wh]
            df_best.to_csv(out_file, header=True, index=False, sep='\t')
            L = []
            for feature in sorted(df_sig.iloc[:, 0].unique()):
                df_sub = df_sig[df_sig.iloc[:, 0] == feature]
                n = df_sub['rank'].max() + 1
                L.append([feature, n])
            df_ranks = pd.DataFrame(L)
            df_ranks.columns = ['feature', 'n_independent_signals']
            df_ranks.sort_values(by='n_independent_signals', inplace=True)
            df_ranks.to_csv(out_file_ranks, header=True, index=False, sep='\t')
    
    def add_extra_info(self, in_file, vcf_file, do_liftover=True, params={'liftover_from':'hg38', 'liftover_to':'hg19'}):
        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f'{vcf_file} not found')

        converter = None
        if do_liftover:
            converter = liftover.get_lifter(params['liftover_from'], params['liftover_to'])

        D = {}
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = fields[1]
                var_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                maf = '.'
                info = fields[7].split(';')
                for item in info:
                    if item.startswith('MAF='):
                        maf = item.split('=')[-1]
                        break
                ch_pos = '.'
                try:
                    if converter is not None:
                        res = converter[chrom][int(pos)]
                        if res:
                            chrom_new, pos_new = res[0][0], res[0][1]
                        ch_pos = f'{chrom_new}_{pos_new}'
                except Exception as e:
                    print(f'Error in liftover for {chrom}:{pos}, {e}')
                D.setdefault(var_id, [])
                D[var_id].append([alt, ref, maf, ch_pos])

        out_file = in_file.split('.txt')[0] + '_extraInfo'
        out_file_txt = out_file + '.txt'
        with gzip.open(in_file, 'rt') as f, open(out_file, 'w') as out:
            head = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                var_id = fields[7]
                n_var = 0
                if var_id in D:
                    n_var = len(D[var_id])
                    info = D[var_id][0] + [n_var]
                else:
                    info = ['.', '.', '.', '.', n_var]
                out.write(line.strip() + '\t' + '\t'.join([str(x) for x in info]) + '\n')

        with open(out_file_txt, 'w') as out:
            out.write(head + '\t' + '\t'.join(['effective_allele', 'non_effective_allele', 'MAF', f'chr_pos_{params["liftover_to"]}', 'n_var_in_vcf']) + '\n')
        cmd = f'sort -k2,2V -k3,3n {out_file} >> {out_file_txt}; bgzip -f {out_file_txt}; tabix -s 2 -b 3 -e 3 -S 1 {out_file_txt}.gz; rm {out_file}'
        subprocess.run(cmd, shell=True)

