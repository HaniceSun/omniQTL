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

    def get_QTLtools_script(self, pheno_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed_peakFiltered.bed.gz', geno_file='genotype_imputed.vcf.gz', cov_file='ATACseq_peakCounts_closestGene_TPM_subsetRenamed_peakFiltered_PC25.txt', out_suffix='PC25', qtl_type='caQTL', qtl_pass=['nominal', 'permute', 'conditional'], n_chunks=30, with_normal=True, with_std_err=True, with_cov=True, window_size=None, fdr_script=None, params={'nominal':1.0, 'permute':1000, 'conditional':0.05, 'seed':42, 'nominal_p_threshold':5e-8, 'nominal_p_col':12}):
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
        wh1 = sum(geno_samples == cov_samples) == len(geno_samples)
        if not wh1 or not wh2:
            raise ValueError(f'Sample names in genotype, phenotype and covariate files do not match')

        for qtl in qtl_pass:
            run_script = f'run_{qtl_type}_{qtl}_{out_suffix}.sh'
            merge_script = f'merge_{qtl_type}_{qtl}_{out_suffix}.sh'
            if window_size is None:
                if qtl_type in ['caQTL']:
                    window_size = 1e4
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
                        mapping_file = f'{qtl_type}_permute-{params["permute"]}_w{window_size_name_str}.thresholds.txt'
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
                        cmd = f'Rscript {fdr_script} {out_file_merged}.gz {params['conditional']} {out_file_merged}'
                        if self.QTLtools_env is not None: 
                            cmd = f'conda run -n {self.QTLtools_env} ' + cmd
                    else:
                        raise FileNotFoundError(f'FDR script {fdr_script} not found')
                    out.write(cmd + '\n')
                if qtl == 'nominal':
                    out_file_merged_sig = out_file_merged.replace('.txt', '_sig.txt') 
                    p_threshold = params['nominal_p_threshold']
                    p_col = params['nominal_p_col']
                    cmd = f"zcat {out_file_merged}.gz | awk 'NR==1{{print $0}} NR>1{{if(${p_col}<{p_threshold}) print $0}}' > {out_file_merged_sig}"
                    out.write(cmd + '\n')
                    cmd = f'gzip -f {out_file_merged_sig}'
                    out.write(cmd + '\n')
                elif qtl == 'permute':
                    pass

    def add_extra_info(self, in_file, vcf_file):
        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f'{vcf_file} not found')

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
                if var_id not in D:
                    D[var_id] = [alt, ref, maf]

        out_file = in_file.split('.txt')[0] + '_extraInfo.txt.gz'
        with gzip.open(in_file, 'rt') as f, gzip.open(out_file, 'wt') as out:
            head = f.readline().strip()
            out.write(head + '\t' + '\t'.join(['EA', 'NEA', 'MAF']) + '\n')
            for line in f:
                fields = line.strip().split('\t')
                var_id = fields[7]
                if var_id in D:
                    info = D[var_id]
                else:
                    info = ['.', '.', '.']
                out.write(line.strip() + '\t' + '\t'.join(info) + '\n')
