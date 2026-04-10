from .utils import *
from .qtl import QTL
from .qc import SeqQC

class CAQTL(QTL, SeqQC):
    def __init__(self, atac_seq_pipeline_env='atac-seq-pipeline', QTLtools_env='QTLtools'):
        super().__init__()
        self.atac_seq_pipeline_env = atac_seq_pipeline_env
        self.QTLtools_env = QTLtools_env

    def run_atac_seq_pipeline_encode(self, fq_file_list='bulkATACseqStanford.txt', wdl_file=None, config_tsv=None, ref_genome='hg38', image='singularity', n_task=4, out_dir='working', out_script = 'run_pipeline.sh'):
        out_dir = os.path.abspath(out_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        if wdl_file is None:
            wdl_file = BASE.parent.parent / 'vendor' / 'atac-seq-pipeline' / 'atac.wdl'
        if config_tsv is None:
            config_tsv = BASE.parent.parent / 'vendor' / 'atac-seq-pipeline' / 'data/' / ref_genome / f'{ref_genome}.tsv'
        if not os.path.exists(wdl_file) or not os.path.exists(config_tsv):
            raise FileNotFoundError('WDL file or config TSV file not found.')

        D = {}
        with open(fq_file_list) as f:
            for line in f:
                line = line.strip()
                if os.path.isfile(line):
                    sample = line.split('/')[-1].split('_')[0]
                    D.setdefault(sample, {'R1': [], 'R2': [], 'R':[]})
                    if line.endswith('_1.fq.gz') or line.endswith('_R1.fq.gz'):
                        D[sample]['R1'].append(line)
                    elif line.endswith('_2.fq.gz') or line.endswith('_R2.fq.gz'):
                        D[sample]['R2'].append(line)
                    else:
                        D[sample]['R'].append(line)
                else:
                    print(f'File {line} not found, skipping.')

        with open(out_script, 'w') as out:
            for sample in sorted(D):
                out_dir_sample = os.path.join(out_dir, sample)
                os.makedirs(out_dir_sample, exist_ok=True)
                out_json = f'{out_dir}/{sample}.json'
                out.write(f'cd {out_dir_sample}\n')
                cmd = f'caper run {wdl_file} -i {out_json} --{image} --max-concurrent-tasks {n_task}'
                if self.atac_seq_pipeline_env is not None:
                    cmd = f'conda run -n {self.atac_seq_pipeline_env} ' + cmd
                out.write(cmd + '\n')

                config = {}
                config['atac.pipeline_type'] = 'atac'
                config['atac.genome_tsv'] = str(config_tsv)
                config['atac.auto_detect_adapter'] = True
                config['atac.enable_xcor'] = False
                config['atac.title'] = sample

                R1 = sorted(D[sample]['R1'])
                R2 = sorted(D[sample]['R2'])
                R = sorted(D[sample]['R'])
                if R1 and R2 and not R:
                    print(f'{sample}: paired-end')
                    config['atac.paired_end'] = True
                    config['atac.fastqs_rep1_R1'] = D[sample]['R1']
                    config['atac.fastqs_rep1_R2'] = D[sample]['R2']

                elif R and not R1 and not R2:
                    print(f'{sample}: single-end')
                    config['atac.paired_end'] = False
                    config['atac.fastqs_rep1_R1'] = D[sample]['R']

                else:
                    print(f'check fastq files for sample {sample}.')

                with open(out_json, 'w') as out2:
                    json.dump(config, out2, indent=4)

    def get_pipeline_output(self, in_dir='working', out_file='get_pipeline_output.sh', out_type=['peaks_pvalue', 'peaks_idr', 'peaks_qvalue', 'bams', 'qc'], params={'qvalue_threshold': 0.05}):
        fs = glob.glob(f'{in_dir}/**/*', recursive=True)
        with open(out_file, 'w') as out:
            for ot in out_type:
                print(f'Processing output: {ot}', flush=True)
                os.makedirs(ot, exist_ok=True)
                if ot == 'peaks_pvalue':
                    for f in fs:
                        if f.endswith('.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz') and f.find('glob-') == -1:
                            sample = f.split('/')[1]
                            out.write(f'ln {f} {ot}/{sample}.narrowPeak.gz\n')

                if ot == 'peaks_qvalue':
                    for f in fs:
                        if f.endswith('.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak.gz') and f.find('glob-') == -1:
                            sample = f.split('/')[1]
                            df = pd.read_table(f, header=0, sep='\t')
                            df.dropna(subset=df.columns[8], inplace=True)
                            wh = np.power(10, -df.iloc[:, 8].values) < params['qvalue_threshold']
                            df = df[wh]
                            out_f = f'{ot}/{sample}.narrowPeak.gz'
                            df.to_csv(out_f, header=False, index=False, sep='\t')

                elif ot == 'peaks_idr':
                    for f in fs:
                        if f.endswith('idr.conservative_peak.narrowPeak.gz') and f.find('glob-') == -1:
                            sample = f.split('/')[1]
                            out.write(f'ln {f} {ot}/{sample}.narrowPeak.gz\n')
                elif ot == 'bams':
                    for f in fs:
                        if f.endswith('.bam') and f.find('glob-') == -1 and f.find('/call-align/') != -1:
                            sample = f.split('/')[1]
                            out.write(f'ln {f} {ot}/{sample}.bam\n')
                            out.write(f'ln {f}.bai {ot}/{sample}.bam.bai\n')
                elif ot == 'qc':
                    for f in fs:
                        if f.endswith('qc.json') and f.find('glob-') == -1:
                            sample = f.split('/')[1]
                            out.write(f'ln {f} {ot}/{sample}.json\n')

    def get_consensus_peaks(self, in_dir='peaks_qvalue', out_file='ATACseq_consensus_peaks.bed', threshold=0.05):
        fs = sorted([x for x in os.listdir(in_dir) if x.endswith('.narrowPeak.gz')])
        dfs = []
        for f in fs:
            df = pd.read_table(f'{in_dir}/{f}', header=None, sep='\t')
            df = df.iloc[:, 0:3]
            df.columns = ["Chromosome","Start","End"]
            dfs.append(pyranges.PyRanges(df))
        df_concated = pyranges.concat(dfs)
        df_merged = df_concated.merge()

        D = {}
        for n in range(len(dfs)):
            df_overlapped = df_merged.overlap(dfs[n])
            keys = df_overlapped.df.iloc[:, 0:3].astype(str).agg('_'.join, axis=1)
            for k in keys:
                D.setdefault(k, 0)
                D[k] += 1

        df = df_merged.df.copy()
        df['peakID'] = df.iloc[:, 0:3].astype(str).agg('_'.join, axis=1)
        df['count'] = df['peakID'].map(D).fillna(0).astype(int)
        wh = df['count'] >= len(dfs) * threshold
        df = df[wh].iloc[:, 0:3]
        df.to_csv(out_file, header=False, index=False, sep='\t')

    def get_summit_extended_fixed_width_peaks(self, in_dir='peaks_qvalue', out_file='ATACseq_summitExtended_peaks.bed', chrom=None, half_width=250, params={'chrom_col':0, 'start_col':1, 'summit_col':9, 'pval_col':7}):
        '''
        if the process is too slow on a large number of peaks, consider using the chrom parameter to run the function on each chromosome in parallel and then merge the results. Without the chrom parameter, the function will loop on all the chroms one by one.
        '''

        fs = sorted([x for x in os.listdir(in_dir) if x.endswith('.narrowPeak.gz')])
        chrom_col = params['chrom_col']
        start_col = params['start_col']
        summit_col = params['summit_col']
        pval_col = params['pval_col']

        dfs = []
        for f in fs:
            df = pd.read_table(f'{in_dir}/{f}', header=None, sep='\t')
            dfs.append(df)
        df_concated = pd.concat(dfs)

        if chrom is None:
            chroms = sorted(df_concated[chrom_col].unique())
        else:
            chroms = [chrom]
            out_dir = out_file.split('.bed')[0] + f'_chroms'
            os.makedirs(out_dir, exist_ok=True)

        L = []
        for ch in chroms:
            print(ch, flush=True)
            df = df_concated[df_concated[chrom_col] == ch]
            # Step 1: expand peaks
            df['start'] = df[start_col] + df[summit_col] - half_width + 1
            df['end'] = df[start_col] + df[summit_col] + half_width + 1
            # Step 2: sort by significance (smallest p-value first, which is largest -log10(p-value))
            df = df.sort_values(pval_col, ascending=False)
            selected = []
            while len(df) > 0:
                # Step 3: pick most significant peak
                top = df.iloc[0, [chrom_col, -2, -1]]
                selected.append(top)
                # Step 4: remove overlapping peaks
                non_overlap_mask = (df['start'] > top['end']) | (df['end'] < top['start'])
                df = df[non_overlap_mask]
            df_selected = pd.DataFrame(selected)
            if chrom is not None:
                df_selected.columns = ['chr', 'start', 'end']
                df_selected['start'] = df_selected['start'] - 1
                df_selected = df_selected[df_selected['start'] >= 0]
                df_selected.sort_values(by=['chr', 'start', 'end'], inplace=True)
                out_file_chrom = out_dir + '/' + out_file.split('.bed')[0] + f'_{chrom}.bed'
                df_selected.to_csv(out_file_chrom, header=False, index=False, sep='\t')
            else:
                L.append(df_selected)
        if L:
            df_merged = pd.concat(L)
            df_merged.columns = ['chr', 'start', 'end']
            df_merged['start'] = df_merged['start'] - 1
            df_merged = df_merged[df_merged['start'] >= 0]
            df_merged.sort_values(by=['chr', 'start', 'end'], inplace=True)
            df_merged.to_csv(out_file, header=False, index=False, sep='\t')

    def get_peak_counts(self, in_file='ATACseq_consensus_peaks.bed', bam_dir='bams', out_file='featureCounts.sh', strand='+', n_threads=4, min_quality=30):
        saf = in_file.replace('.bed', '.saf')
        with open(in_file) as f, open(saf, 'w') as out:
            for line in f:
                line = line.strip()
                chrom, start, end = line.split('\t')
                start = str(int(start) + 1)
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom 
                k = '_'.join([chrom, start, end])
                out.write('\t'.join([k, chrom, start, end, strand]) + '\n')

        fs = sorted([x for x in os.listdir(bam_dir) if x.endswith('.bam')])
        with open(out_file, 'w') as out:
            for f in fs:
                sample = f.split('.bam')[0]
                bam = bam_dir + '/' + f
                counts_file = f'{bam_dir}/{saf.split(".saf")[0]}_{sample}_peakCounts.txt'
                cmd = f'featureCounts -T {n_threads} -O -p -F SAF -a {saf} -Q {min_quality} -s 0 -o {counts_file} {bam}'
                out.write(cmd + '\n')

    def merge_counts_tables(self, counts_tables='counts_tables.txt', out_file='ATACseq_peakCounts.txt'):
        df_tables = pd.read_table(counts_tables, header=None, sep='\t')
        L = []
        for n in range(df_tables.shape[0]):
            f = df_tables.iloc[n, 0]
            sample = f.split('/')[-1].split('_peakCounts.txt')[0].split('_peaks_')[-1]
            df = pd.read_table(f, header=0, sep='\t', comment='#')
            if n == 0:
                df_peak = df.iloc[:, 0:1]
                df_peak.columns = ['peakID']
                L.append(df_peak)
            df_sample = df.iloc[:, 6:7]
            df_sample.columns = [sample]
            L.append(df_sample)
        df = pd.concat(L, axis=1)
        df.to_csv(out_file, index=False, sep='\t')

    def get_closest_genes(self, in_file='ATACseq_peakCounts.txt', gtf_file='Homo_sapiens.GRCh38.115.gtf', bed_cols=['chr', 'start', 'end'], params={'bed_gene_cols':[2, 3, 4, 5, 0, 1, 6], 'bed_gene_cols_name':['chr', 'start', 'end', 'strand', 'GeneID', 'GeneName', 'GeneBiotype'], 'filter_by_gene_name':True}):
        gene_table = gtf_file.replace('.gtf', '_GenePosType.txt')
        if not os.path.exists(gene_table):
            raise ValueError(f'gene table {gene_table} is not found, run gtf_to_GenePosType in utils.py on the gtf file first')
        bed_peak = in_file.split('.txt')[0] + '_tmp.bed'
        bed_gene = os.path.basename(gene_table).split('.txt')[0] + '_tmp.bed'
        df_peak = pd.read_table(in_file, header=0, sep='\t')
        df_gene = pd.read_table(gene_table, header=None, sep='\t')
        bed = df_peak['peakID'].str.split('_', expand=True)
        bed.columns = bed_cols
        bed['start'] = bed['start'].astype(int)
        bed['end'] = bed['end'].astype(int)
        bed.sort_values(by=bed_cols, inplace=True)

        bed2 = df_gene[params['bed_gene_cols']].copy()
        bed2.columns = params['bed_gene_cols_name']
        if bed2['chr'].iloc[0].find('chr') == -1:
            bed2['chr'] = 'chr' + bed2['chr']
        if params['filter_by_gene_name']:
            wh = bed2['GeneName'].str.startswith('ENS')
            bed2 = bed2.loc[~wh]
        bed2.sort_values(by=bed_cols, inplace=True)
        bed.to_csv(bed_peak, header=False, index=False, sep='\t')
        bed2.to_csv(bed_gene, header=False, index=False, sep='\t')

        bed_closest = in_file.split('.txt')[0] + '_closest.txt'
        cmd = f'bedtools closest -D ref -a {bed_peak} -b {bed_gene} > {bed_closest}'
        print('Running command:', cmd)
        subprocess.run(cmd, shell=True)

        df_closest = pd.read_table(bed_closest, header=None, sep='\t')
        D = {}
        for i in range(df_closest.shape[0]):
            chrom = df_closest.iloc[i, 0]
            start = df_closest.iloc[i, 1]
            end = df_closest.iloc[i, 2]
            peakID = f'{chrom}_{start}_{end}'
            geneID = df_closest.iloc[i, 7]
            geneName = df_closest.iloc[i, 8]
            geneBiotype = df_closest.iloc[i, 9]
            distance = df_closest.iloc[i, -1]
            D.setdefault(peakID, [])
            D[peakID].append(geneName)

        L = []
        for i in range(df_peak.shape[0]):
            peakID =  df_peak.iloc[i, 0]
            gene = D.get(peakID, ['.'])
            L.append(','.join(sorted(set(gene))))
        df_peak.insert(1, 'closestGene', L)
        df_peak.to_csv(in_file.replace('.txt', '_closestGene.txt'), index=False, sep='\t')
        os.remove(bed_peak)
        os.remove(bed_gene)
        os.remove(bed_closest)

    def counts_to_tpm(self, in_file='ATACseq_peakCounts_closestGene.txt', norm_base=1e6):
        out_file = in_file.replace('.txt', '_TPM.txt')
        df = pd.read_table(in_file, header=0, sep='\t')
        mat = df.iloc[:, 2:]
        mat_peak = df.iloc[:, 0:2]

        L = []
        for n in range(df.shape[0]):
            fields = df.iloc[n, 0].split('_')
            lh = int(fields[2]) - int(fields[1]) + 1
            L.append(lh)
        Length = np.array(L)

        matTotalRaw = mat.sum(axis=0)
        print(f'Total Reads (million):\n{matTotalRaw/norm_base}')
        matLength = (mat.T/Length).T
        matTotal = matLength.sum(axis=0)
        M = matLength/matTotal*norm_base
        df = pd.concat([mat_peak, M], axis=1)
        df.to_csv(out_file, header=True, index=False, sep='\t', float_format='%.4f')

if __name__ == '__main__':
    qtl = CAQTL()
    qtl.run_atac_seq_pipeline_encode()
