from .utils import *

class ArrayQC:
    def check_missingness(self, params={'mind': 0.05, 'geno': 0.05}):
        self.log(self.bfile)
        cmd = f'plink --bfile {self.bfile} --missing --out {self.bfile}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.bfile} --mind {params["mind"]} --geno {params["geno"]} --make-bed --out {self.output_prefix}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        self.log()

    def check_sex(self):
        cmd = f'plink --bfile {self.output_prefix} --check-sex --out {self.output_prefix}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        sexcheck_failed_file = f'{self.output_prefix}.sexcheck.failed'
        cmd = f'''awk '$5 == "PROBLEM"' {self.output_prefix}.sexcheck > {sexcheck_failed_file}'''
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --remove {sexcheck_failed_file} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def check_heterozygosity(self, params={'prune':[50, 5, 0.2], 'het':3}):
        p1, p2, p3 = params['prune']
        cmd = f'plink --bfile {self.output_prefix} --indep-pairwise {p1} {p2} {p3} --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --extract {self.output_prefix}.prune.in --het --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

        df = pd.read_csv(f'{self.output_prefix}.het', sep=r"\s+")
        mean_F = df["F"].mean()
        std_F = df["F"].std()
        p = params['het']
        het_failed_file = f'{self.output_prefix}.het.failed'
        outliers = df[(df["F"] > mean_F + p*std_F) | (df["F"] < mean_F - p*std_F)]
        outliers[["FID","IID"]].to_csv(het_failed_file, sep="\t", index=False, header=False)
        cmd = f'plink --bfile {self.output_prefix} --remove {het_failed_file} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def check_relatedness(self, params={'prune':[50, 5, 0.2], 'rel':0.185}):
        p1, p2, p3 = params['prune']
        cmd = f'plink --bfile {self.output_prefix} --indep-pairwise {p1} {p2} {p3} --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --extract {self.output_prefix}.prune.in --genome --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --extract {self.output_prefix}.prune.in --rel-cutoff {params["rel"]} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def check_hwe(self, params={'hwe':1e-6}):
        cmd = f'plink --bfile {self.output_prefix} --hwe {params["hwe"]} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def log(self, bfile=None):
        if bfile is None:
            bfile = self.output_prefix
        print('----------------------------------')
        cmd = f'wc -l {bfile}.fam {bfile}.bim'
        subprocess.run(cmd, shell=True)
        print('----------------------------------')

class SeqQC:
    def bam_flagstat(self, out_file='bam_flagstat.sh', bam_dir='bams'):
        bams = sorted([x for x in os.listdir(bam_dir) if x.endswith('.bam')])
        with open(out_file, 'w') as f:
            for bam in bams:
                in_file = os.path.join(bam_dir, bam)
                out_file = in_file.replace('.bam', '.flagstat')
                cmd = f'samtools flagstat {in_file} > {out_file}'
                f.write(cmd + '\n')

    def get_number_mapped_reads(self, bam_dir='bams', out_file='number_mapped_reads.txt', flag='primary mapped'):
        fs = sorted([x for x in os.listdir(bam_dir) if x.endswith('.flagstat')])
        L = []
        for f in fs:
            sample = f.split('.flagstat')[0]
            with open(os.path.join(bam_dir, f)) as f_in:
                for line in f_in:
                    line = line.strip()
                    if line.find(flag) != -1:
                        try:
                            n = int(line.split()[0])
                            L.append([sample, n])
                        except:
                            pass
        df = pd.DataFrame(L, columns=['sample', 'number_of_reads'])
        df.sort_values(by='number_of_reads', inplace=True)
        df.to_csv(out_file, index=False, sep='\t')

    def get_percent_mapped_reads(self, bam_dir='bams', out_file='percent_mapped_reads.txt', flag='primary mapped'):
        fs = sorted([x for x in os.listdir(bam_dir) if x.endswith('.flagstat')])
        L = []
        for f in fs:
            sample = f.split('.flagstat')[0]
            with open(os.path.join(bam_dir, f)) as f_in:
                for line in f_in:
                    line = line.strip()
                    if line.find(flag) != -1:
                        try:
                            p = float(line.split('%')[0].split('(')[-1])
                            L.append([sample, p])
                        except:
                            pass
        df = pd.DataFrame(L, columns=['sample', 'percent_mapped_reads'])
        df.sort_values(by='percent_mapped_reads', inplace=True)
        df.to_csv(out_file, index=False, sep='\t')

    def plot_number_mapped_reads(self, in_file='number_mapped_reads.txt', batch_file=None, figsize=(6, 4), base=1e6, paired=True, cmap='Set2'):
        out_file = in_file.split('.txt')[0] + '_plot.pdf'
        df = pd.read_table(in_file, sep='\t')
        batch = False
        hue_order = None
        if batch_file is not None:
            if os.path.exists(batch_file):
                batch_df = pd.read_table(batch_file, sep='\t')
                df = df.merge(batch_df, on='sample')
                batch = True
                hue_order = sorted(df['batch'].unique())
        if paired:
            df['number_of_reads'] = df['number_of_reads'].astype(int)/base/2

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if batch:
            sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax, hue='batch', palette=cmap, hue_order=hue_order)
        else:
            sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax)
        if base == 1e6:
            ax.set_ylabel('Number of mapped reads (million)')
        else:
            ax.set_ylabel('Number of mapped reads')
        ax.set_xticks([])
        plt.tight_layout()
        plt.savefig(out_file)

    def plot_percent_mapped_reads(self, in_file='number_mapped_reads.txt', batch_file=None, figsize=(6, 4), cmap='Set2'):
        out_file = in_file.split('.txt')[0] + '_plot.pdf'
        df = pd.read_table(in_file, sep='\t')
        batch = False
        hue_order = None
        if batch_file is not None:
            if os.path.exists(batch_file):
                batch_df = pd.read_table(batch_file, sep='\t')
                df = df.merge(batch_df, on='sample')
                batch = True
                hue_order = sorted(df['batch'].unique())

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if batch:
            sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax, hue='batch', palette=cmap, hue_order=hue_order)
        else:
            sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax)
        ax.set_ylabel('Percent of mapped reads')
        ax.set_xticks([])
        plt.tight_layout()
        plt.savefig(out_file)

    def get_mbv_script(self, bam_dir='bams', vcf_file='variants.vcf.gz', out_file='run_mbv.sh', chrom=None, quality=10, QTLtools_env='QTLtools'):
        bams = sorted([x for x in os.listdir(bam_dir) if x.endswith('.bam')])
        with open(out_file, 'w') as f:
            for bam in bams:
                sample = bam.split('.bam')[0]
                vcf = vcf_file.split('.vcf')[0]
                in_file = os.path.join(bam_dir, bam)
                out = os.path.join(bam_dir, f'{bam}_{vcf}_mbv.txt')
                cmd = f'QTLtools mbv --filter-mapping-quality {quality} --bam {in_file} --vcf {vcf_file} --out {out}'
                if chrom is not None:
                    cmd += f' --reg {chrom}'
                if QTLtools_env is not None:
                    cmd = f'conda run -n {QTLtools_env} ' + cmd
                f.write(cmd + '\n')

    def get_mbv_results(self, mbv_dir='bams', out_file='mbv_merged.txt', params={'perc_het_consistent':0.8, 'perc_hom_consistent':0.8}):
        fs = sorted([x for x in os.listdir(mbv_dir) if x.endswith('_mbv.txt')])
        L = []
        for f in fs:
            sample = f.split('.bam')[0]
            df = pd.read_csv(os.path.join(mbv_dir, f), sep=r'\s+')
            df.insert(0, 'sample', sample)
            L.append(df)
        df = pd.concat(L)
        df.to_csv(out_file, index=False, sep='\t')

        wh1 = df['perc_het_consistent'] > params['perc_het_consistent']
        wh2 = df['perc_hom_consistent'] > params['perc_hom_consistent']
        df_sub = df[wh1 & wh2]
        df_sub.sort_values(by=['perc_het_consistent', 'perc_hom_consistent'], ascending=False, inplace=True)
        out_file = out_file.replace('.txt', '_filtered.txt')
        df_sub.to_csv(out_file, index=False, sep='\t')

    def get_tss_score(self, qc_dir='qc', out_file='ATACseq_tss_score.txt', score_threshold=4):
        out_file_low = out_file.replace('.txt', '_low.txt')
        fs = sorted([x for x in os.listdir(qc_dir) if x.endswith('.json')])
        L = []
        for f in fs:
            sample = f.split('.json')[0]
            try:
                with open(os.path.join(qc_dir, f)) as f_in:
                    data = json.load(f_in)
                    scores = []
                    tss = data['align_enrich']['tss_enrich']
                    for k in tss:
                        scores.append(tss[k]['tss_enrich'])
                    L.append([sample, np.mean(scores), ','.join([str(x) for x in scores])])
            except Exception as e:
                print(f'Error processing {f} {e}')
        if L:
            df = pd.DataFrame(L, columns=['sample', 'mean_tss_score', 'tss_scores'])
            df.sort_values(by='mean_tss_score', inplace=True)
            df.to_csv(out_file, index=False, sep='\t')
            df_sub = df[df['mean_tss_score'] < score_threshold]
            df_sub.to_csv(out_file_low, index=False, sep='\t')

    def plot_tss_score(self, in_file='ATACseq_tss_score.txt', batch_file=None, figsize=(6, 4), cmap='Set2'):
        out_file = in_file.split('.txt')[0] + '_plot.pdf'
        df = pd.read_table(in_file, sep='\t')
        batch = False
        hue_order = None
        if batch_file is not None:
            if os.path.exists(batch_file):
                batch_df = pd.read_table(batch_file, sep='\t')
                df = df.merge(batch_df, on='sample')
                batch = True
                hue_order = sorted(df['batch'].unique())

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if batch:
            sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax, hue='batch', palette=cmap, hue_order=hue_order)
        else:
            sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax)
        ax.set_ylabel('TSS enrichment score')
        ax.set_xticks([])
        plt.tight_layout()
        plt.savefig(out_file)
