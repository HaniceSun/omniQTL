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

    def get_numer_reads(self, bam_dir='bams', flag='primary mapped'):
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
        out_file = f'{bam_dir}/bam_number_reads.txt'
        df = pd.DataFrame(L, columns=['sample', 'num_reads'])
        df.to_csv(out_file, index=False, sep='\t')

    def plot_number_reads(self, in_file='bams/number_reads.txt'):
        pass

    def get_percent_mapped_reads(self, bam_dir='bams', flag='primary mapped'):
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
        out_file = f'{bam_dir}/bam_percent_mapped_reads.txt'
        df = pd.DataFrame(L, columns=['sample', 'percent_mapped_reads'])
        df.to_csv(out_file, index=False, sep='\t')

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

    def mbv_summary(self, mbv_dir='bams', out_file='mbv_summary.txt'):
        fs = sorted([x for x in os.listdir(mbv_dir) if x.endswith('_mbv.txt')])
        L = []
        for f in fs:
            sample = f.split('_mbv.txt')[0]
            with open(os.path.join(mbv_dir, f)) as f_in:
                for line in f_in:
                    line = line.strip()
                    if line.startswith('MBV'):
                        try:
                            n_total = int(line.split()[2])
                            n_concordant = int(line.split()[3])
                            n_discordant = int(line.split()[4])
                            p_concordant = n_concordant / n_total if n_total > 0 else 0
                            L.append([sample, n_total, n_concordant, n_discordant, p_concordant])
                        except:
                            pass
        df = pd.DataFrame(L, columns=['sample', 'n_total', 'n_concordant', 'n_discordant', 'p_concordant'])
        df.to_csv(out_file, index=False, sep='\t')

    def tss_enrichment_atacseq(self):
        pass

