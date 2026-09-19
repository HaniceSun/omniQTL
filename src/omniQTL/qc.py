"""Quality-control routines for genotyping arrays and sequencing data.

``ArrayQC`` implements standard PLINK-based sample/variant QC steps for SNP
array genotyping data and is intended to be mixed into a class (such as
:class:`omniQTL.genotyping.Genotyping`) that sets ``self.bfile`` and
``self.output_prefix``. ``SeqQC`` provides QC helpers for sequencing-derived
data (BAM flagstats, QTLtools ``mbv`` genotype concordance, ATAC-seq TSS
enrichment) shared by the caQTL and eQTL pipelines.
"""

from typing import Any

from .utils import *


class ArrayQC:
    """PLINK-based QC steps for genotyping array data.

    This class expects to be mixed into a subclass that sets
    ``self.bfile`` (the input PLINK binary fileset prefix) and
    ``self.output_prefix`` (the prefix used for QC'd output and
    intermediate files), e.g. :class:`omniQTL.genotyping.Genotyping`.
    """

    def check_missingness(self, params: dict[str, Any] = {'mind': 0.05, 'geno': 0.05}) -> None:
        """Remove samples and variants with excess missingness via PLINK.

        Runs ``plink --missing`` for reporting, then filters samples and
        variants using ``--mind``/``--geno`` thresholds, writing the
        filtered fileset to ``self.output_prefix``.

        Args:
            params: Dict with keys ``mind`` (per-sample missingness
                threshold) and ``geno`` (per-variant missingness threshold).
        """
        self.log(self.bfile)
        cmd = f'plink --bfile {self.bfile} --missing --out {self.bfile}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.bfile} --mind {params["mind"]} --geno {params["geno"]} --make-bed --out {self.output_prefix}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        self.log()

    def check_sex(self) -> None:
        """Remove samples that fail PLINK's genetic sex check.

        Runs ``plink --check-sex``, extracts samples flagged ``PROBLEM``,
        and removes them from ``self.output_prefix`` in place.
        """
        cmd = f'plink --bfile {self.output_prefix} --check-sex --out {self.output_prefix}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        sexcheck_failed_file = f'{self.output_prefix}.sexcheck.failed'
        cmd = f'''awk '$5 == "PROBLEM"' {self.output_prefix}.sexcheck > {sexcheck_failed_file}'''
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --remove {sexcheck_failed_file} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def check_heterozygosity(self, params: dict[str, Any] = {'prune': [50, 5, 0.2], 'het': 3}) -> None:
        """Remove heterozygosity-rate outlier samples via PLINK.

        LD-prunes variants, computes per-sample heterozygosity (``F``) with
        ``plink --het``, flags samples more than ``params['het']`` standard
        deviations from the mean ``F``, and removes them from
        ``self.output_prefix`` in place.

        Args:
            params: Dict with keys ``prune`` (a ``[window, step, r2]`` triple
                passed to ``plink --indep-pairwise``) and ``het`` (the
                standard-deviation cutoff for the heterozygosity outlier
                filter).
        """
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
        outliers = df[(df["F"] > mean_F + p * std_F) | (df["F"] < mean_F - p * std_F)]
        outliers[["FID", "IID"]].to_csv(het_failed_file, sep="\t", index=False, header=False)
        cmd = f'plink --bfile {self.output_prefix} --remove {het_failed_file} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def check_relatedness(self, params: dict[str, Any] = {'prune': [50, 5, 0.2], 'rel': 0.185}) -> None:
        """Remove one of each pair of related samples via PLINK.

        LD-prunes variants, computes pairwise IBD with ``plink --genome``,
        and applies ``--rel-cutoff`` to drop samples so that no remaining
        pair exceeds the relatedness threshold, writing the result to
        ``self.output_prefix``.

        Args:
            params: Dict with keys ``prune`` (a ``[window, step, r2]``
                triple passed to ``plink --indep-pairwise``) and ``rel``
                (the relatedness cutoff passed to ``--rel-cutoff``).
        """
        p1, p2, p3 = params['prune']
        cmd = f'plink --bfile {self.output_prefix} --indep-pairwise {p1} {p2} {p3} --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --extract {self.output_prefix}.prune.in --genome --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        cmd = f'plink --bfile {self.output_prefix} --extract {self.output_prefix}.prune.in --rel-cutoff {params["rel"]} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def check_hwe(self, params: dict[str, Any] = {'hwe': 1e-6}) -> None:
        """Remove variants that fail Hardy-Weinberg equilibrium via PLINK.

        Args:
            params: Dict with key ``hwe`` (the HWE exact test p-value
                threshold passed to ``plink --hwe``).
        """
        cmd = f'plink --bfile {self.output_prefix} --hwe {params["hwe"]} --make-bed --out {self.output_prefix}'
        subprocess.run(cmd, shell=True)
        self.log()

    def log(self, bfile: str | None = None) -> None:
        """Print the sample and variant counts of a PLINK fileset.

        Args:
            bfile: PLINK binary fileset prefix to report on. Defaults to
                ``self.output_prefix`` when not given.
        """
        if bfile is None:
            bfile = self.output_prefix
        print('----------------------------------')
        cmd = f'wc -l {bfile}.fam {bfile}.bim'
        subprocess.run(cmd, shell=True)
        print('----------------------------------')


class SeqQC:
    """QC helpers for sequencing-derived (RNA-seq/ATAC-seq) data.

    Covers BAM mapping-rate summaries, QTLtools ``mbv`` genotype
    concordance checks, and ATAC-seq TSS enrichment scoring, along with
    simple bar-plot visualizations for each metric.
    """

    def bam_flagstat(self, out_file: str = 'bam_flagstat.sh', bam_dir: str = 'bams') -> None:
        """Write a shell script running ``samtools flagstat`` on every BAM.

        Args:
            out_file: Path to write the generated shell script.
            bam_dir: Directory containing ``.bam`` files to process.
        """
        bams = sorted([x for x in os.listdir(bam_dir) if x.endswith('.bam')])
        with open(out_file, 'w') as f:
            for bam in bams:
                in_file = os.path.join(bam_dir, bam)
                out_file = in_file.replace('.bam', '.flagstat')
                cmd = f'samtools flagstat {in_file} > {out_file}'
                f.write(cmd + '\n')

    def get_number_mapped_reads(
        self, bam_dir: str = 'bams', out_file: str = 'number_mapped_reads.txt', flag: str = 'primary mapped'
    ) -> None:
        """Tabulate the number of mapped reads per sample from flagstat files.

        Args:
            bam_dir: Directory containing ``.flagstat`` files (as produced
                by :meth:`bam_flagstat`).
            out_file: Path to write the resulting tab-delimited table of
                sample vs. number of reads, sorted ascending.
            flag: The flagstat line label to extract the read count from.
        """
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

    def get_percent_mapped_reads(
        self, bam_dir: str = 'bams', out_file: str = 'percent_mapped_reads.txt', flag: str = 'primary mapped'
    ) -> None:
        """Tabulate the percent of mapped reads per sample from flagstat files.

        Args:
            bam_dir: Directory containing ``.flagstat`` files (as produced
                by :meth:`bam_flagstat`).
            out_file: Path to write the resulting tab-delimited table of
                sample vs. percent mapped, sorted ascending.
            flag: The flagstat line label to extract the percentage from.
        """
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

    def plot_number_mapped_reads(
        self,
        in_file: str = 'number_mapped_reads.txt',
        batch_file: str | None = None,
        figsize: tuple[float, float] = (4, 4),
        base: float = 1e6,
        paired: bool = True,
        cmap: str = 'Set2',
        paired_count: int = 2,
        line_plot: bool = False,
    ) -> None:
        """Plot the number of mapped reads per sample as a bar or line plot.

        Args:
            in_file: Tab-delimited table as produced by
                :meth:`get_number_mapped_reads`.
            batch_file: Optional tab-delimited file with a ``sample`` column
                and a ``batch`` column, merged in to color/facet by batch.
            figsize: Figure size in inches, as ``(width, height)``.
            base: Divisor applied to the read counts (e.g. ``1e6`` for
                millions).
            paired: If True, further divide the (already ``base``-scaled)
                counts by ``paired_count`` to account for paired-end reads.
            cmap: Seaborn palette name used when ``batch_file`` is given.
            paired_count: Number of mates per read pair, used to convert
                paired-end read counts to fragment/pair counts.
            line_plot: If True, draw a line plot instead of a bar plot.
        """
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
            df['number_of_reads'] = df['number_of_reads'].astype(int) / base / paired_count

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if batch:
            if line_plot:
                sns.lineplot(x='sample', y=df.columns[1], data=df, ax=ax, hue='batch', palette=cmap, hue_order=hue_order)
            else:
                sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax, hue='batch', palette=cmap, hue_order=hue_order)
        else:
            if line_plot:
                sns.lineplot(x='sample', y=df.columns[1], data=df, ax=ax)
            else:
                sns.barplot(x='sample', y=df.columns[1], data=df, ax=ax)
        if base == 1e6:
            ax.set_ylabel('Number of mapped reads (million)')
        else:
            ax.set_ylabel('Number of mapped reads')
        ax.set_xticks([])
        plt.tight_layout()
        plt.savefig(out_file)

    def plot_percent_mapped_reads(
        self,
        in_file: str = 'number_mapped_reads.txt',
        batch_file: str | None = None,
        figsize: tuple[float, float] = (4, 4),
        cmap: str = 'Set2',
    ) -> None:
        """Plot the percent of mapped reads per sample as a bar plot.

        Args:
            in_file: Tab-delimited table with a per-sample metric in its
                second column (e.g. as produced by
                :meth:`get_percent_mapped_reads`).
            batch_file: Optional tab-delimited file with a ``sample`` column
                and a ``batch`` column, merged in to color by batch.
            figsize: Figure size in inches, as ``(width, height)``.
            cmap: Seaborn palette name used when ``batch_file`` is given.
        """
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

    def get_mbv_script(
        self,
        bam_dir: str = 'bams',
        vcf_file: str = 'variants.vcf.gz',
        out_file: str = 'run_mbv.sh',
        chrom: str | None = None,
        quality: int = 10,
        QTLtools_env: str | None = 'QTLtools',
    ) -> None:
        """Write a shell script running QTLtools ``mbv`` for every BAM.

        ``mbv`` checks genotype concordance between a BAM file and a VCF,
        useful for detecting sample swaps.

        Args:
            bam_dir: Directory containing ``.bam`` files to process.
            vcf_file: VCF file of genotypes to check concordance against.
            out_file: Path to write the generated shell script.
            chrom: If given, restrict ``mbv`` to this region via ``--reg``.
            quality: Minimum mapping quality passed to
                ``--filter-mapping-quality``.
            QTLtools_env: Conda environment containing the ``QTLtools``
                executable, or None to call it directly.
        """
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

    def get_mbv_results(
        self,
        mbv_dir: str = 'bams',
        out_file: str = 'mbv_merged.txt',
        params: dict[str, Any] = {'perc_het_consistent': 0.8, 'perc_hom_consistent': 0.8},
    ) -> None:
        """Merge QTLtools ``mbv`` output files and flag concordant samples.

        Args:
            mbv_dir: Directory containing ``*_mbv.txt`` files (as produced
                by :meth:`get_mbv_script`).
            out_file: Path to write the merged table of all ``mbv`` results.
                A second file with ``_filtered`` appended to the stem is
                also written, containing only samples passing the
                concordance thresholds.
            params: Dict with keys ``perc_het_consistent`` and
                ``perc_hom_consistent``, the minimum heterozygous/homozygous
                concordance fractions required to keep a sample.
        """
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

    def get_tss_score(
        self, qc_dir: str = 'qc', out_file: str = 'ATACseq_tss_score.txt', score_threshold: float = 4
    ) -> None:
        """Summarize per-sample TSS enrichment scores from ENCODE QC JSONs.

        Args:
            qc_dir: Directory containing ENCODE ATAC-seq pipeline ``.json``
                QC reports.
            out_file: Path to write the tab-delimited table of sample,
                mean TSS enrichment score, and the per-replicate scores. A
                second file with ``_low`` appended to the stem is also
                written, containing only samples below ``score_threshold``.
            score_threshold: Mean TSS enrichment score below which a sample
                is written to the ``_low`` output file.
        """
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

    def plot_tss_score(
        self,
        in_file: str = 'ATACseq_tss_score.txt',
        batch_file: str | None = None,
        figsize: tuple[float, float] = (4, 4),
        cmap: str = 'Set2',
    ) -> None:
        """Plot per-sample TSS enrichment scores as a bar plot.

        Args:
            in_file: Tab-delimited table as produced by :meth:`get_tss_score`.
            batch_file: Optional tab-delimited file with a ``sample`` column
                and a ``batch`` column, merged in to color by batch.
            figsize: Figure size in inches, as ``(width, height)``.
            cmap: Seaborn palette name used when ``batch_file`` is given.
        """
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
