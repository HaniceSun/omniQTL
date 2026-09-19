"""caQTL mapping utilities for bulk ATAC-seq data.

This module runs the ENCODE ATAC-seq pipeline, collects its outputs, builds
consensus/summit-extended peak sets, and quantifies/annotates/TPM-normalizes
ATAC-seq peak counts for caQTL mapping.
"""

from typing import Any

from .utils import *
from .qtl import QTL
from .qc import SeqQC


class CAQTL(QTL, SeqQC):
    """Chromatin accessibility QTL (caQTL) pipeline for bulk ATAC-seq data.

    Wraps the ENCODE ATAC-seq pipeline to go from raw FASTQ files to a
    peak-by-sample count matrix suitable for caQTL mapping: launching the
    pipeline, harvesting its outputs, building consensus and summit-extended
    peak sets, counting reads per peak, annotating peaks with their closest
    gene, and normalizing counts to TPM.

    Attributes:
        atac_seq_pipeline_env: Name of the conda environment used to run the
            ENCODE ATAC-seq pipeline (``caper run``), or ``None`` to run
            without activating a conda environment.
        QTLtools_env: Name of the conda environment used to run QTLtools.
    """

    def __init__(
        self,
        atac_seq_pipeline_env: str = 'atac-seq-pipeline',
        QTLtools_env: str = 'QTLtools',
    ) -> None:
        """Initialize the CAQTL pipeline with its conda environment names.

        Args:
            atac_seq_pipeline_env: Name of the conda environment used to run
                the ENCODE ATAC-seq pipeline.
            QTLtools_env: Name of the conda environment used to run
                QTLtools.
        """
        super().__init__()
        self.atac_seq_pipeline_env = atac_seq_pipeline_env
        self.QTLtools_env = QTLtools_env

    def run_atac_seq_pipeline_encode(
        self,
        fq_file_list: str = 'bulkATACseqStanford.txt',
        wdl_file: str | None = None,
        config_tsv: str | None = None,
        ref_genome: str = 'hg38',
        image: str = 'singularity',
        n_task: int = 4,
        out_dir: str = 'working',
        out_script: str = 'run_pipeline.sh',
    ) -> None:
        """Generate per-sample config JSONs and a shell script to run the ENCODE ATAC-seq pipeline.

        Reads a list of FASTQ file paths, groups them by sample based on
        filename conventions (``_1``/``_R1`` and ``_2``/``_R2`` suffixes for
        paired-end, otherwise single-end), writes a ``caper run`` command per
        sample to ``out_script``, and writes a matching ``atac.wdl`` input
        JSON per sample under ``out_dir``.

        Args:
            fq_file_list: Path to a text file listing FASTQ file paths, one
                per line.
            wdl_file: Path to the ``atac.wdl`` workflow file. Defaults to the
                bundled ``vendor/atac-seq-pipeline/atac.wdl`` file if not
                provided.
            config_tsv: Path to the genome config TSV file for the pipeline.
                Defaults to the bundled
                ``vendor/atac-seq-pipeline/data/<ref_genome>/<ref_genome>.tsv``
                file if not provided.
            ref_genome: Reference genome build used to locate the default
                ``config_tsv`` when ``config_tsv`` is not provided.
            image: Container technology flag passed to ``caper run`` (e.g.
                ``'singularity'``).
            n_task: Maximum number of concurrent tasks passed to
                ``caper run``.
            out_dir: Directory in which per-sample working directories and
                input JSON files are created.
            out_script: Path to the shell script that will be written with
                one ``caper run`` command per sample.

        Raises:
            FileNotFoundError: If the resolved ``wdl_file`` or ``config_tsv``
                does not exist.
        """
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
                    D.setdefault(sample, {'R1': [], 'R2': [], 'R': []})
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

    def get_pipeline_output(
        self,
        in_dir: str = 'working',
        out_file: str = 'get_pipeline_output.sh',
        out_type: list[str] = ['peaks_pvalue', 'peaks_idr', 'peaks_qvalue', 'bams', 'qc'],
        params: dict[str, Any] = {'qvalue_threshold': 0.05},
    ) -> None:
        """Collect and link/filter ENCODE ATAC-seq pipeline outputs by type.

        Scans the pipeline's working directory tree and, for each requested
        output type, either writes ``ln`` commands (to ``out_file``) linking
        the located files into a same-named output directory, or (for
        ``'peaks_qvalue'``) filters narrowPeak files by q-value threshold and
        writes the filtered peaks directly.

        Args:
            in_dir: Root working directory produced by
                :meth:`run_atac_seq_pipeline_encode` to scan recursively for
                pipeline outputs.
            out_file: Path to the shell script that will be written with
                ``ln`` commands for the linked output types.
            out_type: List of output types to collect. Supported values are
                ``'peaks_pvalue'``, ``'peaks_idr'``, ``'peaks_qvalue'``,
                ``'bams'``, and ``'qc'``. For each type, a same-named
                subdirectory is created for its outputs.
            params: Extra parameters. ``'qvalue_threshold'`` is the q-value
                cutoff used to filter narrowPeak files when
                ``'peaks_qvalue'`` is in ``out_type``.
        """
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

    def get_consensus_peaks(
        self,
        in_dir: str = 'peaks_qvalue',
        bed_merged: str | None = None,
        out_file: str = 'ATACseq_consensus_peaks.bed',
        threshold: float = 0.05,
    ) -> None:
        """Build a consensus peak set across samples' narrowPeak files.

        Implemented according to the definition of consensus peak from
        Currin et al., AJHG 2021 (https://pubmed.ncbi.nlm.nih.gov/34038741/).
        If ``bed_merged`` is provided, the function uses the provided bed
        file as the merged peaks instead of merging the peaks from the input
        narrowPeak files; this can be used on the summit-extended peaks
        obtained from :meth:`get_summit_extended_fixed_width_peaks`. If
        ``bed_merged`` is not provided, the function merges the peaks from
        the input narrowPeak files and then applies the consensus peak
        definition on the merged peaks.

        Args:
            in_dir: Directory containing per-sample ``*.narrowPeak.gz``
                files to merge and evaluate for consensus.
            bed_merged: Optional path to a pre-merged bed file to use in
                place of merging ``in_dir``'s narrowPeak files.
            out_file: Path to the output bed file of consensus peaks.
            threshold: Minimum fraction of samples in which a merged peak
                must be observed (via overlap) to be retained as a
                consensus peak.
        """
        fs = sorted([x for x in os.listdir(in_dir) if x.endswith('.narrowPeak.gz')])
        dfs = []
        for f in fs:
            df = pd.read_table(f'{in_dir}/{f}', header=None, sep='\t')
            df = df.iloc[:, 0:3]
            df.columns = ["Chromosome", "Start", "End"]
            dfs.append(pyranges.PyRanges(df))

        if bed_merged is None:
            df_concated = pyranges.concat(dfs)
            df_merged = df_concated.merge()
        else:
            df_merged = pyranges.read_bed(bed_merged)

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

    def get_summit_extended_fixed_width_peaks(
        self,
        in_dir: str = 'peaks_qvalue',
        out_file: str = 'ATACseq_summitExtended_peaks.bed',
        chrom: str | None = None,
        half_width: int = 250,
        params: dict[str, Any] = {'chrom_col': 0, 'start_col': 1, 'summit_col': 9, 'pval_col': 7},
    ) -> None:
        """Build fixed-width, summit-centered peaks via greedy overlap removal.

        Implemented according to the description of the ``merge_peaks``
        function in snapatac2
        (https://scverse.org/SnapATAC2/api/_autosummary/snapatac2.tl.merge_peaks.html).
        For each chromosome, narrowPeak intervals are expanded to
        ``2 * half_width`` around their summit, sorted by significance
        (most significant first, with chromosome/start as tiebreakers for
        reproducibility), and then greedily selected: the most significant
        remaining peak is kept and all peaks overlapping it are discarded,
        repeated until no peaks remain. If the process is too slow on a
        large number of peaks, consider using the ``chrom`` parameter to run
        the function on each chromosome in parallel and then concatenate the
        results. Without providing the ``chrom`` parameter, the function
        loops over all chromosomes one by one.

        Args:
            in_dir: Directory containing per-sample ``*.narrowPeak.gz``
                files to combine and process.
            out_file: Path to the output bed file. When ``chrom`` is
                provided, this is only used to derive the per-chromosome
                output directory/file naming and no combined file is
                written.
            chrom: If provided, restrict processing to this single
                chromosome and write its result to its own file under a
                ``<out_file stem>_chroms`` directory instead of writing (and
                returning) a combined file. If ``None``, all chromosomes are
                processed and concatenated into ``out_file``.
            half_width: Number of base pairs to extend on each side of the
                peak summit.
            params: Column-index mapping into the narrowPeak files.
                ``'chrom_col'``, ``'start_col'``, ``'summit_col'``, and
                ``'pval_col'`` give the 0-based column indices for
                chromosome, peak start, summit offset, and -log10(p-value),
                respectively.
        """
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
            # Step 1: expand peaks, +1 because narrowPeak format is 0-based half-open like bed format
            df['start'] = df[start_col] + df[summit_col] - half_width + 1
            df['end'] = df[start_col] + df[summit_col] + half_width + 1
            # Step 2: sort by significance (smallest p-value first, which is largest -log10(p-value))
            # sorted by chrom and start as well to ensure reproducibility when there are ties in p-value
            df = df.sort_values([pval_col, chrom_col, start_col], ascending=False)
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

    def get_peak_counts(
        self,
        in_file: str = 'ATACseq_consensus_peaks.bed',
        bam_dir: str = 'bams',
        out_file: str = 'featureCounts.sh',
        strand: str = '+',
        n_threads: int = 4,
        min_quality: int = 30,
    ) -> None:
        """Write a SAF annotation file and a featureCounts script for peaks.

        Converts the input bed file of peaks into a SAF (Simplified
        Annotation Format) file, then writes one ``featureCounts`` command
        per BAM file found in ``bam_dir`` to ``out_file``.

        Args:
            in_file: Path to the input bed file of peaks (e.g. consensus or
                summit-extended peaks). A sibling ``.saf`` file is written
                next to it.
            bam_dir: Directory containing per-sample ``*.bam`` files to
                count.
            out_file: Path to the shell script that will be written with one
                ``featureCounts`` command per BAM file.
            strand: Strand value written into the SAF file for every peak.
            n_threads: Number of threads passed to ``featureCounts`` (``-T``).
            min_quality: Minimum mapping quality passed to ``featureCounts``
                (``-Q``).
        """
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

    def merge_counts_tables(
        self,
        counts_tables: str = 'counts_tables.txt',
        out_file: str = 'ATACseq_peakCounts.txt',
    ) -> None:
        """Merge per-sample featureCounts output tables into one count matrix.

        Args:
            counts_tables: Path to a text file listing the paths of
                per-sample featureCounts output files, one per line.
            out_file: Path to the output tab-separated peak-by-sample count
                matrix.
        """
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

    def get_closest_genes(
        self,
        in_file: str = 'ATACseq_peakCounts.txt',
        gtf_file: str = 'Homo_sapiens.GRCh38.115.gtf',
        bed_cols: list[str] = ['chr', 'start', 'end'],
        params: dict[str, Any] = {
            'bed_gene_cols': [2, 3, 4, 5, 0, 1, 6],
            'bed_gene_cols_name': ['chr', 'start', 'end', 'strand', 'GeneID', 'GeneName', 'GeneBiotype'],
            'filter_by_gene_name': True,
        },
    ) -> None:
        """Annotate each peak with its closest gene using ``bedtools closest``.

        Reads a peak count matrix and a gene position/type table (produced
        by ``gtf_to_GenePosType`` in ``utils.py``), converts both to bed
        files, runs ``bedtools closest``, and writes an annotated copy of
        the count matrix with a ``closestGene`` column inserted.

        Args:
            in_file: Path to the tab-separated peak count matrix, with a
                ``peakID`` column formatted as ``chr_start_end``.
            gtf_file: Path to the source GTF file; used to locate the
                companion gene position/type table
                (``<gtf_file>`` with ``.gtf`` replaced by
                ``_GenePosType.txt``).
            bed_cols: Column names to assign to the peak bed file derived
                from splitting ``peakID`` on ``'_'``.
            params: Extra parameters. ``'bed_gene_cols'`` selects and orders
                columns from the gene table into bed format,
                ``'bed_gene_cols_name'`` names those columns, and
                ``'filter_by_gene_name'`` (if true) drops gene entries whose
                ``GeneName`` still starts with ``'ENS'``.

        Raises:
            ValueError: If the companion gene position/type table for
                ``gtf_file`` does not exist.
        """
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
            peakID = df_peak.iloc[i, 0]
            gene = D.get(peakID, ['.'])
            L.append(','.join(sorted(set(gene))))
        df_peak.insert(1, 'closestGene', L)
        df_peak.to_csv(in_file.replace('.txt', '_closestGene.txt'), index=False, sep='\t')
        os.remove(bed_peak)
        os.remove(bed_gene)
        os.remove(bed_closest)

    def counts_to_tpm(
        self,
        in_file: str = 'ATACseq_peakCounts_closestGene.txt',
        norm_base: float = 1e6,
    ) -> None:
        """Normalize a peak count matrix to TPM (transcripts per million).

        Divides each peak's counts by its width (in base pairs) derived from
        its ``peakID``, then rescales each sample's column so its total sums
        to ``norm_base``.

        Args:
            in_file: Path to the tab-separated peak count matrix, with a
                ``peakID`` column formatted as ``chr_start_end``.
            norm_base: Target sum for each sample's normalized column (e.g.
                ``1e6`` for TPM).
        """
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
