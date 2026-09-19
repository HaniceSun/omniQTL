"""RNA-seq processing utilities for eQTL/sQTL analyses.

This module provides the :class:`EQTL` class, which generates the shell
commands used to run a pancreatic islet RNA-seq pipeline: STAR genome
indexing and read mapping, gene- and exon-level read counting with
featureCounts, intron-retention/exon-ratio (IR/ER) and percent-spliced-in
(PSI) calculations for splicing (sQTL) analyses, and TPM normalization of
count tables.
"""

from typing import Any

from .utils import *
from .qtl import QTL
from .qc import SeqQC


class EQTL(QTL, SeqQC):
    """Generate shell commands for an RNA-seq eQTL/sQTL processing pipeline.

    This class wraps STAR (genome indexing and read alignment),
    featureCounts (gene- and exon-level counting), and a set of helper
    methods that derive intron-retention/exon-ratio (IR/ER) counts,
    percent-spliced-in (PSI) values, merged count tables, gene/exon name
    annotations, and TPM-normalized expression tables. Most methods do not
    execute commands directly; they write shell scripts (or other output
    tables) to disk for submission on the cluster.

    Attributes:
        QTLtools_env: Name of the QTLtools conda/module environment used
            by inherited QTL functionality.
        samples: List of sample names discovered from the input fastq
            files. Only set after :meth:`ranseq_mapping` has been called.
        paired_end: Whether the fastq files processed by
            :meth:`ranseq_mapping` are paired-end. Only set after
            :meth:`ranseq_mapping` has been called.
    """

    def __init__(self, QTLtools_env: str = 'QTLtools') -> None:
        """Initialize the EQTL pipeline helper.

        Args:
            QTLtools_env: Name of the QTLtools environment/module to use
                for QTLtools-related commands.
        """
        super().__init__()
        self.QTLtools_env = QTLtools_env

    def genome_indexing(
        self,
        fa_file: str = 'GRCh38.fa',
        gtf_file: str = 'GRCh38.115.gtf',
        genome_dir: str = 'GRCh38_STAR',
        out_file: str = 'indexing_genome.sh',
        n_threads: int = 4,
        sjdbOverhang: int = 100,
    ) -> None:
        """Write a shell script that builds a STAR genome index.

        Args:
            fa_file: Path to the reference genome FASTA file.
            gtf_file: Path to the reference annotation GTF file.
            genome_dir: Output directory for the STAR genome index.
            out_file: Path of the shell script to write.
            n_threads: Number of threads passed to STAR (--runThreadN).
            sjdbOverhang: Splice junction database overhang length passed
                to STAR (--sjdbOverhang), typically read length minus 1.
        """
        cmd = (
            f'STAR --runThreadN {n_threads} --runMode genomeGenerate --genomeDir {genome_dir} '
            f'--genomeFastaFiles {fa_file} --sjdbGTFfile {gtf_file} --sjdbOverhang {sjdbOverhang}'
        )
        with open(out_file, 'w') as out:
            out.write(cmd + '\n')

    def ranseq_mapping(
        self,
        out_file: str = 'mapping.sh',
        fq_dir: str = '.',
        flag: str = 'NH HI AS nM XS',
        n_threads: int = 4,
        genome_dir: str = 'GRCh38_STAR',
    ) -> None:
        """Write a shell script that aligns fastq files with STAR.

        Scans ``fq_dir`` for ``.fq.gz``/``.fastq.gz`` files and groups them
        into samples by stripping the trailing ``_1``, ``_2``, ``_R1``, or
        ``_R2`` suffix from the file name. Files containing ``_1`` or
        ``_R1`` are treated as mate 1, files containing ``_2`` or ``_R2``
        are treated as mate 2, and any other files are treated as
        single-end reads. Sets ``self.samples`` (list of discovered sample
        names, in the order STAR commands are written) and
        ``self.paired_end`` (True/False, based on the last sample
        processed) as a side effect.

        Args:
            out_file: Path of the shell script to write.
            fq_dir: Directory containing the input fastq files.
            flag: Space-separated list of extra SAM attributes passed to
                STAR (--outSAMattributes).
            n_threads: Number of threads passed to STAR (--runThreadN).
            genome_dir: Path to the STAR genome index directory.

        Raises:
            ValueError: If a sample has fastq files that are neither a
                consistent set of paired mates nor single-end reads.
        """
        self.samples = []
        fqs = sorted(
            [x for x in os.listdir(fq_dir) if x.endswith('.fq.gz') or x.endswith('.fastq.gz')]
        )
        D: dict[str, dict[str, list[str]]] = {}
        for fq in fqs:
            sample = fq.split('_1')[0].split('_2')[0].split('_R1')[0].split('_R2')[0]
            D.setdefault(sample, {})
            D[sample].setdefault('fq', [])
            D[sample].setdefault('fq1', [])
            D[sample].setdefault('fq2', [])
            if fq_dir != '.':
                fq = f'{fq_dir}/{fq}'
            if fq.find('_1') != -1 or fq.find('_R1') != -1:
                D[sample]['fq1'].append(fq)
            elif fq.find('_2') != -1 or fq.find('_R2') != -1:
                D[sample]['fq2'].append(fq)
            else:
                D[sample]['fq'].append(fq)

        with open(out_file, 'w') as outfile:
            for sample in D:
                if (
                    len(D[sample]['fq1']) > 0
                    and len(D[sample]['fq2']) > 0
                    and len(D[sample]['fq1']) == len(D[sample]['fq2'])
                ):
                    fs = ','.join(D[sample]['fq1']) + ' ' + ','.join(D[sample]['fq2'])
                    self.paired_end = True
                elif len(D[sample]['fq']) > 0:
                    fs = ','.join(D[sample]['fq'])
                    self.paired_end = False
                else:
                    raise ValueError(f'check fastq file names of {sample}')
                cmd = (
                    f'STAR --runThreadN {n_threads} --runMode alignReads --genomeDir {genome_dir} '
                    f'--readFilesIn {fs} --readFilesCommand zcat --outFileNamePrefix {sample}_ '
                    f'--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx '
                    f'--outSAMattributes {flag}; samtools index {sample}_Aligned.sortedByCoord.out.bam'
                )
                outfile.write(cmd + '\n')
                self.samples.append(sample)

    def counting_genes(
        self,
        out_file: str = 'counting_genes.sh',
        strand: int = 0,
        min_quality: int = 30,
        n_threads: int = 4,
        gtf_file: str = 'GRCh38.115.gtf',
    ) -> None:
        """Write a shell script that runs featureCounts at the gene level.

        Iterates over ``self.samples`` (set by :meth:`ranseq_mapping`) and
        writes one featureCounts command per sample, using the STAR
        coordinate-sorted BAM file produced by :meth:`ranseq_mapping`.

        Args:
            out_file: Path of the shell script to write.
            strand: Strandedness passed to featureCounts (-s).
            min_quality: Minimum mapping quality passed to featureCounts (-Q).
            n_threads: Number of threads passed to featureCounts (-T).
            gtf_file: Path to the reference annotation GTF file (-a).
        """
        with open(out_file, 'w') as outfile:
            for sample in self.samples:
                bam = f'{sample}_Aligned.sortedByCoord.out.bam'
                if self.paired_end:
                    cmd = (
                        f'featureCounts -p -T {n_threads} -s {strand} -Q {min_quality} '
                        f'-a {gtf_file} -o {sample}_geneCounts.txt {bam}'
                    )
                else:
                    cmd = (
                        f'featureCounts -T {n_threads} -s {strand} -Q {min_quality} '
                        f'-a {gtf_file} -o {sample}_geneCounts.txt {bam}'
                    )
                outfile.write(cmd + '\n')

    def counting_exons(
        self,
        out_file: str = 'counting_exons.sh',
        gtf_file: str = 'GRCh38.115.gtf',
        strand: int = 0,
        min_quality: int = 30,
        n_threads: int = 4,
    ) -> None:
        """Write a shell script that runs featureCounts at the exon level.

        Iterates over ``self.samples`` (set by :meth:`ranseq_mapping`) and
        builds one featureCounts command per sample, using ``-J -f -t exon
        -O`` to count reads and junctions per exon feature (used later by
        :meth:`get_exonIRER`).

        Args:
            out_file: Path of the shell script to write.
            gtf_file: Path to the reference annotation GTF file (-a).
            strand: Strandedness passed to featureCounts (-s).
            min_quality: Minimum mapping quality passed to featureCounts (-Q).
            n_threads: Number of threads passed to featureCounts (-T).
        """
        with open(out_file, 'w') as outfile:
            for sample in self.samples:
                bam = f'{sample}_Aligned.sortedByCoord.out.bam'
                if self.paired_end:
                    cmd = (
                        f'featureCounts -J -f -t exon -O -p -T {n_threads} -s {strand} '
                        f'-Q {min_quality} -a {gtf_file} -o {sample}_exonCounts.txt {bam}'
                    )
                else:
                    cmd = (
                        f'featureCounts -J -f -t exon -O -T {n_threads} -s {strand} '
                        f'-Q {min_quality} -a {gtf_file} -o {sample}_exonCounts.txt {bam}'
                    )

    def get_exonIRER(
        self,
        exon_counts_dir: str = 'IRER',
        exon_table: str = 'Homo_sapiens.GRCh38.115.Exons',
        chrom: str = 'chr',
        suffix: str = '.bam',
        knownJuncsOnly: bool = True,
    ) -> None:
        """Compute per-exon intron-retention (IR) and exon-ratio (ER) counts.

        For every ``*_exonCounts.txt`` file in ``exon_counts_dir`` (produced
        by :meth:`counting_exons`), finds the matching junction-counts table
        (same path with ``_exonCounts.txt`` replaced by ``_jcounts.txt``)
        and combines them with the exon annotation table (``exon_table``,
        tab-separated with columns ``[ExonID, Chr, Start, End, GeneID,
        GeneName, ...]``) to produce a per-exon IR/ER table.

        IR is the per-exon, per-sample featureCounts read count (reads
        overlapping/included in the exon interval), taken directly from the
        exon counts table. ER is the per-exon, per-sample count of junction
        reads that splice across (skip) the exon -- i.e. junctions whose
        donor position is upstream of the exon start and whose acceptor
        position is downstream of the exon end, on the same chromosome and
        associated with the exon's gene. Junction records are read from the
        jcounts table, where columns 0/1 hold comma-separated gene IDs for
        the two splice sites, columns 2-4 describe the first splice site
        (chromosome, position, strand), columns 5-7 describe the second
        splice site, and columns 8 onward hold per-sample junction read
        counts. If ``knownJuncsOnly`` is True, only junctions whose two
        splice-site positions both match a known exon boundary (start or
        end position) found in ``exon_table`` are kept; otherwise all
        junctions are used. IR and ER are later combined into a
        percent-spliced-in (PSI) value by :meth:`get_exonPSI`.

        For each input ``*_exonCounts.txt`` file, writes a companion
        ``*_exonIRER.txt`` file with columns ``ExonID, Chr, Start, End,
        GeneID, GeneName`` followed by ``<sample>_IR`` and ``<sample>_ER``
        columns for every sample; exons missing an IR and/or ER value are
        written with ``-1`` placeholders.

        Args:
            exon_counts_dir: Directory containing ``*_exonCounts.txt`` and
                matching ``*_jcounts.txt`` tables produced by
                :meth:`counting_exons`.
            exon_table: Path to the tab-separated exon annotation table
                with columns ``[ExonID, Chr, Start, End, GeneID,
                GeneName, ...]``.
            chrom: Chromosome-name prefix used when building junction
                position keys and when comparing junction/exon
                chromosomes; also used as a truthiness flag to decide
                whether to prepend ``chr`` to junction chromosome names.
            suffix: Suffix stripped from BAM-derived column names (e.g.
                ``.bam``) to recover sample names from the counts/jcounts
                table headers.
            knownJuncsOnly: If True, restrict junction counts to those
                whose splice-site positions match known exon boundaries in
                ``exon_table``.

        Raises:
            ValueError: If the sample order/names inferred from the
                jcounts table header do not match those inferred from the
                exon counts table header.
        """
        exonCountsTables = [
            f'{exon_counts_dir}/{x}' for x in os.listdir(exon_counts_dir) if x.endswith('_exonCounts.txt')
        ]
        for exonCountsTable in sorted(exonCountsTables):
            jcountsTable = exonCountsTable.replace('_exonCounts.txt', '_jcounts.txt')
            print(f'processing {exonCountsTable} and {jcountsTable}', flush=True)
            IR: dict[str, list[str]] = {}
            ER: dict[str, Any] = {}
            JC: dict[str, list[list[str]]] = {}
            EL: list[list[str]] = []
            knownJuncs: dict[str, bool] = {}
            if knownJuncsOnly:
                with open(exon_table) as f:
                    for line in f:
                        line = line.strip()
                        fields = line.split('\t')
                        k1 = fields[1] + ':' + fields[2]
                        k2 = fields[1] + ':' + fields[3]
                        knownJuncs[k1] = True
                        knownJuncs[k2] = True

            with open(jcountsTable) as f:
                head = f.readline().strip().split('\t')
                sampleExon = [x.split(suffix)[0] for x in head[8:]]
                for line in f:
                    line = line.strip()
                    fields = line.split('\t')
                    if knownJuncsOnly:
                        k1 = chrom + fields[2] + ':' + fields[3]
                        k2 = chrom + fields[2] + ':' + fields[6]
                        if k1 in knownJuncs and k2 in knownJuncs:
                            if fields[0] != 'NA':
                                for k in fields[0].split(','):
                                    JC.setdefault(k, [])
                                    JC[k].append(fields)
                            if fields[1] != 'NA':
                                for k in fields[1].split(','):
                                    JC.setdefault(k, [])
                                    JC[k].append(fields)
                    else:
                        if fields[0] != 'NA':
                            for k in fields[0].split(','):
                                JC.setdefault(k, [])
                                JC[k].append(fields)
                        if fields[1] != 'NA':
                            for k in fields[1].split(','):
                                JC.setdefault(k, [])
                                JC[k].append(fields)

            with open(exonCountsTable) as f:
                for line in f:
                    line = line.strip()
                    fields = line.split('\t')
                    if line[0] != '#':
                        if fields[0] == 'Geneid':
                            sampleExon2 = [x.split(suffix)[0] for x in fields[6:]]
                        else:
                            k = '\t'.join(['chr' + fields[1], fields[2], fields[3]])
                            IR[k] = fields[6:]

            if sampleExon != sampleExon2:
                raise ValueError('samples inconsistent!')

            with open(exon_table) as f:
                for line in f:
                    line = line.strip()
                    fields = line.split('\t')
                    EL.append(fields)
                    exonID = fields[0]
                    geneID = fields[4]
                    ch = fields[1]
                    start = int(fields[2])
                    end = int(fields[3])
                    JCL = np.array([0] * len(sampleExon))
                    if geneID in JC:
                        for item in JC[geneID]:
                            s1_ch = item[2]
                            s2_ch = item[5]
                            if chrom:
                                s1_ch = 'chr' + item[2]
                                s2_ch = 'chr' + item[5]
                            s1_pos = int(item[3])
                            s1_strand = item[4]
                            s2_pos = int(item[6])
                            s2_strand = item[7]
                            SampleCounts = np.array([int(x) for x in item[8:]])
                            if s1_ch == ch and s2_ch == ch:
                                if s1_pos < start and s2_pos > end:
                                    JCL += SampleCounts
                    ER[exonID] = JCL

            out_file = exonCountsTable.replace('_exonCounts.txt', '_exonIRER.txt')
            with open(out_file, 'w') as f:
                f.write(
                    'ExonID\tChr\tStart\tEnd\tGeneID\tGeneName'
                    + '\t'
                    + '\t'.join(['%s_IR\t%s_ER' % (x, x) for x in sampleExon])
                    + '\n'
                )
                for E in EL:
                    k = '\t'.join([E[1], E[2], E[3]])
                    if k in IR and E[0] in ER:
                        f.write(
                            '\t'.join(E)
                            + '\t'
                            + '\t'.join(
                                [str(IR[k][n] + '\t' + str(ER[E[0]][n])) for n in range(0, len(IR[k]))]
                            )
                            + '\n'
                        )
                    elif k in IR:
                        f.write(
                            '\t'.join(E)
                            + '\t'
                            + '\t'.join([str(IR[k][n] + '\t' + '-1') for n in range(0, len(IR[k]))])
                            + '\n'
                        )
                    elif E[0] in ER:
                        f.write(
                            '\t'.join(E)
                            + '\t'
                            + '\t'.join(['-1' + '\t' + str(ER[E[0]][n]) for n in range(0, len(ER[E[0]]))])
                            + '\n'
                        )
                    else:
                        f.write(
                            '\t'.join(E) + '\t' + '\t'.join(['-1' + '\t' + '-1' for x in sampleExon]) + '\n'
                        )

    def merge_counts_tables(
        self,
        counts_tables: str = 'counts_tables.txt',
        out_file: str = 'eQTL_geneCounts.txt',
    ) -> None:
        """Merge per-sample gene/exon/exonIRER count tables into one table.

        Reads a list of per-sample table paths from ``counts_tables`` (one
        path per line, no header) and merges them column-wise into a single
        table. The type of each table (``gene``, ``exon``, or ``exonIRER``)
        is inferred from its file name, and the sample name is derived from
        the portion of the file name preceding the type-specific suffix
        (e.g. ``_geneCounts``). For ``gene`` and ``exon`` tables, only the
        featureCounts count column is kept per sample (plus the shared
        feature/annotation columns from the first table). For ``exonIRER``
        tables, both the ``_IR`` and ``_ER`` columns are kept per sample.
        Exon tables have duplicate rows (by chromosome/start/end/strand)
        dropped before writing.

        Args:
            counts_tables: Path to a text file listing the per-sample
                counts table paths to merge, one per line.
            out_file: Path of the merged, tab-separated output table.

        Raises:
            ValueError: If a listed table's file name does not match the
                expected ``geneCounts``/``exonCounts``/``exonIRER`` naming,
                or if the merged tables do not all have the same number of
                feature rows.
        """
        df_tables = pd.read_table(counts_tables, header=None, low_memory=False)
        L = []
        n_features = []
        for n in range(df_tables.shape[0]):
            f = df_tables.iloc[n, 0]
            if f.split('/')[-1].find('geneCounts') != -1:
                counts_type = 'gene'
            elif f.split('/')[-1].find('exonCounts') != -1:
                counts_type = 'exon'
            elif f.split('/')[-1].find('exonIRER') != -1:
                counts_type = 'exonIRER'
            else:
                raise ValueError('check type of the counts tables')

            sample = f.split('/')[-1].split('_' + counts_type)[0]
            if counts_type == 'gene':
                df2 = pd.read_table(f, header=0, comment='#', low_memory=False)
                n_features.append(df2.shape[0])
                if n == 0:
                    df3 = df2.iloc[:, [0, 6]]
                    df3.columns = ['GeneID', sample]
                    L.append(df3)
                else:
                    df3 = df2.iloc[:, 6]
                    df3.name = sample
                    L.append(df3)
            elif counts_type == 'exon':
                df2 = pd.read_table(f, header=0, comment='#', low_memory=False)
                n_features.append(df2.shape[0])
                if n == 0:
                    df3 = df2.iloc[:, 0:7]
                    df3.columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand', 'Length', sample]
                    L.append(df3)
                else:
                    df3 = df2.iloc[:, 6]
                    df3.name = sample
                    L.append(df3)
            elif counts_type == 'exonIRER':
                df2 = pd.read_table(f, header=0, comment='#', low_memory=False)
                n_features.append(df2.shape[0])
                if n == 0:
                    df3 = df2.iloc[:, 0:8]
                    df3.columns = [
                        'ExonID', 'Chr', 'Start', 'End', 'GeneID', 'GeneName', f'{sample}_IR', f'{sample}_ER',
                    ]
                    L.append(df3)
                else:
                    df3 = df2.iloc[:, 6:8]
                    df3.columns = [f'{sample}_IR', f'{sample}_ER']
                    L.append(df3)
        df4 = pd.concat(L, axis=1)
        if np.sum(np.array(n_features) != n_features[0]) > 0:
            raise ValueError('number of features in the counts tables are different')
        if counts_type == 'exon':
            df4.drop_duplicates(subset=['Chr', 'Start', 'End', 'Strand'], inplace=True)
        df4.to_csv(out_file, sep='\t', index=False)

    def get_exonPSI(self, in_file: str = 'sQTL_exonIRER.txt') -> None:
        """Compute percent-spliced-in (PSI) values from an exonIRER table.

        Reads an ``*_exonIRER.txt`` table (as produced by
        :meth:`get_exonIRER`, with per-sample ``_IR``/``_ER`` column pairs)
        and computes, for each exon and sample, ``PSI = IR / (IR + ER)``
        (using ``PSI = (IR + 1) / (IR + ER + 1)`` when ``IR + ER == 0`` to
        avoid division by zero). Writes two output files: one with all
        original columns plus a ``_PSI`` column per sample (path derived by
        appending ``_PSI`` before the ``.txt`` extension of ``in_file``),
        and a compact one with only ``ExonID``, ``GeneName``, and the
        per-sample PSI values (path derived by replacing ``exonIRER.txt``
        with ``exonPSI.txt`` in ``in_file``).

        Args:
            in_file: Path to the input ``*_exonIRER.txt`` table.
        """
        out_file = in_file.split('.txt')[0] + '_PSI.txt'
        out_file2 = in_file.replace('exonIRER.txt', 'exonPSI.txt')
        with open(in_file) as f, open(out_file, 'w') as fout, open(out_file2, 'w') as fout2:
            head = f.readline().strip().split('\t')
            H = head
            H2 = ['ExonID', 'GeneName']
            for n in range(6, len(head), 2):
                H.append(head[n].split('_IR')[0] + '_PSI')
                H2.append(head[n].split('_IR')[0])
            fout.write('\t'.join(H) + '\n')
            fout2.write('\t'.join(H2) + '\n')

            for line in f:
                line = line.strip()
                fields = line.split('\t')
                L = fields
                L2 = ['_'.join(fields[1:4] + fields[0:1] + fields[4:5]), fields[5]]
                for n in range(6, len(fields), 2):
                    IR = float(fields[n])
                    ER = float(fields[n + 1])
                    if IR + ER == 0:
                        PSI = (IR + 1) / (IR + ER + 1)
                    else:
                        PSI = IR / (IR + ER)
                    L.append(f'{PSI:.4f}')
                    L2.append(f'{PSI:.4f}')
                fout.write('\t'.join(L) + '\n')
                fout2.write('\t'.join(L2) + '\n')

    def annotate_gene_name(
        self,
        in_file: str = 'eQTL_geneCounts.txt',
        gene_table: str = 'GRCh38.115_GenePosType.txt',
        sep: str = '\t',
    ) -> None:
        """Add a GeneName column to a gene counts table.

        Looks up gene names from ``gene_table`` (a two-column, headerless
        table mapping gene ID to gene name, e.g. as produced by
        ``gtf_to_GenePosType`` in ``utils.py``) by matching the
        version-stripped Ensembl gene ID (everything before the first
        ``.``) in the first column of ``in_file``. Genes without a match
        keep their (version-stripped) gene ID as the gene name. Writes the
        result to a new file with ``_geneName.txt`` appended to the base
        name of ``in_file``.

        Args:
            in_file: Path to the input counts table; its first column must
                hold gene IDs.
            gene_table: Path to a headerless, tab-separated table with gene
                ID in column 0 and gene name in column 1.
            sep: Field separator used to parse ``in_file``.

        Raises:
            ValueError: If ``gene_table`` does not exist.
        """
        if not os.path.exists(gene_table):
            raise ValueError(f'{gene_table} is not found, run gtf_to_GenePosType in utils.py on the gtf file first')

        out_file = in_file.split('.txt')[0].split('.csv')[0] + '_geneName.txt'
        df = pd.read_table(gene_table, header=None)
        D = dict(zip(df[0], df[1]))
        with open(in_file, 'r') as f, open(out_file, 'w') as fout:
            head = f.readline().strip().split(sep)
            fout.write('\t'.join(head[0:1] + ['GeneName'] + head[1:]) + '\n')
            for line in f:
                line = line.strip()
                fields = line.split(sep)
                gene_id = fields[0].split('.')[0]
                gene_name = D.get(gene_id, gene_id)
                fout.write('\t'.join([gene_id, gene_name] + fields[1:]) + '\n')

    def annotate_exon_name(
        self,
        exon_table: str = 'GRCh38.115_Exons.txt',
        in_file: str = 'eQTL_exonCounts.txt',
    ) -> None:
        """Add ExonID/GeneName columns to a merged exon counts table.

        Builds a lookup from ``exon_table`` (tab-separated, columns
        ``[ExonID, Chr, Start, End, GeneID, GeneName, ...]``) keyed by
        ``Chr_Start_End`` to recover the exon's original ID and gene name,
        and a second lookup from gene ID to gene name. For each row of
        ``in_file`` (a merged exon counts table, gene ID in column 0 and
        chromosome/start/end in columns 1-3), matches on
        ``chr<Chr>_<Start>_<End>`` and writes ``<ExonID>, GeneName`` plus
        the remaining count columns (from column 6 onward); rows without a
        match are written with an ExonID of ``<chr_start_end>_NA_<geneID>``
        and the gene name (or gene ID, if unknown) instead. Writes the
        result to a new file with ``.txt`` replaced by ``_geneName.txt``.

        Args:
            exon_table: Path to the tab-separated exon annotation table
                with columns ``[ExonID, Chr, Start, End, GeneID,
                GeneName, ...]``.
            in_file: Path to the merged exon counts table to annotate.

        Raises:
            ValueError: If ``exon_table`` does not exist.
        """
        if not os.path.exists(exon_table):
            raise ValueError(f'{exon_table} is not found')

        out_file = in_file.replace('.txt', '_geneName.txt')
        D: dict[str, list[str]] = {}
        G: dict[str, str] = {}
        with open(exon_table) as f:
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                exonID = '_'.join(fields[1:4])
                geneID = fields[4]
                geneName = fields[5]
                if exonID not in D:
                    D[exonID] = fields[1:4] + fields[0:1] + fields[4:6]
                if geneID not in G:
                    G[geneID] = geneName

        with open(in_file, 'r') as f, open(out_file, 'w') as fout:
            head = f.readline().strip().split('\t')
            fout.write('ExonID\tGeneName\t' + '\t'.join(head[6:]) + '\n')
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                exonID = 'chr' + '_'.join(fields[1:4])
                geneID = fields[0]
                geneName = G.get(geneID, geneID)
                if exonID in D:
                    fout.write('\t'.join(['_'.join(D[exonID][0:-1]), D[exonID][-1]] + fields[6:]) + '\n')
                else:
                    fout.write('\t'.join([exonID + '_NA_' + geneID, geneName] + fields[6:]) + '\n')

    def counts_to_tpm(
        self,
        counts_table: str = 'eQTL_geneCounts_geneName.txt',
        counts_sample: str = 'sample_geneCounts.txt',
        sample_start_idx: int = 2,
        norm_base: float = 1e6,
    ) -> None:
        """Normalize a counts table to TPM.

        Equivalent to the edgeR convention
        ``TPM <- t(t(RPKM) / colSums(RPKM)) * 1e6``: counts are first
        divided by feature length, then each sample column is rescaled so
        its length-normalized counts sum to ``norm_base``. Feature length
        is taken from the ``Length`` column of ``counts_sample`` (a raw
        featureCounts table matching ``counts_table`` row-for-row) if it
        exists and matches; otherwise, if ``counts_table`` has an
        ``ExonID`` column formatted as ``..._<start>_<end>_...``, length is
        derived as ``end - start + 1``; otherwise every feature is treated
        as length 1 (i.e. no length normalization, only library-size
        normalization).

        Args:
            counts_table: Path to the input counts table; columns before
                ``sample_start_idx`` are treated as feature/annotation
                columns and columns from ``sample_start_idx`` onward are
                treated as per-sample counts.
            counts_sample: Path to a raw featureCounts table (with a
                ``Length`` column) matching ``counts_table`` row order,
                used to recover feature lengths. Ignored if it does not
                exist.
            sample_start_idx: Index of the first sample count column in
                ``counts_table``.
            norm_base: Target sum for length-normalized counts per sample
                (1e6 for TPM).

        Raises:
            ValueError: If ``counts_table`` and ``counts_sample`` have the
                same number of rows but their feature identifiers (first
                column) do not match.
        """
        Length: Any = []
        df1 = pd.read_table(counts_table, header=0)
        if counts_sample and os.path.exists(counts_sample):
            df2 = pd.read_table(counts_sample, header=0, comment='#')
            if df1.shape[0] == df2.shape[0]:
                wh = df1.iloc[:, 0] == df2.iloc[:, 0]
                if wh.all():
                    if 'Length' in df2.columns:
                        Length = df2['Length']
                else:
                    raise ValueError('Feature and Length are not the same version')
        elif 'ExonID' in df1.columns:
            Length = [int(x.split('_')[2]) - int(x.split('_')[1]) + 1 for x in df1['ExonID']]
            print('Normalized by the Exon Length calculated from the ExonID column')
        else:
            Length = [1] * df1.shape[0]
            print('Warning: Not normalized by Length!')

        out_file = counts_table.split('.txt')[0] + '_TPM.txt'

        mat = df1.iloc[:, sample_start_idx:]
        mat2 = df1.iloc[:, 0:sample_start_idx]

        matTotalRaw = mat.sum(axis=0)
        print(f'Total Reads (million):\n{matTotalRaw / norm_base}')
        matLength = (mat.T / Length).T
        matTotal = matLength.sum(axis=0)
        M = matLength / matTotal * norm_base
        df = pd.concat([mat2, M], axis=1)
        df.to_csv(out_file, header=True, index=False, sep='\t', float_format='%.4f')

    def filter_PSI_by_variance(
        self,
        in_file: str = 'sQTL_exonPSI.txt',
        var_threshold: float = 0.001,
    ) -> None:
        """Drop low-variance exons from a PSI table.

        Reads an ``*_exonPSI.txt`` table (as produced by
        :meth:`get_exonPSI`, feature columns followed by per-sample PSI
        values starting at column index 2), computes the variance of each
        row's PSI values across samples, and keeps only rows whose variance
        exceeds ``var_threshold``. Writes the filtered table to a new file
        with ``.txt`` replaced by ``_exonFiltered.txt``.

        Args:
            in_file: Path to the input PSI table.
            var_threshold: Minimum per-exon PSI variance required to keep
                a row.

        Raises:
            FileNotFoundError: If ``in_file`` does not exist.
        """
        if os.path.exists(in_file):
            df = pd.read_table(in_file, header=0, sep='\t')
        else:
            raise FileNotFoundError(f'{in_file} not found.')
        wh = []
        for n in range(df.shape[0]):
            var = np.var(df.iloc[n, 2:].values.astype(float))
            if var > var_threshold:
                wh.append(True)
            else:
                wh.append(False)
        df_filtered = df.loc[wh, ]
        out_file = in_file.replace('.txt', f'_exonFiltered.txt')
        df_filtered.to_csv(out_file, header=True, index=False, sep='\t')
