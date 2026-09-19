"""Summary reporting utilities for the omniQTL pipeline.

This module provides the ``Summary`` class, a collection of methods for
generating QQ plots, significant-locus/independent-signal bar plots, and
upset plots of shared genes/donors across caQTL/eQTL/pQTL analyses. It also
builds donor QC/metadata summary tables (merging REDCap records, genetic
ancestry/sex calls, and RRID lookups), summarizes peak count and length
distributions, computes variant functional-consequence and regulatory-region
enrichment via Fisher's exact test, and compares overlap/correlation of
significant associations against external datasets such as GTEx, UK Biobank
plasma pQTL, InsPIRE, and PRS variant lists.
"""
from typing import Any

from .utils import *


class Summary:
    """Collection of QTL summary, QC, and enrichment reporting methods.

    This class bundles a large set of loosely related, mostly
    file-in/file-out helper methods used to build summary tables and
    publication-style plots (QQ plots, bar plots, heatmaps, upset plots,
    and correlation scatter plots) for caQTL/eQTL/pQTL mapping results,
    donor QC/metadata, and cross-study overlap/enrichment analyses. Most
    methods read one or more delimited text files, compute a derived
    table or figure, and write the result back to disk rather than
    returning a value.
    """

    def __init__(self) -> None:
        """Initialize the Summary object (no state is stored)."""
        pass

    def get_table_for_qq_plot(
        self,
        in_file: str = 'pQTL_permute-1000_w1M_PC25_extraInfo.txt.gz',
        out_file: str = 'pQTL_qq_plot_table.txt',
        p_col: str = 'adj_beta_pval',
    ) -> None:
        """Build a table of expected vs. observed -log10(p) values for a QQ plot.

        For each phenotype (``phe_id``) group in ``in_file``, the minimum
        p-value in ``p_col`` is taken, sorted ascending, and paired with the
        corresponding expected p-value under the uniform null to produce
        ``expected``/``observed`` -log10 columns.

        Args:
            in_file: Path to the input tab-delimited association file.
            out_file: Path to write the resulting QQ-plot table.
            p_col: Name of the p-value column to summarize per phenotype.
        """
        L = []
        df = pd.read_table(in_file, header=0, sep='\t')
        for gi, g in df.groupby('phe_id'):
            g_sub = g.dropna(subset=[p_col])
            if g_sub.shape[0]:
                L.append([gi, g_sub[p_col].min()])
        df = pd.DataFrame(L, columns=['phe_id', 'min_p'])
        df.sort_values('min_p', inplace=True)
        n = df.shape[0]
        df['expected'] = -np.log10(np.arange(1, n + 1) / (n + 1))
        df['observed'] = -np.log10(df['min_p'])
        df.to_csv(out_file, index=False, sep='\t')

    def qq_plot(
        self,
        in_file: str,
        title: str = 'QQ plot',
        scatter_size: float = 4,
        color: str = 'C1',
        figsize: tuple[float, float] = (4, 4),
    ) -> None:
        """Draw a QQ plot from a table produced by ``get_table_for_qq_plot``.

        Args:
            in_file: Path to the tab-delimited table with ``expected`` and
                ``observed`` columns; the output PDF path is derived by
                replacing ``.txt`` with ``.pdf``.
            title: Plot title.
            scatter_size: Marker size for the scatter points.
            color: Color for the scatter points.
            figsize: Figure size in inches, as ``(width, height)``.
        """
        out_file = in_file.replace('.txt', '.pdf')
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.scatterplot(x='expected', y='observed', s=scatter_size, ax=ax, data=df, color=color)
        ax.plot([0, df['expected'].max()], [0, df['expected'].max()], linestyle="--", color='C0')
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(out_file)

    def bar_plot_significant_loci(
        self,
        in_file: str,
        axes: list[float] = [0.3, 0.4, 0.6, 0.5],
        cmap: str = 'Dark2',
        show_numbers: bool = True,
        figsize: tuple[float, float] = (4, 4),
        ylabel: str = 'Number of significant signals',
    ) -> None:
        """Plot a bar chart of the number of significant loci per study.

        Draws bars for each study/sample-size combination (grouped visually
        into caQTL, eQTL, and pQTL brackets via manual bracket lines), with
        optional value labels on top of each bar.

        Args:
            in_file: Path to a tab-delimited file with ``Study``,
                ``SampleSize``, and ``Number of significant loci`` columns;
                the output PDF path replaces the first ``.txt`` segment
                with ``.pdf``.
            axes: Axes rectangle ``[left, bottom, width, height]`` in figure
                coordinates used for ``fig.add_axes``.
            cmap: Seaborn/matplotlib color palette name.
            show_numbers: Whether to annotate each bar with its value.
            figsize: Figure size in inches, as ``(width, height)``.
            ylabel: Y-axis label; if falsy, no label is set.
        """
        out_file = in_file.split('.txt')[0] + '.pdf'
        cmap = sns.color_palette(cmap)

        df = pd.read_table(in_file, header=0, sep='\t')
        df['StudySampleSize'] = [
            f"{df['Study'].iloc[n]}\n(N={df['SampleSize'].iloc[n]})" for n in range(df.shape[0])
        ]

        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes(axes)
        sns.barplot(
            y='Number of significant loci',
            x='StudySampleSize',
            hue='StudySampleSize',
            data=df,
            palette=[cmap[0], cmap[-1], cmap[1], cmap[-1], cmap[2]],
            legend=False,
            linewidth=1,
            edgecolor='black',
        )
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_xlabel('')
        if ylabel:
            ax.set_ylabel(ylabel)
        if show_numbers:
            for i, row in df.iterrows():
                N = row['Number of significant loci']
                ax.text(i, N, N, ha='center', va='bottom')
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], ylim[1] * 1.1)

        y1 = -0.6
        y2 = -0.65
        ax.plot([0, 0, 1, 1], [y1, y2, y2, y1], transform=ax.get_xaxis_transform(), lw=1, color='k', clip_on=False)
        ax.plot([2, 2, 3, 3], [y1, y2, y2, y1], transform=ax.get_xaxis_transform(), lw=1, color='k', clip_on=False)
        ax.plot([4, 4], [y1, y2], transform=ax.get_xaxis_transform(), lw=1, color='k', clip_on=False)
        ax.text(0.5, y2 - 0.03, 'caQTL', ha='center', va='top', transform=ax.get_xaxis_transform())
        ax.text(2.5, y2 - 0.03, 'eQTL', ha='center', va='top', transform=ax.get_xaxis_transform())
        ax.text(4, y2 - 0.03, 'pQTL', ha='center', va='top', transform=ax.get_xaxis_transform())

        #plt.tight_layout()
        plt.savefig(out_file)

    def get_table_for_upset_plot(
        self,
        in_files: list[str] = ['caQTL_permute-1000_w1k_qvalue.significant.txt'],
        out_file: str = 'QTL_upset_plot_table.txt',
    ) -> None:
        """Summarize the unique gene sets targeted by each QTL type for an upset plot.

        For each input file, the QTL type is inferred from the filename
        prefix. For ``caQTL`` files, comma-separated gene lists embedded in
        the first column are split and pooled; for other QTL types, the
        last underscore-delimited token of the first column is treated as
        the gene id.

        Args:
            in_files: List of paths to significant-association files, one
                per QTL type, whose filenames start with the QTL type
                (e.g. ``caQTL_...``).
            out_file: Path to write the resulting summary table.
        """
        L = []
        for f in in_files:
            qtl = f.split('_')[0]
            genes = []
            df = pd.read_table(f, header=None, sep=r'\s+')
            if qtl in ['caQTL']:
                for item in df.iloc[:, 0]:
                    x = item.split('_')[-1].split(',')
                    genes += x
            else:
                for item in df.iloc[:, 0]:
                    genes.append(item.split('_')[-1])
            genes_uniq = sorted(set(genes))
            L.append([qtl, len(genes_uniq), ','.join(genes_uniq)])
        df = pd.DataFrame(L, columns=['qtl', 'numer_of_genes', 'genes'])
        df.to_csv(out_file, header=True, index=False, sep='\t')

    def plot_upset_qtl(
        self,
        in_file: str = 'QTL_upset_plot_table.txt',
        cmap: str = 'deep',
    ) -> None:
        """Draw an upset plot of shared genes across QTL types.

        Reads the table produced by ``get_table_for_upset_plot``, builds a
        membership dict of gene sets keyed by QTL type, and colors specific
        subsets (single QTL types, pairwise intersections, and the
        three-way intersection) with distinct palette colors.

        Args:
            in_file: Path to the tab-delimited upset-plot table; the output
                PDF path appends ``_upset.pdf`` after stripping ``.txt``.
            cmap: Seaborn color palette name used to style subsets.
        """
        from upsetplot import from_contents
        from upsetplot import UpSet
        cmap = sns.color_palette(cmap)
        out_file = in_file.split('.txt')[0] + '_upset.pdf'

        df = pd.read_table(in_file, header=0, sep='\t')
        D = {}
        for n in range(df.shape[0]):
            D[df.iloc[n, 0]] = df.iloc[n, -1].split(',')

        D2 = from_contents(D)
        ax = UpSet(D2, subset_size="count", facecolor='C0', show_counts=True)
        ax.style_subsets(present='caQTL', facecolor=cmap[0])
        ax.style_subsets(present='eQTL', facecolor=cmap[1])
        ax.style_subsets(present='pQTL', facecolor=cmap[2])
        ax.style_subsets(present=('caQTL', 'eQTL'), facecolor=cmap[4])
        ax.style_subsets(present=('caQTL', 'pQTL'), facecolor=cmap[4])
        ax.style_subsets(present=('eQTL', 'pQTL'), facecolor=cmap[4])
        ax.style_subsets(present=('caQTL', 'eQTL', 'pQTL'), facecolor=cmap[5])
        ax.plot()
        plt.savefig(out_file)

    def plot_upset_donors(
        self,
        in_files: list[str] = [],
        out_file: str = 'donors_upset.pdf',
        color: str = 'C0',
    ) -> None:
        """Draw an upset plot of shared donors across sample lists.

        Args:
            in_files: List of paths to single-column, headerless donor-id
                files; the group name for each is derived from the part of
                the filename before the first underscore.
            out_file: Path to write the resulting PDF.
            color: Facecolor used for all upset-plot subsets.
        """
        from upsetplot import from_contents
        from upsetplot import UpSet

        D = {}
        for f in in_files:
            df = pd.read_table(f, header=None, sep='\t')
            typ = f.split('_')[0]
            D[typ] = df.iloc[:, 0]

        D2 = from_contents(D)
        ax = UpSet(D2, subset_size="count", facecolor=color, show_counts=True)
        ax.plot()
        plt.savefig(out_file)

    def get_summary_table_donor_qc(
        self,
        in_files: list[str] = [
            'caQTL_samples.txt', 'eQTL_samples.txt', 'pQTL_samples.txt', 'GSIS_samples.txt',
        ],
        in_files2: list[str] = [
            'ATACseq_number_mapped_reads.txt', 'RNAseq_number_mapped_reads.txt',
            'ATACseq_tss_score.txt', 'RNAseq_data_source.txt',
        ],
        in_files3: list[str] = ['ATACseq_number_peaks_by_qvalue.txt'],
        out_file: str = 'donor_qc_summary.txt',
        cols: list[str] = [
            'donor_id', 'caQTL', 'eQTL', 'pQTL', 'GSIS', 'ATACseq_number_mapped_reads',
            'ATACseq_tss_score', 'ATACseq_number_peaks_by_qvalue', 'RNAseq_number_mapped_reads',
            'RNAseq_data_source',
        ],
    ) -> None:
        """Build a combined donor QC/inclusion summary table across QTL modalities.

        Merges donor membership lists (``in_files``, e.g. per-QTL-type
        sample lists marked ``Yes``), per-donor QC metrics (``in_files2``),
        and per-donor peak counts (``in_files3``) into one wide table keyed
        by donor id, sorted by inclusion across the four QTL/GSIS modalities
        then by donor id. Also writes a filtered subset of new,
        non-duplicated samples that are included in at least one modality.

        Args:
            in_files: Paths to headerless per-modality donor-id list files
                (e.g. ``caQTL_samples.txt``); presence marks ``Yes``.
            in_files2: Paths to headered per-donor QC metric files whose
                first two columns are donor id and value.
            in_files3: Paths to headered per-donor peak-count files whose
                second and third columns are donor id and value.
            out_file: Path to write the combined summary table; a second
                file with suffix ``_new_samples.txt`` is also written.
            cols: Output column names, with ``cols[0]`` as the donor-id
                column and the remainder matching the modality/metric keys
                derived from ``in_files``/``in_files2``/``in_files3``.
        """
        D = {}
        S = []
        for f in in_files:
            k = f.split('_')[0]
            df = pd.read_table(f, header=None, sep='\t')
            D[k] = {s: 'Yes' for s in df.iloc[:, 0]}
            S += D[k]
        for f in in_files2:
            k = f.split('.txt')[0]
            df = pd.read_table(f, header=0, sep='\t')
            D[k] = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
            S += df.iloc[:, 0].tolist()
        for f in in_files3:
            k = f.split('.txt')[0]
            df = pd.read_table(f, header=0, sep='\t')
            D[k] = dict(zip(df.iloc[:, 1], df.iloc[:, 2]))
            S += df.iloc[:, 1].tolist()
        S = sorted(set(S))

        Ls = []
        for sample in S:
            L = []
            for col in cols[1:]:
                val = 'N/A'
                if col in D:
                    if sample in D[col]:
                        val = D[col][sample]
                L.append(val)
            Ls.append([sample] + L)
        df = pd.DataFrame(Ls, columns=cols)
        df['sort'] = df.iloc[:, 1:5].apply(lambda x: ''.join(['1' if i == 'Yes' else '0' for i in x]), axis=1)
        df.sort_values(['sort', 'donor_id'], ascending=[False, True], inplace=True)
        df.drop(columns=['sort'], inplace=True)
        df.to_csv(out_file, index=False, sep='\t')

        wh1 = (df['eQTL'] == 'Yes') & df['RNAseq_data_source'].str.contains('new')
        wh2 = df['caQTL'] == 'Yes'
        wh3 = df['pQTL'] == 'Yes'
        wh4 = df['GSIS'] == 'Yes'
        wh5 = df['donor_id'].str.contains('ISLET')
        print(f'eQTL: {sum(wh1)}')
        print(f'caQTL: {sum(wh2)}')
        print(f'pQTL: {sum(wh3)}')
        print(f'GSIS: {sum(wh4)}')
        print(f'duplicated: {sum(wh5)}')
        df_new = df[(wh1 | wh2 | wh3 | wh4) & (~wh5)].iloc[:, 0:5]
        df_new.to_csv(out_file.replace('.txt', '_new_samples.txt'), index=False, sep='\t')
        print(f'final number of unique samples: {df_new.shape[0]}')

    def get_donor_info(
        self,
        in_file: str,
        gap_file: str = 'GAP_Oxford-Stanford_conditional_batch18_extraInfo.txt',
        human_islets_file: str = 'human_islets.txt',
        redcap_file: str = 'REDCap.csv',
        extra_file1: str = '',
        extra_file2: str = '',
    ) -> None:
        """Assemble a per-donor metadata table from multiple source files.

        Starting from ``in_file`` (per-donor QTL/genotyping inclusion
        flags), this adds: a ``gwas`` column copied from ``genotyped``;
        national-id and RRID cross-references parsed from ``redcap_file``;
        sex/age/BMI parsed from an optional ``extra_file1`` REDCap-like
        extract, keyed by national id; genetic sex/ancestry/source parsed
        from ``gap_file`` and merged in by donor id (with a warning printed
        if a genetic-sex value contains multiple space-separated tokens);
        a recomputed ``genotyped`` flag that is true if any of
        gwas/rna/atac/protein/genetic-sex is present; and RRID/sex/age/
        BMI/HbA1c pulled from ``human_islets_file`` by donor id, backfilled
        from the REDCap/extra-file cross-references (matched via national
        id) when missing, with a special-cased age of ``53`` hardcoded for
        donor ``H522`` when otherwise unknown.

        Args:
            in_file: Path to the base per-donor inclusion-flag table; the
                output path replaces the first ``.txt`` segment with
                ``_donor_info.txt``.
            gap_file: Path to the genetic ancestry/sex/source table (GAP
                pipeline output).
            human_islets_file: Path to the Human Islets REDCap-style export
                with RRID, sex, age, height, weight, BMI, and HbA1c.
            redcap_file: Path to a REDCap CSV export with national id and
                RRID per donor; ignored if empty or missing.
            extra_file1: Path to an optional supplementary tab-delimited
                file with ``National ID``, ``Gender``, ``Age``, and ``BMI``
                columns; ignored if empty or missing.
            extra_file2: Unused placeholder for an additional supplementary
                file path.
        """
        out_file = in_file.split('.txt')[0] + '_donor_info.txt'
        df = pd.read_table(in_file, header=0, sep='\t')

        # change genotyped column to gwas column, and update the value for subset of samples
        df['gwas'] = df['genotyped'].copy()

        # alternative id for some donors
        NID = {}
        RRID = {}
        SEX = {}
        AGE = {}
        BMI = {}
        if redcap_file and os.path.exists(redcap_file):
            df_redcap = pd.read_csv(redcap_file, header=0)
            for n in range(df_redcap.shape[0]):
                donor_id = df_redcap['sample_id'].iloc[n]
                national_id = str(df_redcap['national_id'].iloc[n])
                rrid = str(df_redcap['rrid'].iloc[n])
                if national_id != 'nan':
                    NID.setdefault(donor_id, set())
                    NID[donor_id].add(national_id)
                if rrid != 'nan':
                    RRID.setdefault(donor_id, set())
                    RRID[donor_id].add(rrid)
        df['national_id'] = [','.join(sorted(NID.get(x, ['.']))) for x in df['donor_id']]

        if extra_file1 and os.path.exists(extra_file1):
            df_extra1 = pd.read_table(extra_file1, header=0, sep='\t')
            for n in range(df_extra1.shape[0]):
                donor_id = df_extra1['National ID'].iloc[n]
                sex = str(df_extra1['Gender'].iloc[n]).lower()
                age = float(df_extra1['Age'].iloc[n])
                bmi = float(df_extra1['BMI'].iloc[n])
                if rrid != 'nan':
                    RRID.setdefault(donor_id, set())
                    RRID[donor_id].add(rrid)
                if sex != 'nan':
                    if sex.lower() in ['m']:
                        sex = 'male'
                    elif sex.lower() in ['f']:
                        sex = 'female'
                    SEX.setdefault(donor_id, set())
                    SEX[donor_id].add(sex)
                if str(age) != 'nan':
                    AGE.setdefault(donor_id, set())
                    AGE[donor_id].add(str(age))
                if str(bmi) != 'nan':
                    BMI.setdefault(donor_id, set())
                    BMI[donor_id].add(str(bmi))

        # genetic sex, ancestry, source
        df_gap = pd.read_table(gap_file, header=0, sep='\t')
        gap = {}
        for n in range(df_gap.shape[0]):
            sample = df_gap['SampleName'].iloc[n]
            sex = df_gap['Sex'].iloc[n]
            ancestry = df_gap['Superpopulation'].iloc[n]
            source = df_gap['Source'].iloc[n]
            gap[sample] = (sex, ancestry, source)

        E = {}
        E['sex'] = []
        E['ancestry'] = []
        E['source'] = []
        for n in range(df.shape[0]):
            sample = df['donor_id'].iloc[n]
            sex = '.'
            ancestry = '.'
            source = 'undetermined'
            if sample in gap:
                sex = gap[sample][0]
                if sex.find(' ') != -1:
                    sex = sex.split()[0]
                    print(f'warning: genetic sex of sample {sample} is {sex}')
                ancestry = gap[sample][1]
                source = gap[sample][2]
            E['sex'].append(sex)
            E['ancestry'].append(ancestry)
            E['source'].append(source)
        df['genetic_sex'] = E['sex']
        df['genetic_ancestry'] = E['ancestry']
        #df['center'] = E['source']

        # update the genotyped column
        whs = {}
        for col in df.columns[1:]:
            whs[col] = df[col] == 1
        whs['genetic_sex'] = df['genetic_sex'] != '.'
        whs['genetic_ancestry'] = df['genetic_ancestry'] != '.'

        wh = whs['gwas'] | whs['rna'] | whs['atac'] | whs['protein'] | whs['genetic_sex']
        df['genotyped'] = wh.astype(int)
        whs['genotyped'] = wh

        # rrid
        df_hi = pd.read_table(human_islets_file, header=0, sep='\t')
        hi = {}
        for n in range(df_hi.shape[0]):
            sample = df_hi['record_id'].iloc[n]
            rrid = df_hi['rrid'].iloc[n]
            age = df_hi['donorage'].iloc[n].astype(float)
            sex = df_hi['donorsex'].iloc[n]
            height = df_hi['donorheight'].iloc[n].astype(float)
            weight = df_hi['donorweight'].iloc[n].astype(float)
            bmi = df_hi['bodymassindex'].iloc[n].astype(float)
            hba1c = df_hi['hba1c'].iloc[n]
            hi[sample] = (rrid, sex, age, bmi, hba1c, height, weight)

        v = ['.'] * 7
        df['RRID'] = [hi.get(sample, v)[0] for sample in df['donor_id']]
        df['sex'] = [hi.get(sample, v)[1] for sample in df['donor_id']]
        df['age'] = [hi.get(sample, v)[2] for sample in df['donor_id']]
        df['bmi'] = [hi.get(sample, v)[3] for sample in df['donor_id']]
        df['hba1c'] = [hi.get(sample, v)[4] for sample in df['donor_id']]
        #df['height'] = [hi.get(sample, v)[5] for sample in df['donor_id']]
        #df['weight'] = [hi.get(sample, v)[6] for sample in df['donor_id']]
        print(df)

        for n in range(df.shape[0]):
            sample = df['national_id'].iloc[n]
            donor = df['donor_id'].iloc[n]
            rrid = df['RRID'].iloc[n]
            sex = df['sex'].iloc[n]
            age = df['age'].iloc[n]
            bmi = df['bmi'].iloc[n]
            if rrid == '.' and sample in RRID:
                df.at[n, 'RRID'] = ','.join(sorted(RRID[sample]))
            if sex == '.' and sample in SEX:
                df.at[n, 'sex'] = ','.join(sorted(SEX[sample]))
            if age == '.' and sample in AGE:
                df.at[n, 'age'] = ','.join(sorted(AGE[sample]))
            if bmi == '.' and sample in BMI:
                df.at[n, 'bmi'] = ','.join(sorted(BMI[sample]))
            if age == '.' and donor == 'H522':
                df.at[n, 'age'] = '53'
        df.to_csv(out_file, index=False, sep='\t')

    def summarize_donor_info(self, in_file: str = 'donor_info.txt') -> None:
        """Summarize donor counts per QC/inclusion criterion and flag sex mismatches.

        Counts total donors and, for each column, the number of donors
        with a positive flag (``== 1`` for the first several columns) or a
        non-missing value (``!= '.'`` otherwise). Also reports donors
        included in any QTL modality, any GWAS/QTL modality, and those
        with a mismatch between reported (``sex``) and genetically
        inferred (``genetic_sex``) sex (printed, not written to file).

        Args:
            in_file: Path to the donor-info table produced by
                ``get_donor_info``; the output path replaces ``.txt`` with
                ``_summary.txt``.
        """
        df = pd.read_table(in_file, header=0, sep='\t')
        out_file = in_file.replace('.txt', '_summary.txt')
        L = []
        L.append(['total_donors', df.shape[0], ''])
        for col in df.columns[1:]:
            if col in list(df.columns)[1:10]:
                wh = df[col] == 1
            else:
                wh = df[col] != '.'
            L.append([col, sum(wh)])

        wh1 = df['rna'] == 1
        wh2 = df['atac'] == 1
        wh3 = df['protein'] == 1
        wh4 = df['gwas'] == 1
        L.append(['qtl', sum(wh1 | wh2 | wh3)])
        L.append(['gwas/qtl', sum(wh1 | wh2 | wh3 | wh4)])
        L.append(['deposited to EGA', sum(wh1 | wh2 | wh3 | wh4)])

        wh1 = df['genetic_sex'] == '.'
        wh2 = df['sex'] == '.'
        df_sub = df[~(wh1 | wh2)]
        wh = df_sub['sex'].str.lower() == df_sub['genetic_sex']
        df_sub = df_sub[~wh]
        print(['sex_mismatch', df_sub.shape[0], ','.join(df_sub['donor_id'].tolist())])

        df_summary = pd.DataFrame(L, columns=['info', 'number of donors', 'comments'])
        df_summary.to_csv(out_file, index=False, sep='\t')

    def get_number_raw_peaks(
        self,
        in_dirs: list[str],
        out_file: str = 'caQTL_number_raw_peaks.txt',
    ) -> None:
        """Count raw narrowPeak calls per sample across one or more peak-call directories.

        Args:
            in_dirs: Directories to scan for ``*.narrowPeak.gz`` files; the
                peak type for each directory is derived from the text after
                the last underscore in the directory name.
            out_file: Path to write the per-sample peak-count table.
        """
        L = []
        for in_dir in in_dirs:
            for f in os.listdir(in_dir):
                if f.endswith('.narrowPeak.gz'):
                    peak_type = in_dir.split('_')[-1]
                    sample = f.split('.narrowPeak')[0]
                    n = 0
                    with gzip.open(os.path.join(in_dir, f)) as fin:
                        for line in fin:
                            n += 1
                    L.append([peak_type, sample, n])
        df = pd.DataFrame(L, columns=['peak_type', 'sample', 'number_of_peaks'])
        df.to_csv(out_file, index=False, sep='\t')

    def plot_number_raw_peaks(
        self,
        in_file: str = 'caQTL_number_raw_peaks.txt',
        cmap: str = 'Blues',
        ylabel: str = 'Number of raw peaks (million)',
        add_stripplot: bool = False,
        figsize: tuple[float, float] = (4, 4),
    ) -> None:
        """Draw a boxplot of raw peak counts per peak type.

        Args:
            in_file: Path to the table produced by ``get_number_raw_peaks``;
                the output PDF path replaces ``.txt`` with ``_boxplot.pdf``.
            cmap: Color palette for the boxes.
            ylabel: Y-axis label; if it contains ``'million'``, peak counts
                are divided by 1e6 before plotting.
            add_stripplot: Whether to overlay a strip plot of individual
                sample values.
            figsize: Figure size in inches, as ``(width, height)``.
        """
        out_file = in_file.split('.txt')[0] + '_boxplot.pdf'
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if ylabel.find('million') != -1:
            df['number_of_peaks'] = df['number_of_peaks'] / 1e6
        sns.boxplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_type', palette=cmap, legend=False)
        if add_stripplot:
            sns.stripplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, color='C0', size=4)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_number_merged_peaks(
        self,
        in_files: list[str],
        out_file: str = 'caQTL_number_merged_peaks.txt',
        params: dict[str, str] = {'consensus': 'consensus peaks', 'summitExtended': 'summit extended peaks'},
    ) -> None:
        """Count lines (peaks) in each merged-peak BED-like file.

        Args:
            in_files: Paths to merged peak files, named like
                ``<x>_<peak_type>_<peak_method>_...``; peak type and method
                are parsed from the second and third underscore-delimited
                filename tokens.
            out_file: Path to write the resulting peak-count table.
            params: Mapping used to rename raw peak-method tokens (e.g.
                ``consensus`` -> ``consensus peaks``) for display.
        """
        L = []
        for f in in_files:
            peak_type = f.split('_')[1]
            peak_method = f.split('_')[2]
            peak_method = params.get(peak_method, peak_method)
            n = 0
            with open(f) as fin:
                for line in fin:
                    n += 1
            L.append([peak_type, peak_method, n])
        df = pd.DataFrame(L, columns=['peak_type', 'peak_method', 'number_of_peaks'])
        df.to_csv(out_file, index=False, sep='\t')

    def plot_number_merged_peaks(
        self,
        in_file: str = 'caQTL_number_merged_peaks.txt',
        cmap: str = 'Blues',
        ylabel: str = 'Number of merged peaks (million)',
        figsize: tuple[float, float] = (4, 4),
    ) -> None:
        """Draw a bar chart of merged peak counts per peak type (and method, if multiple).

        Args:
            in_file: Path to the table produced by ``get_number_merged_peaks``;
                the output PDF path replaces ``.txt`` with ``_barplot.pdf``.
            cmap: Color palette for the bars.
            ylabel: Y-axis label; if it contains ``'million'``, peak counts
                are divided by 1e6 before plotting.
            figsize: Figure size in inches, as ``(width, height)``.
        """
        out_file = in_file.split('.txt')[0] + '_barplot.pdf'
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if ylabel.find('million') != -1:
            df['number_of_peaks'] = df['number_of_peaks'] / 1e6
        if df['peak_method'].nunique() > 1:
            sns.barplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_method', palette=cmap)
            ax.legend(title=None, loc='upper left', prop={'size': 7})
        else:
            sns.barplot(
                x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_type', palette=cmap, legend=False,
            )
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def plot_correlation_number_peaks_and_reads(
        self,
        in_file_reads: str = 'caQTL_number_mapped_reads.txt',
        in_file_peaks: str = 'caQTL_number_raw_peaks_qvalue.txt',
        out_file: str = 'correlation_number_peaks_and_reads.pdf',
        cmap: str = 'Blues',
        xlabel: str = 'Number of mapped reads (million)',
        ylabel: str = 'Number of raw peaks (million)',
        figsize: tuple[float, float] = (4, 4),
        base: float = 1e6,
    ) -> None:
        """Plot a regression of peak counts against mapped-read counts per sample.

        Args:
            in_file_reads: Path to a per-sample mapped-read-count table.
            in_file_peaks: Path to a per-sample peak-count table; merged
                with ``in_file_reads`` on ``sample``.
            out_file: Path to write the resulting PDF.
            cmap: Unused color-palette parameter kept for signature
                consistency with sibling plotting methods.
            xlabel: X-axis label.
            ylabel: Y-axis label.
            figsize: Figure size in inches, as ``(width, height)``.
            base: Divisor applied to both read and peak counts before
                plotting (e.g. to express counts in millions).
        """
        df_reads = pd.read_table(in_file_reads, header=0, sep='\t')
        df_peaks = pd.read_table(in_file_peaks, header=0, sep='\t')
        df = pd.merge(df_reads, df_peaks, on='sample')
        df['number_of_peaks'] = df['number_of_peaks'] / base
        df['number_of_reads'] = df['number_of_reads'] / base
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.regplot(x='number_of_reads', y='number_of_peaks', data=df, ax=ax, color='C0')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_length_distribution_merged_peaks(
        self,
        in_files: list[str] = [],
        out_file: str = 'caQTL_length_distribution_merged_peaks.txt',
        params: dict[str, str] = {'consensus': 'consensus peaks', 'summitExtended': 'summit extended peaks'},
    ) -> None:
        """Tabulate the length of every interval in a set of merged-peak BED files.

        Args:
            in_files: Paths to BED-like peak files; the peak threshold is
                parsed from the second underscore-delimited filename token
                and the peak type from the second-to-last token (with
                ``.bed`` stripped), then renamed via ``params``.
            out_file: Path to write the per-interval length table.
            params: Mapping used to rename raw peak-type tokens (e.g.
                ``consensus`` -> ``consensus peaks``) for display.
        """
        L = []
        for f in in_files:
            peak_threshold = f.split('_')[1]
            peak_type = f.split('_')[-2].split('.bed')[0]
            peak_type = params.get(peak_type, peak_type)
            with open(f) as fin:
                for line in fin:
                    items = line.strip().split('\t')
                    length = int(items[2]) - int(items[1])
                    L.append([peak_type, peak_threshold, length])
        df = pd.DataFrame(L, columns=['peak_type', 'peak_threshold', 'length'])
        df.to_csv(out_file, index=False, sep='\t')

    def plot_length_distribution_merged_peaks(
        self,
        in_file: str = 'caQTL_length_distribution_merged_peaks.txt',
        cmap: str = 'Set2',
        xlabel: str = 'Length of merged peaks (bp)',
        figsize: tuple[float, float] = (4, 4),
        peak_types: list[str] = ['consensus peaks', 'summit extended peaks'],
        show_summit: bool = True,
    ) -> None:
        """Plot a KDE of merged-peak lengths, optionally marking the summit-extension length.

        Args:
            in_file: Path to the table produced by
                ``get_length_distribution_merged_peaks``; the output PDF
                path replaces ``.txt`` with ``_hist.pdf`` (and
                ``_with_summit.pdf`` if ``show_summit`` is set).
            cmap: Seaborn color palette name.
            xlabel: X-axis label.
            figsize: Figure size in inches, as ``(width, height)``.
            peak_types: Two-element list naming the "main" peak type (index
                0, plotted as a KDE per threshold) and the "summit
                extended" peak type (index 1, used for the reference line).
            show_summit: Whether to draw a vertical dashed line at the
                fixed length of the first ``peak_types[1]`` interval.
        """
        out_file = in_file.split('.txt')[0] + '_hist.pdf'
        cmap = sns.color_palette(cmap)
        df = pd.read_table(in_file, header=0, sep='\t')
        df1 = df[df['peak_type'] == peak_types[0]]
        df2 = df[df['peak_type'] == peak_types[1]]
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.kdeplot(x='length', ax=ax, data=df1, hue='peak_threshold', palette=cmap, legend=True)
        ax.get_legend().set_title("")
        if show_summit:
            x = df2['length'].iloc[0]
            y = ax.get_ylim()[1] * 0.9
            ax.plot([x, x], [0, y], color=cmap[1], label=peak_types[1], linestyle='--', lw=2)
            out_file = out_file.split('.pdf')[0] + '_with_summit.pdf'
        ax.set_xlabel(xlabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_number_independent_signals(
        self,
        in_files: list[str] = [],
        out_file: str = 'QTL_number_of_independent_signals.txt',
    ) -> None:
        """Tabulate the distribution of independent-signal counts per QTL type.

        Args:
            in_files: Paths to per-QTL-type conditional-analysis result
                files with an ``n_independent_signals`` column; the QTL
                type is taken from the text before the first underscore in
                each filename.
            out_file: Path to write the concatenated value-count table.
        """
        L = []
        for f in in_files:
            df = pd.read_table(f, header=0, sep='\t')
            counts = df['n_independent_signals'].value_counts().to_frame().reset_index()
            counts['qtl_type'] = f.split('_')[0]
            L.append(counts)
        df = pd.concat(L, axis=0)
        df.to_csv(out_file, index=False, sep='\t')

    def plot_number_independent_signals(
        self,
        in_file: str = 'QTL_number_of_independent_signals.txt',
        show_numbers: bool = True,
        ylim: list[float] = [0, 10000],
        title: str = 'QTL conditional analysis',
        cmap: str = 'Dark2',
    ) -> None:
        """Plot a bar chart of independent-signal counts, colored and grouped by QTL type.

        Args:
            in_file: Path to the table produced by
                ``get_number_independent_signals``; the output PDF path
                replaces ``.txt`` with ``_barplot.pdf``.
            show_numbers: Whether to annotate each bar with its count.
            ylim: Y-axis limits as ``[low, high]``.
            title: Plot title.
            cmap: Seaborn color palette name used to assign one color per
                QTL type.
        """
        out_file = in_file.split('.txt')[0] + '_barplot.pdf'
        df = pd.read_table(in_file, header=0, sep='\t').reset_index()

        cmap = sns.color_palette(cmap)
        color_map = {qtl_type: cmap[i] for i, qtl_type in enumerate(df['qtl_type'].unique())}
        palette = list(df['qtl_type'].map(color_map))

        fig = plt.figure()
        ax = fig.add_subplot()
        sns.barplot(x='index', y='count', data=df, ax=ax, palette=palette, hue='index', legend=False)
        ax.set_xticks(range(len(df['n_independent_signals'])))
        ax.set_xticklabels(df['n_independent_signals'])
        ax.set_xlabel('Number of independent signals')
        ax.set_ylabel('Count')
        ax.set_ylim(ylim)
        ax.set_title(title)
        if show_numbers:
            for i, row in df.iterrows():
                ax.text(i, row['count'], row['count'], ha='center', va='bottom')
        legends = []
        for qtl_type in df['qtl_type'].unique():
            legends.append(mpatches.Patch(color=color_map[qtl_type], label=qtl_type))
        ax.legend(handles=legends)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_sig_variants(self, in_file: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz') -> None:
        """Extract the unique set of significant variants (id, chrom, position) from a result file.

        Args:
            in_file: Path to a gzip-compressed tab-delimited association
                file with ``var_id``, ``var_chr``, and ``var_from``
                columns; the output path replaces ``.txt.gz`` with
                ``_variants.txt``.
        """
        out_file = in_file.replace('.txt.gz', '_variants.txt')
        S = set()
        with gzip.open(in_file, 'rt') as f:
            head = f.readline().strip().split('\t')
            id_idx = head.index('var_id')
            chrom_idx = head.index('var_chr')
            pos_idx = head.index('var_from')
            for line in f:
                items = line.strip().split('\t')
                var = (items[id_idx], items[chrom_idx], items[pos_idx])
                S.add(var)
        with open(out_file, 'w') as f:
            for k in sorted(S):
                f.write('\t'.join(k) + '\n')

    def get_non_sig_variants(
        self,
        in_file: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo.txt.gz',
        params: dict[str, Any] = {'p_col': 'nom_pval', 'p_threshold': 0.05},
    ) -> None:
        """Extract variants whose minimum p-value across all tests exceeds a threshold.

        For each variant (identified by id, chromosome, and position), the
        minimum p-value in column ``params['p_col']`` across all its
        occurrences in ``in_file`` is compared to
        ``params['p_threshold']``; variants where the minimum exceeds the
        threshold are considered non-significant and written out.

        Args:
            in_file: Path to a gzip-compressed tab-delimited association
                file; the output path replaces ``.txt.gz`` with
                ``_non_sig_variants.txt``.
            params: Dict with keys ``p_col`` (p-value column name) and
                ``p_threshold`` (float cutoff).
        """
        out_file = in_file.replace('.txt.gz', '_non_sig_variants.txt')
        D = {}
        p_col = params.get('p_col', 'nom_pval')
        p_threshold = params.get('p_threshold', 0.05)
        with gzip.open(in_file, 'rt') as f:
            head = f.readline().strip().split('\t')
            p_idx = head.index(p_col)
            id_idx = head.index('var_id')
            chrom_idx = head.index('var_chr')
            pos_idx = head.index('var_from')
            for line in f:
                items = line.strip().split('\t')
                var = '\t'.join([items[id_idx], str(items[chrom_idx]), str(items[pos_idx])])
                try:
                    p = float(items[idx])
                except:
                    p = 1.0
                D.setdefault(var, [])
                D[var].append(p)
        with open(out_file, 'w') as f:
            for k in sorted(D):
                if min(D[k]) > p_threshold:
                    f.write(k + '\n')

    def get_variant_annotation(
        self,
        in_file: str = 'pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_variants.txt',
        vep_file: str = 'pQTL_genotyping_sampleRenamed_rsID_variantFiltered_vep.vcf.gz',
        canonical_transcript_only: bool = True,
    ) -> None:
        """Annotate variants with VEP consequence terms parsed from a VEP-annotated VCF.

        Parses the ``CSQ`` INFO field format from the VCF header, then for
        each variant collects consequence terms from the ``CSQ`` entries
        (optionally restricted to the canonical transcript), and writes
        each input variant's row with an appended sorted, comma-joined
        annotation column (defaulting to ``intergenic_variant`` when no
        annotation is found).

        Args:
            in_file: Path to a tab-delimited file of variants (first
                column is the variant id matching VCF ``ID``); the output
                path replaces ``.txt`` with ``_annotated.txt``.
            vep_file: Path to a gzip-compressed VEP-annotated VCF.
            canonical_transcript_only: If True, only consequence terms from
                transcripts flagged ``CANONICAL=YES`` are kept.
        """
        D = {}
        with gzip.open(vep_file, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    if line.find('ID=CSQ') != -1:
                        csq_format = line.split('Format: ')[-1].split('"')[0].split('|')
                        idx_canonical = csq_format.index('CANONICAL')
                        idx_consequence = csq_format.index('Consequence')
                else:
                    fields = line.split('\t')
                    var_id = fields[2]
                    D.setdefault(var_id, [])
                    fds = fields[7].split(';')
                    for fd in fds:
                        if fd.startswith('CSQ='):
                            csq = fd.split('CSQ=')[-1].split(',')
                            for item in csq:
                                tm = item.split('|')
                                if canonical_transcript_only:
                                    if tm[idx_canonical] == 'YES':
                                        D[var_id] += tm[idx_consequence].split('&')
                                else:
                                    D[var_id] += tm[idx_consequence].split('&')

        out_file = in_file.replace('.txt', '_annotated.txt')
        with open(in_file) as fin, open(out_file, 'w') as fout:
            for line in fin:
                items = line.strip().split('\t')
                var_id = items[0]
                annotation = 'intergenic_variant'
                if var_id in D and D[var_id]:
                    annotation = ','.join(sorted(set(D[var_id])))
                fout.write('\t'.join(items + [annotation]) + '\n')

    def count_variant_consequence(
        self,
        in_files: list[str] = ['pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_annotated.txt'],
        out_file: str = 'QTL_variants_consequence_count.txt',
        extra_class: list[str] = ['splice', 'UTR', 'inframe', 'frameshift', 'stop', 'start'],
    ) -> None:
        """Count variants per VEP consequence annotation, plus coarse extra-class rollups.

        For each annotated variant file, tallies how many variants carry
        each individual consequence term (from the last, comma-delimited
        annotation column written by ``get_variant_annotation``), and
        additionally accumulates counts into broader categories listed in
        ``extra_class`` whenever a term contains that category substring
        (e.g. any ``*splice*`` term also increments the ``splice`` count).

        Args:
            in_files: Paths to annotated variant files (tab-delimited,
                last column is the annotation string).
            out_file: Path to write the combined in/out count table.
            extra_class: Substrings defining coarse consequence categories
                to additionally tally.
        """
        L = []
        for f in in_files:
            print(f'processing {f}...')
            D = {}
            N = 0
            with open(f) as fin:
                for line in fin:
                    N += 1
                    line = line.strip()
                    fields = line.split('\t')
                    fds = fields[-1].split(',')
                    for item in fds:
                        D.setdefault(item, 0)
                        D[item] += 1
                        if extra_class:
                            for c in extra_class:
                                if item.find(c) != -1:
                                    D.setdefault(c, 0)
                                    D[c] += 1
            for k in sorted(D):
                L.append([f, k, D[k], N - D[k]])
        df = pd.DataFrame(L, columns=['file', 'annotation', 'in_count', 'out_count'])
        df.to_csv(out_file, index=False, sep='\t')

    def split_Ensembl_regulatory_annotation(
        self,
        in_file: str = 'Ensembl_BioMart_RegulatoryAnnotation.txt.gz',
    ) -> None:
        """Split a combined Ensembl BioMart regulatory annotation file into one BED per feature type.

        Args:
            in_file: Path to the gzip-compressed BioMart export with
                ``Feature type``, ``Chromosome/scaffold name``,
                ``Start (bp)``, and ``End (bp)`` columns; one output BED
                file per distinct feature type is written, named by
                replacing ``.txt.gz`` with ``_<feature_type>.bed``.
        """
        # download the regulatory annotation from Ensembl BioMart manually, then split the file into different feature types for downstream analysis
        df = pd.read_table(in_file, header=0, sep='\t', low_memory=False)
        for gi, g in df.groupby('Feature type'):
            out_file = in_file.replace('.txt.gz', f'_{gi}.bed')
            g_sub = g.loc[:, ['Chromosome/scaffold name', 'Start (bp)', 'End (bp)']]
            g_sub.columns = ['ch', 'start', 'end']
            g_sub['ch'] = [f'chr{x}' for x in g_sub['ch']]
            g_sub.to_csv(out_file, header=False, index=False, sep='\t')

    def count_variant_regulatory(
        self,
        in_files: list[str] = ['pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_variants.txt'],
        bed_files: list[str] = ['Ensembl_BioMart_RegulatoryAnnotation_Enhancer.bed'],
        out_file: str = 'QTL_variants_regulatory_count.txt',
    ) -> None:
        """Count variants overlapping each regulatory-annotation BED file using PyRanges.

        Args:
            in_files: Paths to headerless tab-delimited files of
                ``rsID, Chromosome, Start`` variant records (point
                positions).
            bed_files: Paths to headerless BED files of regulatory regions;
                start coordinates are shifted by 1 to convert from BED's
                0-based to 1-based coordinates before intersecting.
            out_file: Path to write the combined in/out overlap count
                table across all ``in_files`` x ``bed_files`` pairs.
        """
        L = []
        for f1 in in_files:
            for f2 in bed_files:
                print(f'processing {f1} and {f2}...')
                df1 = pd.read_table(f1, header=None, sep='\t')
                df1.columns = ['rsID', 'Chromosome', 'Start']
                df1['End'] = df1['Start']

                df2 = pd.read_table(f2, header=None, sep='\t')
                df2 = df2.iloc[:, 0:3]
                df2.columns = ['Chromosome', 'Start', 'End']
                df2['Start'] = df2['Start'] + 1

                pr1 = pyranges.PyRanges(df1)
                pr2 = pyranges.PyRanges(df2)

                pi = pr1.intersect(pr2)
                N = len(pi.rsID.unique())
                L.append([f1, f2, N, df1.shape[0] - N])
        df = pd.DataFrame(L, columns=['file', 'annotation', 'in_count', 'out_count'])
        df.to_csv(out_file, index=False, sep='\t')

    def test_enrichment_using_fisher_exact(
        self,
        in_file: str = 'QTL_variants_regulatory_count.txt',
        idx_qtl: int = 0,
    ) -> None:
        """Test enrichment of significant vs. non-significant variants in each annotation via Fisher's exact test.

        For each QTL type and annotation combination, pairs the
        significant-variant in/out counts with the corresponding
        non-significant-variant counts (identified by filename containing
        ``non_sig``, sorted to put significant rows first) and runs
        ``scipy.stats.fisher_exact`` on the resulting 2x2 contingency
        table to get an odds ratio and p-value. Groups without exactly one
        significant and one non-significant row are skipped with a
        warning.

        Args:
            in_file: Path to the table produced by
                ``count_variant_regulatory`` (or an equivalent count
                table); the output path replaces ``.txt`` with
                ``_enrichment.txt``.
            idx_qtl: Index of the underscore-delimited filename token used
                to derive the QTL type from each row's ``file`` value.
        """
        out_file = in_file.replace('.txt', '_enrichment.txt')
        L = []
        df = pd.read_table(in_file, header=0, sep='\t')
        df['qtl_type'] = [x.split('_')[idx_qtl] for x in df['file']]
        df['order'] = [1 if x.find('non_sig') != -1 else 0 for x in df['file']]
        df.sort_values('order', inplace=True)
        for gi, g in df.groupby(['qtl_type', 'annotation']):
            if g.shape[0] == 2:
                in_count_sig = g['in_count'].iloc[0]
                out_count_sig = g['out_count'].iloc[0]
                in_count_non_sig = g['in_count'].iloc[1]
                out_count_non_sig = g['out_count'].iloc[1]
                table = [[in_count_sig, out_count_sig], [in_count_non_sig, out_count_non_sig]]
                oddsratio, pvalue = scipy.stats.fisher_exact(table)
                L.append([gi[0], gi[1], oddsratio, pvalue])
            else:
                print(f'Annotation {gi} needs to be checked')
        df_out = pd.DataFrame(L, columns=['qtl', 'annotation', 'odds_ratio', 'p_value'])
        df_out.to_csv(out_file, index=False, sep='\t')

    def bar_plot_enrichment(
        self,
        in_file: str = 'QTL_variants_regulatory_count_enrichment.txt',
        subset_renaming_file: str = 'subset_renaming.txt',
        xlim: list[float] = [0, 10],
        cmap: str = 'Dark2',
        title: str = 'Enrichment of QTL significant variants',
    ) -> None:
        """Plot odds ratios of variant-annotation enrichment as a horizontal bar chart.

        Args:
            in_file: Path to the table produced by
                ``test_enrichment_using_fisher_exact``; the output PDF
                path replaces ``.txt`` with ``_barplot.pdf``.
            subset_renaming_file: Optional headerless two-column
                tab-delimited file mapping raw annotation names to display
                names; when present, also filters and reorders annotations
                to those listed.
            xlim: X-axis (odds ratio) limits as ``[low, high]``; skipped if
                falsy.
            cmap: Seaborn color palette name.
            title: Plot title.
        """
        D = {}
        if subset_renaming_file and os.path.exists(subset_renaming_file):
            df_subset = pd.read_table(subset_renaming_file, header=None, sep='\t')
            D = dict(zip(df_subset.iloc[:, 0], df_subset.iloc[:, 1]))

        out_file = in_file.replace('.txt', '_barplot.pdf')
        df = pd.read_table(in_file, header=0, sep='\t')
        if D:
            wh = df['annotation'].isin(D)
            df = df[wh].copy()
            df.sort_values('annotation', key=lambda x: [list(D.keys()).index(i) for i in x], inplace=True)
            df['annotation'] = df['annotation'].map(D)

        fig = plt.figure()
        ax = fig.add_subplot()
        sns.barplot(y='annotation', x='odds_ratio', hue='qtl', data=df, ax=ax, palette=cmap)
        ylim = ax.get_ylim()
        ax.plot([1, 1], ylim, '--', color='orange', lw=2)
        if xlim:
            ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(title)
        ax.set_ylabel('')
        ax.legend(title=None, loc='lower right')
        plt.tight_layout()
        plt.savefig(out_file)

    def get_nominal_sig_associations(
        self,
        in_files: list[str] = [
            'caQTL_nominal-1.0_w1k_qvalue_extraInfo_sig.txt.gz',
            'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz',
            'pQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz',
        ],
        out_file: str = 'QTL_nomnial_sig_associations.txt',
        qtl_types: dict[str, str] = {},
    ) -> None:
        """Pool significant nominal associations from multiple QTL types into one long table.

        For ``caQTL`` files, the phenotype id's trailing comma-separated
        gene list is exploded into one row per gene; for other QTL types
        the phenotype id's trailing token is used directly as the gene.

        Args:
            in_files: Paths to gzip-compressed tab-delimited significant
                nominal-association files, one per QTL type.
            out_file: Path to write the pooled association table.
            qtl_types: Optional mapping from file path to an explicit QTL
                type label; falls back to the first underscore-delimited
                filename token when a file is not present in the mapping.
        """
        L = []
        for f in in_files:
            qtl_type = qtl_types.get(f, f.split('_')[0])
            with gzip.open(f, 'rt') as fin:
                head = fin.readline().strip().split('\t')
                phe_idx = head.index('phe_id')
                var_idx = head.index('var_id')
                p_idx = head.index('nom_pval')
                beta_idx = head.index('slope')
                for line in fin:
                    items = line.strip().split('\t')
                    if qtl_type.find('caQTL') != -1:
                        genes = items[phe_idx].split('_')[-1].split(',')
                    else:
                        genes = [items[phe_idx].split('_')[-1]]
                    for gene in genes:
                        L.append([items[phe_idx], items[var_idx], gene, items[beta_idx], items[p_idx], qtl_type])
        df = pd.DataFrame(L, columns=['phe_id', 'var_id', 'gene', 'beta', 'pval', 'qtl'])
        df.to_csv(out_file, index=False, sep='\t')

    def get_recurrent_associatoins(
        self,
        in_file: str = 'QTL_nominal_sig_associations.txt',
        N: int = 3,
    ) -> None:
        """Find gene-variant associations recurring across at least N QTL types.

        Groups pooled associations (from ``get_nominal_sig_associations``)
        by gene and then by ``(gene, var_id)`` pair, keeping pairs observed
        in at least ``N`` distinct QTL types (deduplicated to the
        lowest-p-value row per QTL type). For each gene, the single
        best-supported pair (by minimum p-value) is retained in the
        output.

        Args:
            in_file: Path to the pooled nominal significant-association
                table; the output path replaces ``.txt`` with
                ``_recurrent{N}.txt``.
            N: Minimum number of distinct QTL types a gene-variant pair
                must appear in to be considered recurrent.
        """
        out_file = in_file.replace('.txt', f'_recurrent{N}.txt')
        df = pd.read_table(in_file, header=0, sep='\t')
        D = {}
        for gi, g in df.groupby('gene'):
            for gi2, g2 in g.groupby(['gene', 'var_id']):
                if g2['qtl'].nunique() >= N:
                    D.setdefault(gi, [])
                    g3 = g2.sort_values('pval')
                    g3.drop_duplicates(subset='qtl', keep='first', inplace=True)
                    D[gi].append(g3)

        for k in D:
            D[k] = sorted(D[k], key=lambda x: x['pval'].min())

        L = []
        for k in sorted(D):
            L.append(D[k][0])
        if L:
            df_out = pd.concat(L, axis=0)
            df_out.columns = df.columns
            df_out.to_csv(out_file, index=False, sep='\t')

    def plot_heatmap_of_recurrent_associatoins(
        self,
        in_file: str = 'QTL_nominal_sig_associations_recurrent3.txt',
        cmap: str = 'coolwarm',
        figsize: tuple[float, float] = (4, 8),
        fontsize: float = 8,
        customize_cbar: bool = True,
        genes_highlight: list[str] = ['PTGFRN', 'STARD10', 'PEPD'],
    ) -> None:
        """Draw a clustered heatmap of effect sizes for recurrent multi-QTL-type associations.

        Pivots the recurrent-association table to a gene/variant x QTL-type
        matrix of beta values (missing combinations filled with 0), draws
        a row-clustered heatmap without dendrograms, bolds highlighted gene
        labels, annotates each row with its variant id, and optionally
        repositions the colorbar above the heatmap.

        Args:
            in_file: Path to the table produced by
                ``get_recurrent_associatoins``; the output PDF path
                replaces ``.txt`` with ``_heatmap.pdf``.
            cmap: Colormap for the heatmap.
            figsize: Figure size in inches, as ``(width, height)``.
            fontsize: Font size for row (gene/variant) labels.
            customize_cbar: Whether to reposition the colorbar to a small
                horizontal bar above the heatmap.
            genes_highlight: Gene names whose y-axis tick labels should be
                bolded.
        """
        out_file = in_file.replace('.txt', '_heatmap.pdf')
        df = pd.read_table(in_file, header=0, sep='\t')
        df_pivot = df.pivot(index=['gene', 'var_id'], columns='qtl', values='beta')
        df_pivot.fillna(0, inplace=True)

        df = pd.DataFrame(df_pivot.values)
        df.columns = df_pivot.columns
        df.index = df_pivot.index.get_level_values('gene')
        df.columns.name = None
        df.index.name = None
        df['variant'] = df_pivot.index.get_level_values('var_id')
        g = sns.clustermap(df.iloc[:, 0:-1], cmap=cmap, yticklabels=True, figsize=figsize)
        g.ax_row_dendrogram.set_visible(False)
        g.ax_col_dendrogram.set_visible(False)
        ax = g.ax_heatmap
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize)
        for tick in ax.get_yticklabels():
            if tick.get_text() in genes_highlight:
                tick.set_fontweight('bold')

        row_indices = g.dendrogram_row.reordered_ind
        for i, idx in enumerate(row_indices):
            label = df['variant'].iloc[idx]
            ax.text(0 - 0.02, i + 0.5, label, ha='right', va='center', fontsize=fontsize)

        if customize_cbar:
            hm_pos = g.ax_heatmap.get_position()
            cax_left = (hm_pos.x0 + hm_pos.x1) / 2 - 0.1
            cax_bottom = hm_pos.y1 + 0.01
            cax_width = 0.2
            cax_height = 0.02

            cax = g.cax
            cax.set_position([cax_left, cax_bottom, cax_width, cax_height])
            cax.clear()
            plt.colorbar(g.ax_heatmap.get_children()[0], cax=cax, orientation='horizontal')
            cax.xaxis.tick_top()
            cax.xaxis.set_label_position('top')
            cax.set_xlabel('beta')

        plt.savefig(out_file)

    def bar_plot_sig_count_by_params(
        self,
        in_file: str,
        contrast: list[str] = ['consensus_peaks', 'summit_extended_peaks'],
        out_suffix: str = 'summit_vs_consensus',
        figsize: tuple[float, float] = (4, 4),
        cmap: str = 'Set2',
        ylim: list[float] = [0, 8000],
        xticklabels: list[str] = [],
        fontsize: float = 8,
    ) -> None:
        """Plot a grouped bar chart of significant-peak counts across two parameter settings.

        Args:
            in_file: Path to a headerless tab-delimited file whose columns
                are ``[count, file, param1, param2]``; the output PDF path
                is derived by appending ``_{out_suffix}.pdf`` after
                stripping ``.txt``.
            contrast: The two ``param1`` values to compare; rows with other
                values are dropped, and rows are ordered to match this
                list.
            out_suffix: Suffix appended to the output PDF filename.
            figsize: Figure size in inches, as ``(width, height)``.
            cmap: Seaborn color palette name.
            ylim: Y-axis limits as ``[low, high]``.
            xticklabels: Optional custom x-axis tick labels.
            fontsize: Font size for the bar-value labels.
        """
        df = pd.read_table(in_file, header=None, sep='\t')
        df.columns = ['Number of significant peaks', 'file', 'param1', 'param2']

        df = df[df['param1'].isin(contrast)]
        df.sort_values(by='param1', inplace=True, key=lambda x: [contrast.index(i) for i in x])
        print(df)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.barplot(x='param1', y='Number of significant peaks', hue='param2', data=df, ax=ax, palette=cmap)
        ax.set_xlabel('')
        ax.legend(title=None)
        for i in ax.containers:
            ax.bar_label(i, fmt='%.0f', padding=3, fontsize=fontsize)
        ax.set_ylim(ylim)
        if xticklabels:
            ax.set_xticklabels(xticklabels)

        plt.tight_layout()
        plt.savefig(f'{in_file.split(".txt")[0]}_{out_suffix}.pdf')

    def prs_overlap_with_qtl(self, in_file: str, in_file2: str, flank: float = 1e6) -> None:
        """Find QTL associations overlapping PRS variants within a flanking window, per PRS variant.

        For each PRS variant in ``in_file``, queries the tabix-indexed
        ``in_file2`` for records within ``flank`` base pairs, then inner
        joins the PRS row with matching QTL rows on rsID/``var_id``.
        Rows that raise an exception (e.g. missing tabix contig) are
        silently skipped.

        Args:
            in_file: Path to a CSV of PRS variants with ``contig_id``,
                ``position_hg38``, and ``rsid`` columns.
            in_file2: Path to a tabix-indexed, gzip-compressed tab-delimited
                QTL association file with a header row and a ``var_id``
                column.
            flank: Half-width, in base pairs, of the genomic window queried
                around each PRS variant's position.
        """
        df = pd.read_table(in_file, header=0, sep=',')
        tb = tabix.open(in_file2)
        df2_cols = pd.read_table(in_file2, header=0, sep='\t', nrows=1).columns
        L = []
        for n in range(df.shape[0]):
            try:
                chrom = df['contig_id'].iloc[n]
                pos = df['position_hg38'].iloc[n]
                rs = df['rsid'].iloc[n]
                start = max(int(pos - flank), 0)
                end = int(pos + flank)
                res = tb.query(chrom, start, end)
                df2 = pd.DataFrame(res)
                if df2.shape[0] > 0:
                    df2.columns = df2_cols
                    df_merged = pd.merge(df.iloc[[n]], df2, left_on='rsid', right_on='var_id', how='inner')
                    L.append(df_merged)
            except Exception as e:
                pass
        if len(L) > 0:
            df = pd.concat(L, axis=0)
            out_file = in_file.split('.csv')[0] + '_' + in_file2.split('.txt')[0] + '_overlap.txt'
            df.to_csv(out_file, index=False, sep='\t')

    def count_overlap_with_prs(
        self,
        in_files: list[str] = ['t2dp_suzuki24_ma_eQTL_nominal-1.0_w1M_PC25_extraInfo_sig_overlap.txt'],
        out_file: str = 'QTL_variants_prs_count_t2dp_suzuki24.txt',
        fsep: str = '_ma_',
    ) -> None:
        """Count, per PRS group, how many significant/non-significant QTL variants overlap PRS hits.

        For each overlap file (produced by ``prs_overlap_with_qtl`` on the
        significant-variant subset), locates the matching significant and
        non-significant variant-id lists and the corresponding full/",
        non-significant overlap files, then for every PRS ``group`` counts
        how many of the significant and non-significant variants are
        represented in the overlap.

        Args:
            in_files: Paths to significant-variant PRS overlap files (named
                ``..._sig_overlap.txt``).
            out_file: Path to write the combined in/out count table.
            fsep: Separator substring used to strip the PRS-specific prefix
                off the filename to recover the base QTL variant-list
                filenames.
        """
        L = []
        for f in in_files:
            sig_variants_file = f.split(fsep)[-1].replace('_sig_overlap.txt', '_sig_variants.txt')
            non_sig_variants_file = f.split(fsep)[-1].replace('_sig_overlap.txt', '_non_sig_variants.txt')
            sig_variants = set(pd.read_table(sig_variants_file, header=None, sep='\t').iloc[:, 0].values)
            non_sig_variants = set(pd.read_table(non_sig_variants_file, header=None, sep='\t').iloc[:, 0].values)

            f_all = f.replace('_sig_overlap.txt', '_overlap.txt')
            f_non_sig = f.replace('_sig_overlap.txt', '_non_sig_overlap.txt')
            df_sig = pd.read_table(f, header=0, sep='\t')
            df_all = pd.read_table(f_all, header=0, sep='\t')

            for g in df_all['group'].unique():
                df_sig_sub = df_sig[(df_sig['group'] == g) & df_sig['var_id'].isin(sig_variants)]
                df_non_sig_sub = df_all[(df_all['group'] == g) & df_all['var_id'].isin(non_sig_variants)]
                n_sig = df_sig_sub['var_id'].nunique()
                n_non_sig = df_non_sig_sub['var_id'].nunique()
                L.append([f, g, n_sig, len(sig_variants) - n_sig])
                L.append([f_non_sig, g, n_non_sig, len(non_sig_variants) - n_non_sig])
        df = pd.DataFrame(L)
        df.columns = ['file', 'annotation', 'in_count', 'out_count']
        df.to_csv(out_file, index=False, sep='\t')

    def count_overlap_sig_pair_eQTL_GTEx(
        self,
        in_file: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz',
        in_files: list[str] = ['GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eQTLs.signif_pairs.parquet'],
        out_file: str = 'QTL_sig_pair_overlap_gtex_count.txt',
        cols: list[str] = ['file', 'eQTL_sig_pair_count', 'GTEx_sig_pair_count', 'overlap_count'],
    ) -> None:
        """Count gene-variant pair overlap between significant islet eQTL and GTEx tissue eQTL.

        Builds a ``gene_chrom_pos_ref_alt`` pair key for the islet eQTL
        table and a ``gene_variant`` pair key (stripping the GTEx build
        suffix) for each GTEx parquet file, then reports the pair counts
        in each set and their overlap.

        Args:
            in_file: Path to the gzip-compressed significant islet eQTL
                table.
            in_files: Paths to GTEx ``signif_pairs.parquet`` files, one per
                tissue.
            out_file: Path to write the overlap-count table.
            cols: Output column names.
        """
        df = pd.read_table(in_file, header=0, sep='\t')
        pairs = []
        for n in range(df.shape[0]):
            gene_id = df['phe_id'].iloc[n].split('_')[0]
            chrom = df['var_chr'].iloc[n]
            pos = df['var_from'].iloc[n]
            ref = df['non_effective_allele'].iloc[n]
            alt = df['effective_allele'].iloc[n]
            pair = f'{gene_id}_{chrom}_{pos}_{ref}_{alt}'
            pairs.append(pair)
        df['pair'] = pairs

        L = []
        for f in in_files:
            print('processing ' + f + '...', flush=True)
            df2 = pd.read_parquet(f)
            pairs = []
            for n in range(df2.shape[0]):
                gene_id = df2['phenotype_id'].iloc[n].split('.')[0]
                var_id = '_'.join(df2['variant_id'].iloc[n].split('_')[0:-1])
                pair = f'{gene_id}_{var_id}'
                pairs.append(pair)
            df2['pair'] = pairs
            L.append([os.path.basename(f), df['pair'].nunique(), df2['pair'].nunique(), df2['pair'].isin(df['pair']).sum()])
        df_out = pd.DataFrame(L, columns=cols)
        df_out.to_csv(out_file, index=False, sep='\t')

    def bar_plot_overlap_eQTL_GTEx(
        self,
        in_file: str = 'QTL_sig_pair_overlap_gtex_count.txt',
        cmap: str = 'colorblind',
        title: str = 'Overlap between eQTL and GTEx',
        subset_renaming_file: str = 'subset_renaming.txt',
        ratio_base: str = 'eQTL_sig_pair_count',
        figsize: tuple[float, float] = (4, 4),
    ) -> None:
        """Plot the percent overlap between islet eQTL and each GTEx tissue as a bar chart.

        Tissues are optionally renamed/filtered via ``subset_renaming_file``,
        otherwise derived from the first underscore-delimited token of each
        filename; bars are sorted by descending mean overlap ratio per
        tissue.

        Args:
            in_file: Path to the table produced by
                ``count_overlap_sig_pair_eQTL_GTEx``.
            cmap: Seaborn color palette name.
            title: Plot title.
            subset_renaming_file: Optional headerless two-column
                tab-delimited file mapping raw filenames to display tissue
                names.
            ratio_base: Column name used as the denominator when computing
                the percent-overlap ratio; also used to name the output
                file.
            figsize: Figure size in inches, as ``(width, height)``.
        """
        out_file = in_file.replace('.txt', f'_on_{ratio_base.split("_")[0]}_barplot.pdf')
        D = {}
        if os.path.exists(subset_renaming_file):
            df_subset = pd.read_table(subset_renaming_file, header=None, sep='\t')
            D = dict(zip(df_subset.iloc[:, 0], df_subset.iloc[:, 1]))
        df = pd.read_table(in_file, header=0, sep='\t')
        if D:
            df['Tissue'] = df['file'].map(D)
            df.dropna(subset=['Tissue'], inplace=True)
        else:
            df['Tissue'] = [x.split('_')[0] for x in df['file']]
        df['ratio'] = df['overlap_count'] / df[ratio_base] * 100

        T = {}
        for gi, g in df.groupby('Tissue'):
            T[gi] = g['ratio'].mean()
        df['ratio_mean'] = df['Tissue'].map(T)
        df.sort_values('ratio_mean', inplace=True, ascending=False)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.barplot(x='Tissue', y='ratio', data=df, ax=ax, palette=cmap, hue='Tissue', legend=False, capsize=0.1)
        ax.set_ylabel('Percent of significant\ngene-variant pairs')
        ax.set_xlabel('')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title(title, fontsize=14)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_overlap_sig_pair_pQTL_UKBBplasma(
        self,
        in_file: str = 'pQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz',
        in_file2: str = 'UKB-PPP_pQTL_Euro_sig_rsID.txt.gz',
    ) -> None:
        """Merge significant islet pQTL with significant UK Biobank plasma pQTL on gene-rsID pairs.

        Args:
            in_file: Path to the gzip-compressed significant islet pQTL
                table.
            in_file2: Path to the gzip-compressed significant UKBB plasma
                pQTL table; the output path concatenates the two input
                basenames (stripping ``.txt.gz`` from ``in_file``).
        """
        df = pd.read_table(in_file, header=0, sep='\t')
        df2 = pd.read_table(in_file2, header=0, sep='\t')

        pairs = []
        for n in range(df.shape[0]):
            gene = df['phe_id'].iloc[n].split('_')[-1]
            rs = df['var_id'].iloc[n]
            pair = f'{gene}_{rs}'
            pairs.append(pair)
        df['pair'] = pairs

        pairs2 = []
        for n in range(df2.shape[0]):
            gene = df2['Gene'].iloc[n]
            rs = df2['rsID'].iloc[n]
            pair = f'{gene}_{rs}'
            pairs2.append(pair)
        df2['pair'] = pairs2

        df = pd.merge(df, df2, on='pair', how='inner')
        out_file = in_file.split('.txt.gz')[0] + '_' + in_file2
        df.to_csv(out_file, index=False, sep='\t')

    def correlation_plot_pQTL_UKBBplasma(
        self,
        in_file: str = 'pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_UKB-PPP_pQTL_Euro_sig_rsID.txt.gz',
        beta_x: str = 'slope',
        beta_y: str = 'BETA',
        cmap: str = 'Blues',
        xlabel: str = 'beta, significant pQTL in islets',
        ylabel: str = 'beta, significant pQTL in plasma\n(pvalue < 5e-6, UKBB)',
        figsize: tuple[float, float] = (4, 4),
        xlim: list[float] = [-2, 2],
        ylim: list[float] = [-2, 2],
        line_params: dict[str, Any] = {
            'hline': [[-1.5, 1.5], [0, 0]],
            'vline': [[0, 0], [-1.5, 1.5]],
            'color': 'orange',
            'ls': '--',
            'lw': 1,
        },
        title: str | None = None,
        color: str = 'C0',
        scatter_size: float = 6,
    ) -> None:
        """Plot a regression of islet pQTL effect sizes against UK Biobank plasma pQTL effect sizes.

        Draws reference horizontal/vertical dashed lines at the origin,
        and if ``title`` is not given, computes and displays the percent
        of shared pQTL that are concordant in direction (both effects
        positive or both negative).

        Args:
            in_file: Path to the merged table produced by
                ``get_overlap_sig_pair_pQTL_UKBBplasma``; the output PDF
                path replaces ``.txt`` with ``_correlation.pdf``.
            beta_x: Column name for the islet pQTL effect size.
            beta_y: Column name for the UKBB plasma pQTL effect size.
            cmap: Unused color-palette parameter kept for signature
                consistency with sibling plotting methods.
            xlabel: X-axis label.
            ylabel: Y-axis label.
            figsize: Figure size in inches, as ``(width, height)``.
            xlim: X-axis limits as ``[low, high]``.
            ylim: Y-axis limits as ``[low, high]``.
            line_params: Dict describing the reference lines, with keys
                ``hline``/``vline`` (each a pair of x/y coordinate lists),
                ``color``, ``ls``, and ``lw``.
            title: Plot title; computed automatically (concordance
                percentage) when None.
            color: Color for the scatter/regression points.
            scatter_size: Marker size for the scatter points.
        """
        out_file = in_file.split('.txt')[0] + '_correlation.pdf'
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.regplot(x=beta_x, y=beta_y, data=df, ax=ax, color=color, scatter_kws={'s': scatter_size})

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks(range(xlim[0], xlim[1] + 1))
        ax.set_yticks(range(ylim[0], ylim[1] + 1))
        ax.plot(
            line_params['hline'][0], line_params['hline'][1],
            color=line_params['color'], ls=line_params['ls'], lw=line_params['lw'],
        )
        ax.plot(
            line_params['vline'][0], line_params['vline'][1],
            color=line_params['color'], ls=line_params['ls'], lw=line_params['lw'],
        )

        if title is None:
            wh1 = (df[beta_x] < 0) & (df[beta_y] < 0)
            wh2 = (df[beta_x] > 0) & (df[beta_y] > 0)
            df_con = df.loc[wh1 | wh2, ]
            title = f'{df_con.shape[0] / df.shape[0] * 100:.1f}% shared pQTL\nare concordant in direction'
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_overlap_sig_pair_eQTL_InsPIRE(
        self,
        in_file: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz',
        in_file2: str = 'InsPIRE_Gene_eQTL.txt',
    ) -> None:
        """Merge significant islet eQTL with the InsPIRE islet eQTL dataset on gene-rsID pairs.

        Args:
            in_file: Path to the gzip-compressed significant islet eQTL
                table.
            in_file2: Path to the InsPIRE gene-level eQTL table; the output
                path concatenates the two input basenames (stripping
                ``.txt.gz`` from ``in_file``).
        """
        df = pd.read_table(in_file, header=0, sep='\t')
        df2 = pd.read_table(in_file2, header=0, sep='\t')

        pairs = []
        for n in range(df.shape[0]):
            gene = df['phe_id'].iloc[n].split('_')[-1]
            rs = df['var_id'].iloc[n]
            pair = f'{gene}_{rs}'
            pairs.append(pair)
        df['pair'] = pairs

        pairs2 = []
        for n in range(df2.shape[0]):
            gene = df2['GeneName'].iloc[n]
            rs = df2['SNPid'].iloc[n]
            pair = f'{gene}_{rs}'
            pairs2.append(pair)
        df2['pair'] = pairs2

        df = pd.merge(df, df2, on='pair', how='inner')
        out_file = in_file.split('.txt.gz')[0] + '_' + in_file2
        df.to_csv(out_file, index=False, sep='\t')

    def correlation_plot_eQTL_InsPIRE(
        self,
        in_file: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig_InsPIRE_Gene_eQTL.txt',
        beta_x: str = 'slope',
        beta_y: str = 'Slope',
        cmap: str = 'Blues',
        xlabel: str = 'beta, significant eQTL in islets',
        ylabel: str = 'beta, significant eQTL in islets\n(InsPIRE)',
        figsize: tuple[float, float] = (4, 4),
        xlim: list[float] = [-2, 2],
        ylim: list[float] = [-2, 2],
        line_params: dict[str, Any] = {
            'hline': [[-1.5, 1.5], [0, 0]],
            'vline': [[0, 0], [-1.5, 1.5]],
            'color': 'orange',
            'ls': '--',
            'lw': 1,
        },
        title: str | None = None,
        color: str = 'C0',
        scatter_size: float = 6,
    ) -> None:
        """Plot a regression of islet eQTL effect sizes against InsPIRE islet eQTL effect sizes.

        Draws reference horizontal/vertical dashed lines at the origin,
        and if ``title`` is not given, computes and displays the percent
        of shared eQTL that are concordant in direction.

        Args:
            in_file: Path to the merged table produced by
                ``get_overlap_sig_pair_eQTL_InsPIRE``; the output PDF path
                replaces ``.txt`` with ``_correlation.pdf``.
            beta_x: Column name for the islet eQTL effect size.
            beta_y: Column name for the InsPIRE eQTL effect size.
            cmap: Unused color-palette parameter kept for signature
                consistency with sibling plotting methods.
            xlabel: X-axis label.
            ylabel: Y-axis label.
            figsize: Figure size in inches, as ``(width, height)``.
            xlim: X-axis limits as ``[low, high]``.
            ylim: Y-axis limits as ``[low, high]``.
            line_params: Dict describing the reference lines, with keys
                ``hline``/``vline`` (each a pair of x/y coordinate lists),
                ``color``, ``ls``, and ``lw``.
            title: Plot title; computed automatically (concordance
                percentage) when None.
            color: Color for the scatter/regression points.
            scatter_size: Marker size for the scatter points.
        """
        out_file = in_file.split('.txt')[0] + '_correlation.pdf'
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.regplot(x=beta_x, y=beta_y, data=df, ax=ax, color=color, scatter_kws={'s': scatter_size})

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks(range(xlim[0], xlim[1] + 1))
        ax.set_yticks(range(ylim[0], ylim[1] + 1))
        ax.plot(
            line_params['hline'][0], line_params['hline'][1],
            color=line_params['color'], ls=line_params['ls'], lw=line_params['lw'],
        )
        ax.plot(
            line_params['vline'][0], line_params['vline'][1],
            color=line_params['color'], ls=line_params['ls'], lw=line_params['lw'],
        )

        if title is None:
            wh1 = (df[beta_x] < 0) & (df[beta_y] < 0)
            wh2 = (df[beta_x] > 0) & (df[beta_y] > 0)
            df_con = df.loc[wh1 | wh2, ]
            title = f'{df_con.shape[0] / df.shape[0] * 100:.1f}% shared eQTL\nare concordant in direction'
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_overlap_sig_pair_eQTLexon_sQTL(
        self,
        in_file: str = 'eQTLexon_nominal-1.0_w100k_PC25_extraInfo_sig.txt.gz',
        in_file2: str = 'sQTL_nominal-1.0_w100k_PC25_extraInfo_sig.txt.gz',
    ) -> None:
        """Merge significant exon-level eQTL with significant sQTL on shared variant and exon id.

        Args:
            in_file: Path to the gzip-compressed significant exon-level
                eQTL table; the exon id is parsed as the fourth
                underscore-delimited token of ``phe_id``.
            in_file2: Path to the gzip-compressed significant sQTL table,
                parsed the same way; the output path concatenates the two
                input basenames (stripping ``.txt.gz`` from ``in_file``).
        """
        df1 = pd.read_table(in_file, header=0, sep='\t', low_memory=False)
        df2 = pd.read_table(in_file2, header=0, sep='\t', low_memory=False)
        df1['ExonID'] = df1['phe_id'].apply(lambda x: x.split('_')[3])
        df2['ExonID'] = df2['phe_id'].apply(lambda x: x.split('_')[3])

        df = pd.merge(df1, df2, on=['var_id', 'ExonID'], suffixes=('_eQTLexon', '_sQTL'))
        out_file = in_file.split('.txt.gz')[0] + '_' + in_file2
        df.to_csv(out_file, index=False, sep='\t')

    def correlation_plot_eQTLexon_sQTL(
        self,
        in_file: str = 'eQTLexon_nominal-1.0_w100k_PC25_extraInfo_sig_sQTL_nominal-1.0_w100k_PC25_extraInfo_sig.txt.gz',
        beta_x: str = 'slope_eQTLexon',
        beta_y: str = 'slope_sQTL',
        cmap: str = 'Blues',
        xlabel: str = 'beta, significant eQTL on exon level',
        ylabel: str = 'beta, significant sQTL',
        figsize: tuple[float, float] = (4, 4),
        xlim: list[float] = [-2, 2],
        ylim: list[float] = [-2, 2],
        line_params: dict[str, Any] = {
            'hline': [[-1.5, 1.5], [0, 0]],
            'vline': [[0, 0], [-1.5, 1.5]],
            'color': 'orange',
            'ls': '--',
            'lw': 1,
        },
        title: str | None = None,
        color: str = 'C1',
        scatter_size: float = 1,
    ) -> None:
        """Plot a regression of exon-level eQTL effect sizes against sQTL effect sizes.

        Draws reference horizontal/vertical dashed lines at the origin,
        and if ``title`` is not given, computes and displays the percent
        of shared exon associations that are concordant in direction.

        Args:
            in_file: Path to the merged table produced by
                ``get_overlap_sig_pair_eQTLexon_sQTL``; the output PDF path
                replaces ``.txt`` with ``_correlation.pdf``.
            beta_x: Column name for the exon-level eQTL effect size.
            beta_y: Column name for the sQTL effect size.
            cmap: Unused color-palette parameter kept for signature
                consistency with sibling plotting methods.
            xlabel: X-axis label.
            ylabel: Y-axis label.
            figsize: Figure size in inches, as ``(width, height)``.
            xlim: X-axis limits as ``[low, high]``.
            ylim: Y-axis limits as ``[low, high]``.
            line_params: Dict describing the reference lines, with keys
                ``hline``/``vline`` (each a pair of x/y coordinate lists),
                ``color``, ``ls``, and ``lw``.
            title: Plot title; computed automatically (concordance
                percentage) when None.
            color: Color for the scatter/regression points.
            scatter_size: Marker size for the scatter points.
        """
        out_file = in_file.split('.txt')[0] + '_correlation.pdf'
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.regplot(x=beta_x, y=beta_y, data=df, ax=ax, color=color, scatter_kws={'s': scatter_size})

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks(range(xlim[0], xlim[1] + 1))
        ax.set_yticks(range(ylim[0], ylim[1] + 1))
        ax.plot(
            line_params['hline'][0], line_params['hline'][1],
            color=line_params['color'], ls=line_params['ls'], lw=line_params['lw'],
        )
        ax.plot(
            line_params['vline'][0], line_params['vline'][1],
            color=line_params['color'], ls=line_params['ls'], lw=line_params['lw'],
        )

        if title is None:
            wh1 = (df[beta_x] < 0) & (df[beta_y] < 0)
            wh2 = (df[beta_x] > 0) & (df[beta_y] > 0)
            df_con = df.loc[wh1 | wh2, ]
            title = f'{df_con.shape[0] / df.shape[0] * 100:.1f}% shared exon\nare concordant in direction'
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(out_file)
