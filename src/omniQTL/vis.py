"""Genome-browser-style visualization utilities for omniQTL.

This module provides gene/transcript track plots rendered from BED12 files,
LocusZoom-style association scatter plots with LD-based coloring, genome-wide
Manhattan plots, genotype-vs-phenotype bar plots, pyGenomeTracks configuration
generation, and UCSC-track-hub / bed-track export helpers for QTL mapping
results and GWAS summary statistics.
"""
from typing import Any

from .utils import *
from matplotlib.patches import Rectangle, Patch


class GeneTxPlot():
    """Gene/transcript track plotting from BED12 annotation files.

    Provides utilities to read a BED12 gene-model file, subset it to a single
    gene plus its flanking region, and render transcript/exon tracks onto a
    matplotlib axis. Also implements the interval-packing helper used to lay
    out overlapping transcripts on separate display rows.
    """

    def __init__(self) -> None:
        """Initialize the plotter (no state is stored)."""
        pass

    def read_bed12(self, bed12: str = 'Homo_sapiens.GRCh38.115.bed12') -> pd.DataFrame:
        """Read a BED12 gene-model annotation file into a DataFrame.

        Args:
            bed12: Path to the BED12 file with gene/transcript annotations.

        Returns:
            DataFrame with named BED12 columns (chrom, chromStart, chromEnd,
            name, score, strand, thickStart, thickEnd, itemRgb, blockCount,
            blockSizes, blockStarts, geneID, geneName, geneBiotype,
            transcriptID, transcriptName, transcriptBiotype, tag).
        """
        df = pd.read_table(bed12, header=None, sep='\t', low_memory=False)
        df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                        'geneID', 'geneName', 'geneBiotype', 'transcriptID', 'transcriptName', 'transcriptBiotype', 'tag']
        return df

    def subset_gene(
        self,
        df: pd.DataFrame,
        gene: str = 'GINS4',
        canonical_only: bool = True,
        flank: float = 1e6,
        geneBiotype_include: list[str] = ['protein_coding'],
        geneName_exclude: list[str] = ['.'],
    ) -> tuple[pd.DataFrame, list]:
        """Subset a BED12 annotation DataFrame to one gene and its flanking window.

        Args:
            df: BED12 annotation DataFrame, as returned by `read_bed12`.
            gene: Gene symbol (matched against `geneName`) to center the window on.
            canonical_only: If True, keep only rows whose `tag` field contains
                'canonical'.
            flank: Number of bases to extend on each side of the gene's
                min/max coordinates to build the display window.
            geneBiotype_include: Gene biotypes to keep in the windowed subset.
            geneName_exclude: Gene names to drop from the windowed subset.

        Returns:
            Tuple of `(df_sub, window)` where `df_sub` is the annotation
            subset within the window and `window` is a 3-element list
            `[chrom, start, end]` (chrom as a `'chr1'`-style string, start/end
            as ints).
        """
        wh = df['geneName'].isin([gene])
        df_sub = df[wh]
        if canonical_only:
            wh = np.array([True if x.find('canonical') != -1 else False for x in df_sub['tag']])
            df_sub = df_sub[wh]

        chrom = str(df_sub['chrom'].iloc[0])
        start = df_sub['chromStart'].min()
        end = df_sub['chromEnd'].max()
        window = ['chr' + chrom, int(start - flank), int(end + flank)]

        wh1 = (df['chrom'] == chrom) & (df['chromStart'] > window[1]) & (df['chromEnd'] < window[2])
        wh2 = df['geneBiotype'].isin(geneBiotype_include)
        wh3 = ~df['geneName'].isin(geneName_exclude)
        df_sub = df[wh1 & wh2 & wh3]
        if canonical_only:
            wh = np.array([True if x.find('canonical') != -1 else False for x in df_sub['tag']])
            df_sub = df_sub[wh]
        return (df_sub, window)

    def plot_gene_tx(
        self,
        df: pd.DataFrame,
        window: list,
        ax: Any,
        show_genes: list[str] = [],
        exon_height: float = 0.2,
        fontsize: int = 12,
        ylabel: str = 'Gene',
        show_x_ticks: bool = True,
        cmap: str = 'pastel',
    ) -> None:
        """Draw a gene/transcript track onto an axis, one row per non-overlapping group.

        Transcripts are packed into display rows via
        `_find_non_overlapping_groups` so overlapping transcripts are drawn on
        separate rows, and each transcript is drawn as a line with exon
        blocks (from BED12 blockSizes/blockStarts) overlaid as rectangles.

        Args:
            df: BED12-style DataFrame of transcripts to draw (e.g. `df_sub`
                from `subset_gene`).
            window: 3-element list `[chrom, start, end]` defining the x-axis
                range to plot.
            ax: Matplotlib axis to draw on.
            show_genes: Gene names for which to draw a text label under the
                transcript.
            exon_height: Height of exon block rectangles in data units.
            fontsize: Font size for axis labels/ticks and gene name text.
            ylabel: Label for the y-axis.
            show_x_ticks: If True, show start/middle(chrom)/end x tick labels;
                otherwise hide x ticks.
            cmap: Seaborn color palette name used to color strands.
        """
        cmap = sns.color_palette(cmap)
        L = df[['chromStart', 'chromEnd', 'transcriptName']].values
        G = self._find_non_overlapping_groups(L)

        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        ax.set_xlim(x_min, x_max)

        if show_x_ticks:
            ax.set_xticks([x_min, int((x_min + x_max) / 2), x_max])
            ax.set_xticklabels([int(x_min), chrom, int(x_max)], fontsize=fontsize)
        else:
            ax.set_xticks([])

        ax.set_yticks([])
        ax.set_ylim(0, len(set(G.values())) + 1)
        ax.set_ylabel(ylabel, fontsize=fontsize)

        for n in range(df.shape[0]):
            start = df['chromStart'].values[n]
            end = df['chromEnd'].values[n]
            strand = df['strand'].values[n]
            strand = 0 if strand == '+' else 1
            tx_name = df['transcriptName'].values[n]
            gene_name = df['geneName'].values[n]
            blockCount = df['blockCount'].values[n]
            blockSizes = df['blockSizes'].values[n].split(',')
            blockStarts = df['blockStarts'].values[n].split(',')
            color = cmap[strand]
            y = G[tx_name]
            ax.plot([start, end], [y, y], color=color)
            for b in range(blockCount):
                block_start = start + int(blockStarts[b])
                block_size = int(blockSizes[b])
                ax.add_patch(Rectangle((block_start, y - exon_height / 2), block_size, exon_height, color=color))
            if gene_name in show_genes:
                ax.text((start + end) / 2, y - 0.15, gene_name, va='top', ha='center', fontsize=fontsize)

    def _find_non_overlapping_groups(
        self,
        ranges: list = [(1, 3), (2, 5), (6, 8), (9, 10), (4, 7)],
    ) -> dict:
        """Pack overlapping ranges into non-overlapping display groups (rows).

        Implements interval-graph-coloring-style packing: ranges are sorted by
        start (then end) and greedily placed into the first existing group
        whose last range ends before the current range starts, otherwise a
        new group is created. Used to assign transcripts to display rows so
        that overlapping transcripts never share a row.

        Args:
            ranges: Sequence of tuples where the first two elements are the
                start and end coordinates and the last element is used as the
                key in the returned mapping (e.g. a transcript name).

        Returns:
            Mapping from each range's last element (e.g. transcript name) to
            its 1-indexed group/row number.
        """
        # Sort the ranges by starting point (and by ending point in case of ties)
        ranges = sorted(ranges, key=lambda x: (x[0], x[1]))
        groups = []
        # Iterate through each range
        for range in ranges:
            placed = False
            # Try to place the current range in an existing group
            for group in groups:
                # Check if the current range overlaps with the last range in the current group
                if group[-1][1] < range[0]:
                    group.append(range)
                    placed = True
                    break
            # If the range couldn't be placed in any existing group, create a new group
            if not placed:
                groups.append([range])
        g = {}
        for i, group in enumerate(groups):
            for x in group:
                g[x[2]] = i + 1
        return (g)


class LocusZoomPlot(GeneTxPlot):
    """LocusZoom-style association scatter plots with LD coloring.

    Inherits gene/transcript track plotting and gene-window subsetting from
    `GeneTxPlot`. Adds a scatter-plot renderer for -log10(p) association
    signals and a tabix-based loader that pulls variant association records
    for a gene window and optionally annotates them with PLINK-computed LD
    (R2) relative to a lead variant.
    """

    def __init__(self) -> None:
        """Initialize the plotter (no state is stored)."""
        pass

    def scatter_plot(
        self,
        df: pd.DataFrame,
        window: list,
        ax: Any,
        x: str = 'pos',
        y: str = 'pv',
        s: int = 10,
        lw: float = 0,
        color: str = 'C0',
        cmap: str = 'flare',
        ylabel: str = '-log10(p)',
        hue: str | None = 'R2_bin',
        hue_title: str = 'R2',
        show_legend: bool = False,
        legend_size: int = 6,
    ) -> None:
        """Draw a LocusZoom-style association scatter plot on an axis.

        Args:
            df: DataFrame of variant records with `x`/`y` (and `hue`) columns.
            window: 3-element list `[chrom, start, end]` defining the x-axis
                range to plot.
            ax: Matplotlib axis to draw on.
            x: Column name to use for the x-axis (typically position).
            y: Column name to use for the y-axis (typically -log10(p)).
            s: Marker size.
            lw: Marker edge line width.
            color: Marker color used when `hue` is falsy.
            cmap: Palette name used to color points when `hue` is set.
            ylabel: Label for the y-axis.
            hue: Column name to color points by (e.g. LD bin), or None/empty
                to use a single `color`.
            hue_title: Legend title when `hue` is used and the legend is shown.
            show_legend: If True, show the (reversed-order) hue legend;
                otherwise remove any legend.
            legend_size: Font size for legend entries.
        """
        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        ax.set_xlim(x_min, x_max)

        if hue:
            sns.scatterplot(x=x, y=y, s=s, data=df, ax=ax, hue=hue, lw=0, palette=cmap)
        else:
            sns.scatterplot(x=x, y=y, s=s, data=df, ax=ax, color=color, lw=0)

        ax.set_xticks([])
        ax.set_xlabel('')
        ax.set_ylim(0, max(ax.get_ylim()[1], 10))
        ax.set_ylabel(ylabel)
        if show_legend:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], title=hue_title, loc='upper right', prop={'size': legend_size})
        else:
            if ax.get_legend():
                ax.get_legend().remove()

    def bgzip_to_df(
        self,
        bgzip_file: str,
        window: list,
        gene: str = 'PTGFRN',
        extra_filter: str | None = None,
        adding_ld: bool = True,
        var_id: str = 'rs1127215',
        bfile: str = 'eQTL_genotyping_sampleRenamed_rsID_variantFiltered',
        pos_col: str = 'var_from',
        pv_col: str = 'nom_pval',
        pv_min: float = 1e-300,
    ) -> pd.DataFrame:
        """Tabix-query a bgzipped nominal association file and optionally add LD.

        Queries `bgzip_file` for records within `window`, optionally restricts
        to a single gene's peak/feature (disambiguating multiple peaks per
        gene via `extra_filter` or by lowest p-value), and optionally computes
        LD (R2) to `var_id` via an external PLINK `--r2` run, binning R2 into
        `R2_bin`. Also derives numeric `pos` and -log10(p) `pv` columns and
        writes the result to `{gene}_data.txt`.

        Args:
            bgzip_file: Path to the bgzip+tabix-indexed nominal association file.
            window: 3-element list `[chrom, start, end]` to query.
            gene: Gene symbol used to filter `phe_id` and to name the output file.
            extra_filter: Substring used to pick among multiple peaks/features
                for the same gene; if None, the peak with the smallest p-value
                is used instead.
            adding_ld: If True, compute LD (R2) to `var_id` using PLINK.
            var_id: Lead variant ID used as the PLINK `--ld-snp` for LD
                calculation.
            bfile: PLINK bfile prefix (without `.bed`) used for LD calculation.
            pos_col: Column name holding variant position.
            pv_col: Column name holding the nominal p-value.
            pv_min: Minimum p-value used to clip before -log10 transform.

        Returns:
            DataFrame of queried variant records annotated with `pos`, `pv`,
            and (if `adding_ld`) `ld_var_id`, `R2`, `R2_bin`.

        Raises:
            FileNotFoundError: If `{bfile}.bed` does not exist and `adding_ld`
                is True.
            ValueError: If `adding_ld` is True and `var_id` is falsy.
        """
        header = pd.read_table(bgzip_file, header=0, nrows=0, sep='\t')
        tb = tabix.open(bgzip_file)
        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        res = tb.query(chrom, x_min, x_max)
        df = pd.DataFrame(res)
        df.iloc[:, -1] = df.iloc[:, -1].str.replace('\r', '')
        df.columns = header.columns
        if gene:
            wh = [True if gene in x.split('_')[-1].split(',') else False for x in df['phe_id']]
            df = df[wh]

        if gene and len(df['phe_id'].unique()) > 1:
            if extra_filter is None:
                # for caQTL, if there are multiple peaks of the same gene, just take the one with the smallest p-value if no extra_filter is provided
                L = []
                print('Multiple peaks of the same gene found:')
                for gi, g in df.groupby('phe_id'):
                    print([gi, g.shape[0]])
                    L.append(g)
                L = sorted(L, key=lambda x: x[pv_col].astype(float).min())
                df = L[0]
            else:
                # for caQTL, if there are multiple peaks of the same gene, take the one with the peak name containing extra_filter
                for gi, g in df.groupby('phe_id'):
                    if gi.find(extra_filter) != -1:
                        df = g
                        break

        if adding_ld:
            if not os.path.exists(bfile + '.bed'):
                raise FileNotFoundError(f'{bfile}.bed not found. Please provide a PLINK bed file for LD calculation.')
            if not var_id:
                raise ValueError('var_id is required for LD calculation. Please provide the variant ID for the lead SNP.')
            D = {}
            try:
                out_file = f'{bfile}_{var_id}'
                cmd = f'plink --bfile {bfile} --r2 --ld-snp {var_id} --ld-window-kb 2000 --ld-window 99999 --ld-window-r2 0.0 --out {out_file}'
                subprocess.run(cmd, shell=True)
                df_ld = pd.read_table(out_file + '.ld', header=0, sep=r'\s+')
                D = dict(zip(df_ld['SNP_B'], df_ld['R2']))
                D[var_id] = 1.0
                df['ld_var_id'] = var_id
                df['R2'] = df['var_id'].map(D).fillna(0)
                df['R2_bin'] = pd.cut(df['R2'], bins=5, labels=[f'{x:.1f}' for x in np.linspace(0.2, 1, 5)])
                print(df['R2_bin'].value_counts())
            except Exception as e:
                df['ld_var_id'] = var_id
                df['R2'] = 0
                df['R2_bin'] = 0
                print(f'Error in LD calculation: {e}')
        df['pos'] = pd.to_numeric(df[pos_col], errors='coerce')
        df['pv'] = pd.to_numeric(df[pv_col], errors='coerce')
        df = df.dropna(subset=['pos', 'pv'])
        df['pv'] = df['pv'].clip(lower=pv_min)
        df['pv'] = -np.log10(df['pv'])
        df.to_csv(f'{gene}_data.txt', sep='\t', index=None)
        return df


class ManhattanPlot():
    """Genome-wide Manhattan plot renderer.

    Loads chromosome sizes to build cumulative genomic coordinates, converts
    per-variant positions into those cumulative coordinates, and renders
    genome-wide -log10(p) scatter plots with alternating chromosome
    background shading.

    Attributes:
        ch_size: Mapping from chromosome label (e.g. `'1'`) to chromosome
            size in bases, in chromosome order.
        ch_cu: List of cumulative chromosome-size sums, one per chromosome,
            used as chromosome boundary positions on the cumulative x-axis.
        ch_cu_middle: List of cumulative midpoint positions, one per
            chromosome, used to place chromosome tick labels.
    """

    def __init__(
        self,
        chrom_size_file: str = 'Homo_sapiens.GRCh38.dna.primary_assembly.OneLine.ChrSizeChr.txt',
        chroms: list[str] = [str(x) for x in range(1, 23)],
    ) -> None:
        """Load chromosome sizes and precompute cumulative genomic coordinates.

        Args:
            chrom_size_file: Path to a two-column (chrom, size) chromosome
                size file.
            chroms: Chromosome labels (without the `chr` prefix) to include,
                in the desired plotting order.

        Raises:
            FileNotFoundError: If `chrom_size_file` does not exist.
        """
        if not os.path.exists(chrom_size_file):
            raise FileNotFoundError(f'{chrom_size_file} not found')
        df = pd.read_table(chrom_size_file, header=None, sep='\t')
        df.columns = ['chrom', 'size']
        df['ch'] = [x.split('chr')[-1] for x in df['chrom']]
        df = df[df['ch'].isin(chroms)]
        self.ch_size = dict(zip(df['ch'], df['size']))
        ch_size_values = list(self.ch_size.values())
        self.ch_cu = [sum(ch_size_values[0:n]) for n in range(1, len(ch_size_values) + 1)]
        self.ch_cu_middle = [sum(ch_size_values[0:n]) + ch_size_values[n] / 2 for n in range(0, len(ch_size_values))]

    def scatter_plot(
        self,
        ax: Any,
        df: pd.DataFrame,
        x: str = 'pos_cu',
        y: str = 'pv',
        s: int = 10,
        lw: float = 0,
        hue: str = 'significant',
        palette: list[str] = ['C0', 'C1'],
        ylabel: str = '-log10(p)',
        sig_threshold: float = 5e-8,
        line_params: str = 'r--',
    ) -> None:
        """Draw a genome-wide Manhattan scatter plot with a significance line.

        Args:
            ax: Matplotlib axis to draw on.
            df: DataFrame of variant records with `x`/`y` columns; a
                `significant` column is added in place.
            x: Column name to use for the x-axis (typically cumulative position).
            y: Column name to use for the y-axis (typically -log10(p)).
            s: Marker size.
            lw: Marker edge line width.
            hue: Column name to color points by (set to `'significant'` here).
            palette: Two-color palette for non-significant/significant points.
            ylabel: Label for the y-axis.
            sig_threshold: P-value threshold used to draw the significance line
                and flag `significant` points.
            line_params: Matplotlib format string for the significance line.
        """
        sig_th = -np.log10(sig_threshold)
        df['significant'] = df[y].astype(float) > sig_th
        sns.scatterplot(x=x, y=y, data=df, s=s, lw=lw, hue=hue, palette=palette, ax=ax, legend=False)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        ax.plot([0, df[x].max()], [sig_th, sig_th], line_params)

    def bgzip_to_df(
        self,
        bgzip_file: str,
        chrom_col: str = 'var_chr',
        pos_col: str = 'var_from',
        pv_col: str = 'adj_beta_pval',
        pv_min: float = 1e-300,
        unique_variants: bool = True,
        var_id: str = 'var_id',
    ) -> pd.DataFrame:
        """Load a genome-wide association file and compute cumulative positions.

        Reads `bgzip_file`, restricts to chromosomes known to `self.ch_size`,
        computes -log10(p) into `pv`, optionally deduplicates to one row per
        variant (keeping the smallest p-value), and computes a cumulative
        genomic position `pos_cu` for Manhattan-plot rendering.

        Args:
            bgzip_file: Path to the (optionally bgzipped) association file.
            chrom_col: Column name holding the chromosome.
            pos_col: Column name holding the variant position.
            pv_col: Column name holding the p-value.
            pv_min: Minimum p-value used to clip before -log10 transform.
            unique_variants: If True, keep only the best (lowest p-value) row
                per `var_id`.
            var_id: Column name holding the variant ID, used for deduplication.

        Returns:
            DataFrame with added `ch`, `pv`, and `pos_cu` columns.

        Raises:
            FileNotFoundError: If `bgzip_file` does not exist.
        """
        if not os.path.exists(bgzip_file):
            raise FileNotFoundError(f'{bgzip_file} not found')

        df = pd.read_table(bgzip_file, header=0, sep='\t')
        df['ch'] = [str(x).split('chr')[-1] for x in df[chrom_col]]
        df = df[df['ch'].isin(self.ch_size)]
        df['pv'] = -np.log10(df[pv_col].astype(float).clip(lower=pv_min))
        if unique_variants:
            df.sort_values([var_id, pv_col], inplace=True)
            df.drop_duplicates(subset=var_id, keep='first', inplace=True)

        L = []
        for n in range(df.shape[0]):
            ch = df['ch'].iloc[n]
            pos = df[pos_col].iloc[n]
            ch_size_keys = list(self.ch_size.keys())
            pos_cu = sum([self.ch_size[k] if ch_size_keys.index(k) < ch_size_keys.index(ch) else 0 for k in ch_size_keys]) + pos
            L.append(pos_cu)
        df['pos_cu'] = L
        return df

    def plot_chroms(
        self,
        ax: Any,
        show_xticklabels: bool = True,
        xticklabels_masked: list[str] = ['17', '19', '21'],
        color_bg: str = 'grey',
        alpha_bg: float = 0.2,
    ) -> None:
        """Set up chromosome tick labels and alternating background shading.

        Args:
            ax: Matplotlib axis to configure.
            show_xticklabels: If True, show chromosome labels (masking labels
                listed in `xticklabels_masked`); otherwise hide all labels.
            xticklabels_masked: Chromosome labels to blank out to reduce
                crowding, when `show_xticklabels` is True.
            color_bg: Background shading color for alternating chromosomes.
            alpha_bg: Alpha (opacity) for the background shading.
        """
        ax.set_xlim(0, self.ch_cu[-1])
        ax.set_xticks(self.ch_cu_middle)
        if show_xticklabels:
            ax.set_xticklabels([x if x not in xticklabels_masked else '' for x in self.ch_size])
        else:
            ax.set_xticklabels([])
        [ax.axvspan(self.ch_cu[n], self.ch_cu[n + 1], facecolor=color_bg, alpha=alpha_bg) for n in range(0, len(self.ch_cu) - 1, 2)]

    def get_top_signals(self, df: pd.DataFrame, pv_col: str = 'pv', top_n: int = 20) -> pd.DataFrame:
        """Return the top-N rows by descending value in `pv_col`.

        Args:
            df: DataFrame of variant records.
            pv_col: Column name to sort by (descending), typically -log10(p).
            top_n: Number of top rows to return.

        Returns:
            DataFrame containing the top `top_n` rows by `pv_col`.
        """
        df_sorted = df.sort_values(pv_col, ascending=False)
        df_top = df_sorted.head(top_n)
        return df_top


class GenoPhenoBarPlot(GeneTxPlot):
    """Genotype-vs-phenotype bar plot builder for a single variant/gene pair.

    Inherits gene-window subsetting and BED12 reading from `GeneTxPlot`. Adds
    tabix-based loaders for a phenotype BED (expression/peak/protein
    quantifications) and a genotype VCF, and a bar-plot renderer that groups
    phenotype values by genotype (0/1/2 copies of the alternate allele).
    """

    def __init__(self) -> None:
        """Initialize the plotter (no state is stored)."""
        pass

    def get_pheno_df(
        self,
        bed_file: str,
        window: list,
        gene: str = 'PTGFRN',
        extra_filter: str | None = None,
        transform: str | None = None,
    ) -> pd.DataFrame:
        """Tabix-query a phenotype BED file and reshape it to one row per sample.

        Args:
            bed_file: Path to the bgzip+tabix-indexed phenotype BED file.
            window: 3-element list `[chrom, start, end]` to query.
            gene: Gene symbol used to filter the `pid` (phenotype/peak ID) column.
            extra_filter: Substring used to pick among multiple
                peaks/features for the same gene; if None, the first matching
                feature is used instead.
            transform: If `'power2'`, back-transform values via `2 ** value`.

        Returns:
            DataFrame with columns `sample`, `feature`, `value` (one row per
            sample).

        Raises:
            ValueError: If more than one phenotype row remains after
                filtering (message references `var_id`, matching existing
                behavior).
        """
        header = pd.read_table(bed_file, header=0, nrows=0, sep='\t')
        tb = tabix.open(bed_file)
        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        res = tb.query(chrom, x_min, x_max)
        df = pd.DataFrame(res)
        df.columns = header.columns
        if gene:
            wh = [True if gene in x.split('_')[-1].split(',') else False for x in df['pid']]
            df = df[wh]

        if gene and len(df['pid'].unique()) > 1:
            if extra_filter is None:
                print(df['pid'].unique())
                print(f'Multiple peaks of the same gene found. Using the first one {df["pid"].iloc[0]}')
                wh = df['pid'] == df['pid'].iloc[0]
                df = df[wh]
            else:
                print(df['pid'].unique())
                print(f'Multiple peaks of the same gene found. Using {extra_filter} for filtering')
                wh = df['pid'].str.contains(extra_filter)
                df = df[wh]

        if df.shape[0] == 1:
            dft = pd.DataFrame()
            dft['sample'] = df.columns[6:]
            dft['feature'] = df['pid'].iloc[0]
            L = []
            for n in range(6, df.shape[1]):
                value = df.iloc[0, n]
                try:
                    value = float(value)
                except:
                    value = np.nan
                L.append(value)
            dft['value'] = L
            if transform == 'power2':
                dft['value'] = np.power(2, dft['value'])
            return dft
        else:
            raise ValueError(f'Check if the var_id {var_id} is correct.')
        return dft

    def get_geno_df(
        self,
        vcf_file: str,
        window: list,
        var_id: str = 'rs1127215',
        extra_filter: str | None = None,
    ) -> pd.DataFrame:
        """Tabix-query a VCF file and reshape genotypes to one row per sample.

        Args:
            vcf_file: Path to the bgzip+tabix-indexed VCF file.
            window: 3-element list `[chrom, start, end]` to query.
            var_id: Variant ID (VCF `ID` field) to select.
            extra_filter: Alternate allele to select when multiple ALT
                alleles exist for `var_id`; if None, the first ALT allele is
                used instead.

        Returns:
            DataFrame with columns `sample`, `REF`, `ALT`, `genotype` (0/1/2
            copies of the alternate allele, or NaN if ungenotyped).

        Raises:
            ValueError: If more than one variant row remains after filtering
                (message references `var_id`).
        """
        tb = tabix.open(vcf_file)
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    break
        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        res = tb.query(chrom, x_min, x_max)
        df = pd.DataFrame(res)
        df.columns = header
        if var_id:
            wh = df['ID'] == var_id
            df = df[wh]

        if var_id and len(df['ALT'].unique()) > 1:
            if extra_filter is not None:
                df = df[df['ALT'] == extra_filter]
            else:
                print('Multiple variants found for the given var_id. Using the first one for now')
                print(df[['ID', 'ALT']].unique())
                df = df[df['ALT'] == df['ALT'].iloc[0]]

        if df.shape[0] == 1:
            dft = pd.DataFrame()
            dft['sample'] = df.columns[9:]
            dft['REF'] = df['REF'].iloc[0]
            dft['ALT'] = df['ALT'].iloc[0]
            L = []
            for n in range(9, df.shape[1]):
                gt = df.iloc[0, n].split(':')[0]
                if gt in ['0/0', '0|0']:
                    g = 0
                elif gt in ['0/1', '1/0', '0|1', '1|0']:
                    g = 1
                elif gt in ['1/1', '1|1']:
                    g = 2
                else:
                    g = np.nan
                L.append(g)
            dft['genotype'] = L
            return dft
        else:
            raise ValueError(f'Check if the var_id {var_id} is correct.')

    def bar_plot(
        self,
        df_geno: pd.DataFrame,
        df_pheno: pd.DataFrame,
        out_file: str = 'geno_pheno_barplot.pdf',
        figsize: tuple[float, float] = (4, 4),
        title: str = 'eQTL',
        label: str = 'eQTL',
        color: str = 'C0',
        capsize: float = 0.05,
    ) -> None:
        """Render and save a genotype-vs-phenotype bar plot for one variant/gene pair.

        Merges genotype and phenotype tables on `sample`, drops rows missing
        `genotype`/`value`, and draws a bar plot of phenotype value by
        genotype (0/1/2 copies of ALT), with x tick labels showing the
        genotype string and sample count per group.

        Args:
            df_geno: Genotype DataFrame as returned by `get_geno_df`.
            df_pheno: Phenotype DataFrame as returned by `get_pheno_df`.
            out_file: Path to save the resulting figure to.
            figsize: Figure size in inches, `(width, height)`.
            title: Plot title.
            label: QTL type label (e.g. `'eQTL'`, `'caQTL'`, `'pQTL'`); its
                prefix (split on `'_'`) selects the y-axis label.
            color: Bar color.
            capsize: Error bar cap size passed to `seaborn.barplot`.
        """
        df = pd.merge(df_geno, df_pheno, on='sample')
        df = df.dropna(subset=['genotype', 'value'])

        LabelName = {'eQTL':'Gene expression (TPM)', 'caQTL':'Peak counts (TPM)', 'pQTL':'Protein abundance (intensity)'}
        ValueCounts = pd.DataFrame(df['genotype'].value_counts()).reset_index()
        ValueCounts = dict(zip(ValueCounts['genotype'], ValueCounts['count']))
        ref, alt = df_geno['REF'].iloc[0], df_geno['ALT'].iloc[0]
        GT = {0: f'{ref}/{ref}', 1: f'{ref}/{alt}', 2: f'{alt}/{alt}'}

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.barplot(x='genotype', y='value', data=df, capsize=capsize, ax=ax, color=color)

        xticks = [0, 1, 2]
        ax.set_xticks(xticks)
        xtick_labels = [f'{GT.get(x, './.')}\nN={ValueCounts.get(x, 0)}' for x in xticks]
        ax.set_xticklabels(xtick_labels)
        ax.set_xlabel('')
        ylabel = LabelName.get(label.split('_')[0], 'Phenotype value')
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        plt.tight_layout()
        plt.savefig(out_file)


class GenomeTracks():
    """pyGenomeTracks configuration builder and runner.

    Builds a `.ini` track configuration from a YAML track-config file (only
    including tracks whose backing files exist), can export a single gene's
    canonical model plus a QTL result set as bed tracks, and invokes the
    external `pyGenomeTracks` command to render a track figure.
    """

    def __init__(self) -> None:
        """Initialize an empty track-config accumulator."""
        self.config_list = []

    def subset_gene(
        self,
        gene: str = 'GINS4',
        bed12: str = 'Homo_sapiens.GRCh38.115.bed12',
        out_file: str | None = None,
        canonical_only: bool = True,
        flank: float = 1e4,
        geneBiotype_include: list[str] = ['protein_coding'],
        geneName_exclude: list[str] = ['.'],
    ) -> list:
        """Subset a BED12 file to one gene's window and write it as a bed track.

        Args:
            gene: Gene symbol (matched against `geneName`) to center the window on.
            bed12: Path to the BED12 gene-model annotation file.
            out_file: Path to write the windowed bed subset to; defaults to
                `gene_{gene}.bed`.
            canonical_only: If True, keep only rows whose `tag` field contains
                'canonical'.
            flank: Number of bases to extend on each side of the gene's
                min/max coordinates to build the display window.
            geneBiotype_include: Gene biotypes to keep in the windowed subset.
            geneName_exclude: Gene names to drop from the windowed subset.

        Returns:
            3-element list `[chrom, start, end]` (chrom as a `'chr1'`-style
            string, start/end as ints) describing the gene's display window.

        Raises:
            FileNotFoundError: If `bed12` does not exist.
        """
        if os.path.exists(bed12):
            df = pd.read_table(bed12, header=None, sep='\t', low_memory=False)
            df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                        'geneID', 'geneName', 'geneBiotype', 'transcriptID', 'transcriptName', 'transcriptBiotype', 'tag']
        else:
            raise FileNotFoundError(f'{bed12} not found.')

        if out_file is None:
            out_file = f'gene_{gene}.bed'

        wh = df['geneName'].isin([gene])
        df_sub = df[wh]
        if canonical_only:
            wh = np.array([True if x.find('canonical') != -1 else False for x in df_sub['tag']])
            df_sub = df_sub[wh]

        chrom = str(df_sub['chrom'].iloc[0])
        start = df_sub['chromStart'].min()
        end = df_sub['chromEnd'].max()
        window = ['chr' + chrom, int(start - flank), int(end + flank)]

        wh1 = (df['chrom'] == chrom) & (df['chromStart'] > window[1]) & (df['chromEnd'] < window[2])
        wh2 = df['geneBiotype'].isin(geneBiotype_include)
        wh3 = ~df['geneName'].isin(geneName_exclude)
        df_sub = df[wh1 & wh2 & wh3]
        df_sub['name'] = df_sub['geneName']
        if canonical_only:
            wh = np.array([True if x.find('canonical') != -1 else False for x in df_sub['tag']])
            df_sub = df_sub[wh]
        df_sub.to_csv(out_file, sep='\t', index=None, header=None)
        return window

    def add_tracks(
        self,
        gene: str = 'PTGFRN',
        config: str = 'config_PTGFRN.yaml',
        ini_file: str | None = None,
        ylim_off: bool = True,
    ) -> None:
        """Build a pyGenomeTracks `.ini` file from a YAML track configuration.

        Loads `config` (a YAML mapping of track name to track options), drops
        tracks whose `file` option points to a missing file, optionally
        strips `min_value`/`max_value` options (to let y-limits auto-scale),
        and writes the resulting `[track]` sections to `self.ini_file`.

        Args:
            gene: Gene symbol (unused for logic, kept for interface symmetry).
            config: Path to the YAML track configuration file.
            ini_file: Path to write the generated `.ini` file to; defaults to
                `config` with its extension replaced by `.ini`.
            ylim_off: If True, drop `min_value`/`max_value` options from every
                track so y-axis limits auto-scale.

        Raises:
            FileNotFoundError: If `config` does not exist.
        """
        if os.path.exists(config):
            with open(config) as f:
                self.config = yaml.safe_load(f)
                print(f'Number of tracks: {len(self.config)}')
        else:
            raise FileNotFoundError(f'{config} not found.')

        if ini_file is None:
            ini_file = config.replace('.yaml', '.ini')
        self.ini_file = ini_file

        config = {}
        for track in self.config:
            if 'file' in self.config[track] and not os.path.exists(self.config[track]['file']):
                pass
            else:
                config[track] = self.config[track]

        for track in config:
            self.config_list.append(f'[{track}]')
            for k, v in self.config[track].items():
                if v != 'None':
                    if ylim_off:
                        if k in ['min_value', 'max_value']:
                            continue
                    self.config_list.append(f'{k} = {v}')

        with open(self.ini_file, 'w') as f:
            f.write('\n'.join(self.config_list))

    def plot_tracks(
        self,
        gene: str,
        window: list,
        width: float,
        out_file: str | None = None,
        conda_env: str | None = 'pygenometracks',
    ) -> None:
        """Invoke pyGenomeTracks to render the configured tracks for a window.

        Args:
            gene: Gene symbol used in the plot title and default output filename.
            window: 3-element list `[chrom, start, end]` defining the region to plot.
            width: Figure width passed to `pyGenomeTracks --width`.
            out_file: Output figure path; defaults to `{gene}_tracks.pdf`.
            conda_env: Conda environment name to run `pyGenomeTracks` in via
                `conda run -n {conda_env}`; if falsy, run directly.
        """
        if out_file is None:
            out_file = f'{gene}_tracks.pdf'

        cmd = f'pyGenomeTracks --tracks {self.ini_file} --region {window[0]}:{window[1]}-{window[2]} --width {width} --title {gene} -o {out_file}'
        if conda_env:
            cmd = f'conda run -n {conda_env} ' + cmd
        print(cmd)
        subprocess.run(cmd, shell=True)

    def qtl_txt_to_bed(
        self,
        in_file: str = 'eQTL.txt.gz',
        window: list = [],
        gene: str = 'PTFRN',
        p_col: str = 'nom_pval',
        p_min: float = 1e-300,
        cols: list[str] = ['var_chr', 'var_from', 'var_to', 'pv'],
    ) -> None:
        """Tabix-query a QTL nominal file and export a per-gene bed track.

        Queries `in_file` for records within `window`, restricts to `gene`,
        deduplicates by `var_id` keeping the smallest p-value, adds a
        -log10(p) `pv` column, and writes/bgzips/tabix-indexes the result as
        `{in_file with .txt.gz -> _{gene}.bed}`.

        Args:
            in_file: Path to the bgzip+tabix-indexed QTL nominal association file.
            window: 3-element list `[chrom, start, end]` to query.
            gene: Gene symbol used to filter records (matched against the
                first column) and to name the output file.
            p_col: Column name holding the nominal p-value.
            p_min: Minimum p-value used to clip before -log10 transform.
            cols: Columns to write to the output bed file, in order.

        Raises:
            FileNotFoundError: If `in_file` does not exist.
        """
        if os.path.exists(in_file):
            header = pd.read_table(in_file, header=0, nrows=0, sep='\t').columns
        else:
            raise FileNotFoundError(f'{in_file} not found.')
        tb = tabix.open(in_file)
        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        res = tb.query(chrom, x_min, x_max)
        out_file = in_file.replace('.txt.gz', f'_{gene}.bed')
        if res:
            df = pd.DataFrame(res)
            df.columns = header
            wh = [True if gene in x.split('_')[-1].split(',') else False for x in df.iloc[:, 0]]
            df = df[wh]
            df.sort_values(p_col, inplace=True)
            df.drop_duplicates(subset=['var_id'], keep='first', inplace=True)
            df['pval'] = np.clip(df[p_col].astype(float), a_min=p_min, a_max=None)
            df['pv'] = -np.log10(df['pval'])
            df.sort_values(cols, inplace=True)
            df[cols].to_csv(out_file, sep='\t', index=None, header=None)
            cmd = f'bgzip -f {out_file}; tabix -f -p bed {out_file}.gz'
            subprocess.run(cmd, shell=True)
        else:
            print(f'No data found in {in_file} for the given window and gene.')


if __name__ == '__main__':
    def plot_locus(gene: str = 'PTGFRN', var_id: str = 'rs1127215') -> None:
        """Render a 4-panel gene-track + caQTL/eQTL/pQTL locus zoom figure."""
        fig = plt.figure()
        ax1 = fig.add_axes([0.15, 0.05, 0.75, 0.2])
        ax2 = fig.add_axes([0.15, 0.27, 0.75, 0.2])
        ax3 = fig.add_axes([0.15, 0.49, 0.75, 0.2])
        ax4 = fig.add_axes([0.15, 0.71, 0.75, 0.2])

        lzp = LocusZoomPlot()

        df = lzp.read_bed12(bed12='Homo_sapiens.GRCh38.115.bed12')
        dfg, window = lzp.subset_gene(df=df, gene=gene)
        lzp.plot_gene_tx(dfg, window=window, ax=ax1, show_genes=[gene])

        dfs = lzp.bgzip_to_df(bgzip_file='caQTL_nominal-1.0_w1k_qvalue_extraInfo.txt.gz', window=window, gene=gene, var_id=var_id, bfile='caQTL_genotyping_sampleRenamed_rsID_variantFiltered')
        lzp.scatter_plot(dfs, window=window, ax=ax4, ylabel='caQTL', show_legend=True)

        dfs = lzp.bgzip_to_df(bgzip_file='eQTL_nominal-1.0_w1M_PC25_extraInfo.txt.gz', window=window, gene=gene, var_id=var_id, bfile='eQTL_genotyping_sampleRenamed_rsID_variantFiltered')
        lzp.scatter_plot(dfs, window=window, ax=ax3, ylabel='-log10(p)\neQTL')

        dfs = lzp.bgzip_to_df(bgzip_file='pQTL_nominal-1.0_w1M_PC25_extraInfo.txt.gz', window=window, gene=gene, var_id=var_id, bfile='pQTL_genotyping_sampleRenamed_rsID_variantFiltered')
        lzp.scatter_plot(dfs, window=window, ax=ax2, ylabel='pQTL')

        ax4.set_title(f'{gene} {var_id}')

        plt.savefig(f'{gene}.pdf')

    def geno_pheo_bar_plot(
        geno_file: str,
        pheno_file: str,
        gene: str,
        var_id: str,
        label: str,
        transform: str | None = None,
        color: str = 'C0',
        extra_filter: str | None = None,
        alt_filter: str | None = None,
    ) -> None:
        """Build and save a genotype-vs-phenotype bar plot for one variant/gene pair.

        Note: this demo function contains pre-existing unreachable code after
        its early logic (it references `gpb` and `df` before they are defined
        later in this `__main__` block); left as-is intentionally.
        """
        dfg, window = gpb.subset_gene(df, gene)

        df_pheno = gpb.get_pheno_df(bed_file=pheno_file, window=window, gene=gene, transform=transform, extra_filter=extra_filter)
        df_geno = gpb.get_geno_df(vcf_file=geno_file, window=window, var_id=var_id, extra_filter=alt_filter)

        gpb.bar_plot(df_geno=df_geno, df_pheno=df_pheno, out_file=f'{gene}_{var_id}_{label}_barplot.pdf', title=f'{gene} {var_id}', label=label, color=color)

        ## locus zoom plot
        plot_locus(gene='PTGFRN', var_id='rs1127215')
        plot_locus(gene='PEPD', var_id='rs79910652')
        plot_locus(gene='STARD10', var_id='rs140130268')

        ## geno pheno bar plot
        gpb = GenoPhenoBarPlot()
        bed12 = 'Homo_sapiens.GRCh38.115.bed12'
        df = gpb.read_bed12(bed12)
        cmap = sns.color_palette('Dark2')

        eqtl_vcf = 'eQTL_genotyping_sampleRenamed_rsID_variantFiltered.vcf.gz'
        eqtl_bed = 'eQTL_geneCounts_geneName_TPM_geneFiltered.bed.gz'
        caqtl_vcf = 'caQTL_genotyping_sampleRenamed_rsID_variantFiltered.vcf.gz'
        caqtl_bed = 'ATACseq_qvalue_peakCounts_closestGene_TPM_peakFiltered.bed.gz'
        pqtl_vcf = 'pQTL_genotyping_sampleRenamed_rsID_variantFiltered.vcf.gz'
        pqtl_bed = 'Proteomics_subsetRenamed_proteinFiltered.bed.gz'

        geno_pheo_bar_plot(geno_file=eqtl_vcf, pheno_file=eqtl_bed, gene='PTGFRN', var_id='rs1127215', label='eQTL', color=cmap[0])
        geno_pheo_bar_plot(geno_file=caqtl_vcf, pheno_file=caqtl_bed, gene='PTGFRN', var_id='rs1127215', label='caQTL_peak_chr1_116988204_116989436', color=cmap[1], extra_filter='chr1_116988204_116989436_PTGFRN')
        geno_pheo_bar_plot(geno_file=pqtl_vcf, pheno_file=pqtl_bed, gene='PTGFRN', var_id='rs1127215', label='pQTL', transform='power2', color=cmap[2])


class GenomeBrowser():
    """UCSC-track-hub and bed-track export helpers.

    Builds a JSON track-hub configuration from a tab-delimited track config
    table, and converts GWAS summary-statistics / QTL nominal-association
    files into sorted, bgzipped, tabix-indexed bed tracks of -log10(p)
    values for genome-browser display.
    """

    def __init__(self) -> None:
        """Initialize the browser helper (no state is stored)."""
        pass

    def get_hub_config(
        self,
        in_file: str = 'track_config.txt',
        out_file: str = 'hub.config.json',
    ) -> None:
        """Build a JSON track-hub configuration from a tab-delimited track table.

        Reads `in_file` (one row per track, with a `url` column and optional
        `options_*` columns), skips rows whose `url` file does not exist, and
        writes the remaining rows as a JSON list of track dicts (with
        `options_*` columns nested under an `options` dict) to `out_file`.

        Args:
            in_file: Path to the tab-delimited track configuration table.
            out_file: Path to write the resulting JSON track-hub config to.
        """
        df = pd.read_table(in_file, sep='\t', header=0, dtype=str)
        L = []
        for n in range(df.shape[0]):
            file = df['url'].iloc[n]
            if os.path.exists(file):
                D = {}
                D.setdefault('options', {})
                for col in df.columns:
                    if col.startswith('options_'):
                        col2 = col.replace('options_', '')
                        val = df[col].iloc[n]
                        if val is not np.nan:
                            if col2 in ['height']:
                                val = int(val)
                            elif col2 in ['category']:
                                val = json.loads(val)
                            D['options'][col2] = val
                    else:
                        val = df[col].iloc[n]
                        D[col] = val
                L.append(D)
            else:
                print(f'File {file} not found.')

        with open(out_file, 'w') as f:
            json.dump(L, f, indent=4)

    def gwas_to_bed(
        self,
        in_file: str = 'T2D_GGI_EUR_sumstat_harmoniser.h.tsv.gz',
        out_file: str = 'T2D_GGI_EUR.bed.gz',
        p_param: float = 1e-300,
    ) -> None:
        """Convert a gzipped GWAS summary-stats file into a sorted, tabix-indexed bed track.

        Reads `in_file` line by line, extracting chromosome (column 0),
        position (column 1), and p-value (column 7), clipping p-values to
        `p_param` and writing `-log10(p)` as a bed4 record; lines that fail
        to parse are silently skipped. The output is sorted, bgzipped, and
        tabix-indexed.

        Args:
            in_file: Path to the gzipped GWAS summary-statistics TSV file.
            out_file: Output bed.gz path (an intermediate uncompressed `.bed`
                is written and then bgzipped/removed).
            p_param: Minimum p-value used to clip before -log10 transform.
        """
        out_file = out_file.replace('.bed.gz', '.bed')
        with gzip.open(in_file, 'rt') as f, open(out_file, 'w') as fo:
            head = f.readline()
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                try:
                    chrom = 'chr' + fields[0]
                    pos = int(fields[1])
                    p = float(fields[7])
                    p_min = max(p, p_param)
                    log_p = -1 * np.log10(p_min)
                    fo.write(f"{chrom}\t{pos-1}\t{pos}\t{log_p}\n")
                except:
                    pass

        subprocess.run(f'sort -k1,1V -k2,2n {out_file} |bgzip > {out_file}.gz', shell=True)
        subprocess.run(f'tabix -p bed {out_file}.gz; rm {out_file}', shell=True)

    def qtl_to_bed(
        self,
        in_file: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo.txt.gz',
        out_file: str = 'eQTL.bed.gz',
        p_param: float = 1e-300,
    ) -> None:
        """Convert a gzipped QTL nominal-association file into a per-variant bed track.

        Reads `in_file`, keeping for each `var_id` the record with the
        smallest p-value across all its rows, then writes a sorted,
        bgzipped, tabix-indexed bed4 track of `-log10(p)` values, one row per
        variant. Lines that fail to parse are silently skipped.

        Args:
            in_file: Path to the gzipped QTL nominal association file.
            out_file: Output bed.gz path (an intermediate uncompressed `.bed`
                is written and then bgzipped/removed).
            p_param: Minimum p-value used to clip before -log10 transform.
        """
        out_file = out_file.replace('.bed.gz', '.bed')
        D = {}
        with gzip.open(in_file, 'rt') as f:
            head = f.readline().strip().split('\t')
            idx_var_id = head.index('var_id')
            idx_chrom = head.index('var_chr')
            idx_pos = head.index('var_from')
            idx_p = head.index('nom_pval')
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                try:
                    var_id = fields[idx_var_id]
                    chrom = fields[idx_chrom]
                    pos = int(fields[idx_pos])
                    p = float(fields[idx_p])
                    D.setdefault(var_id, [])
                    D[var_id].append([chrom, pos, p])
                except:
                    pass

        with open(out_file, 'w') as fo:
            for var_id in sorted(D):
                L = sorted(D[var_id], key=lambda x: x[2])
                chrom, pos, p = L[0]
                p_min = max(p, p_param)
                log_p = -1 * np.log10(p_min)
                fo.write(f"{chrom}\t{pos-1}\t{pos}\t{log_p}\n")

        subprocess.run(f'sort -k1,1V -k2,2n {out_file} |bgzip > {out_file}.gz', shell=True)
        subprocess.run(f'tabix -p bed {out_file}.gz; rm {out_file}', shell=True)
