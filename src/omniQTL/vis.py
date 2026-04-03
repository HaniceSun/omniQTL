from .utils import *
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

class GeneTxPlot():
    def __init__(self):
        pass

    def read_bed12(self, bed12='Homo_sapiens.GRCh38.115.bed12'):
        df = pd.read_table(bed12, header=None, sep='\t', low_memory=False)
        df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts',
                        'geneID', 'geneName', 'geneBiotype', 'transcriptID', 'transcriptName', 'transcriptBiotype', 'tag']
        return df

    def subset_gene(self, df, gene='GINS4', canonical_only=True, flank=1e6, geneBiotype_include=['protein_coding'], geneName_exclude=['.']):
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

    def plot_gene_tx(self, df, window, ax, show_genes=[], exon_height=0.2, fontsize=12, ylabel='Gene', show_x_ticks=True, cmap='pastel'):
        cmap = sns.color_palette(cmap)
        L = df[['chromStart', 'chromEnd', 'transcriptName']].values
        G = self._find_non_overlapping_groups(L)

        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        ax.set_xlim(x_min, x_max)

        if show_x_ticks:
            ax.set_xticks([x_min, int((x_min + x_max)/2), x_max])
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
                ax.add_patch(Rectangle((block_start, y - exon_height/2), block_size, exon_height, color=color))
            if gene_name in show_genes:
                ax.text((start + end)/2, y - 0.15, gene_name, va='top', ha='center', fontsize=fontsize)

    def _find_non_overlapping_groups(self, ranges=[(1, 3), (2, 5), (6, 8), (9, 10), (4, 7)]):
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
        return(g)

class LocusZoomPlot(GeneTxPlot):
    def __init__(self):
        pass

    def scatter_plot(self, df, window, ax, x='pos', y='pv', s=10, lw=0, color='C0', cmap='flare', ylabel='-log10(p)', hue='R2_bin', hue_title='R2', show_legend=False, legend_size=6):
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


    def bgzip_to_df(self, bgzip_file, window, gene='PTGFRN', extra_filter=None, adding_ld=True, var_id='rs1127215', bfile='eQTL_genotyping_sampleRenamed_rsID_variantFiltered', pos_col='var_from', pv_col='nom_pval', pv_min=1e-300):
        header = pd.read_table(bgzip_file, header=0, nrows=0, sep='\t')
        tb = tabix.open(bgzip_file)
        chrom = window[0]
        x_min = window[1]
        x_max = window[2]
        res = tb.query(chrom, x_min, x_max)
        df = pd.DataFrame(res)
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

if __name__ == '__main__':
    pass

