from .utils import *
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']
matplotlib.rcParams['savefig.dpi'] = 300
import pylab as plt
import seaborn as sns

class Visulization:
    def __init__(self):
        pass

    def get_table_for_qq_plot(self, in_files=['pQTL_permute-1000_w1M_PC25_extraInfo.txt.gz'], out_file='QTL_qq_plot_table.txt', p_col='adj_beta_pval'):
        dfs = []
        for f in in_files:
            qtl = f.split('_')[0]
            L = []
            df = pd.read_table(f, header=0, sep='\t')
            for gi, g in df.groupby('phe_id'):
                g_sub = g.dropna(subset=[p_col])
                if g_sub.shape[0]:
                    L.append([gi, g_sub[p_col].min(), qtl])
            df = pd.DataFrame(L, columns=['phe_id', 'min_p', 'qtl'])
            df.sort_values('min_p', inplace=True)
            n = df.shape[0]
            df['expected'] = -np.log10(np.arange(1, n+1) / (n+1))
            df['observed'] = -np.log10(df['min_p'])
            dfs.append(df)
    
        df = pd.concat(dfs)
        df.to_csv(out_file, index=False, sep='\t')

    def qq_plot(self, in_file, title='QQ Plot', markerscale=4, scatter_size=4, cmap='deep'):
        out_file = in_file.replace('.txt', '.pdf')
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure()
        ax = fig.add_subplot()
        sns.scatterplot(x='expected', y='observed', s=scatter_size, ax=ax, data=df, hue='qtl', palette=cmap)
        for gi, g in df.groupby('qtl'):
            expected = g['expected']
            observed = g['observed']
            ax.plot([0, expected.max()], [0, expected.max()], linestyle="--")
    
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        ax.set_title(title)
        ax.legend(title=None, markerscale=markerscale)
        plt.savefig(out_file)

    def get_table_for_upset_plot(self, in_files=['caQTL_permute-1000_w1k_qvalue.significant.txt'], out_file='QTL_upset_plot_table.txt'):
        L = []
        for f in in_files:
            qtl = f.split('_')[0]
            genes = []
            df = pd.read_table(f, header=0, sep=r'\s+')
            if qtl in ['caQTL']:
                for item in df.iloc[:, 0]:
                    x = item.split('_')[-1].split(',')
                    genes += x
            else:
                for item in df.iloc[:, 0]:
                    genes.append(item.split('_')[-1])
            genes_uniq= sorted(set(genes))
            L.append([qtl, len(genes_uniq), ','.join(genes_uniq)])
        df = pd.DataFrame(L, columns=['qtl', 'numer_of_genes', 'genes'])
        df.to_csv(out_file, header=True, index=False, sep='\t')

    def upset_plot(self, in_file='QTL_upset_plot_table.txt', cmap='deep'):
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
