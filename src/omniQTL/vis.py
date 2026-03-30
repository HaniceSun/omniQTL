from .utils import *
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']
matplotlib.rcParams['savefig.dpi'] = 300
import matplotlib.patches as mpatches
import pylab as plt
import seaborn as sns

class Visualization:
    def __init__(self):
        pass

    def get_table_for_qq_plot(self, in_file='pQTL_permute-1000_w1M_PC25_extraInfo.txt.gz', out_file='pQTL_qq_plot_table.txt', p_col='adj_beta_pval'):
        L = []
        df = pd.read_table(in_file, header=0, sep='\t')
        for gi, g in df.groupby('phe_id'):
            g_sub = g.dropna(subset=[p_col])
            if g_sub.shape[0]:
                L.append([gi, g_sub[p_col].min()])
        df = pd.DataFrame(L, columns=['phe_id', 'min_p'])
        df.sort_values('min_p', inplace=True)
        n = df.shape[0]
        df['expected'] = -np.log10(np.arange(1, n+1) / (n+1))
        df['observed'] = -np.log10(df['min_p'])
        df.to_csv(out_file, index=False, sep='\t')

    def qq_plot(self, in_file, title='QQ Plot', markerscale=4, scatter_size=4, color='C1'):
        out_file = in_file.replace('.txt', '.pdf')
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure()
        ax = fig.add_subplot()
        sns.scatterplot(x='expected', y='observed', s=scatter_size, ax=ax, data=df, color=color)
        ax.plot([0, df['expected'].max()], [0, df['expected'].max()], linestyle="--", color='C0')
    
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

    def get_number_raw_peaks(self, in_dirs, out_file='caQTL_number_raw_peaks.txt'):
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

    def plot_number_raw_peaks(self, in_file='caQTL_number_raw_peaks.txt', out_file='caQTL_number_raw_peaks_violin.pdf', ylabel='Number of raw peaks'):
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure()
        ax = fig.add_subplot()
        sns.violinplot(x='peak_type', y='number_of_peaks', data=df, ax=ax)
        sns.stripplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, color='C1', size=4)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        plt.savefig(out_file)

    def get_number_independent_signals(self, in_files=[], out_file='QTL_number_of_independent_signals.txt'):
        L = []
        for f in in_files:
            df = pd.read_table(f, header=0, sep='\t')
            counts = df['n_independent_signals'].value_counts().to_frame().reset_index()
            counts['qtl_type'] = f.split('_')[0]
            L.append(counts)
        df = pd.concat(L, axis=0)
        df.to_csv(out_file, index=False, sep='\t')
    
    def plot_number_independent_signals(self, in_file='QTL_number_of_independent_signals.txt', show_numbers=True, ylim=[0, 10000], title='QTL conditional analysis', cmap='Dark2'):
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
