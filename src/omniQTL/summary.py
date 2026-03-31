from .utils import *

class Summary:
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

    def qq_plot(self, in_file, title='QQ plot', markerscale=4, scatter_size=4, color='C1', figsize=(4, 4)):
        out_file = in_file.replace('.txt', '.pdf')
        qtl_type = in_file.split('_')[0]
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.scatterplot(x='expected', y='observed', s=scatter_size, ax=ax, data=df, color=color)
        ax.plot([0, df['expected'].max()], [0, df['expected'].max()], linestyle="--", color='C0')
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        ax.set_title(title + ' of ' + qtl_type)
        plt.tight_layout()
        plt.savefig(out_file)

    def bar_plot_significant_loci(self, in_file, axes=[0.3, 0.4, 0.6, 0.5], cmap='Dark2', show_numbers=True):
        out_file = in_file.split('.txt')[0] + '.pdf'
        cmap = sns.color_palette(cmap)

        df = pd.read_table(in_file, header=0, sep='\t')
        df['StudySampleSize'] = [f"{df['Study'].iloc[n]}\n(N={df['SampleSize'].iloc[n]})" for n in range(df.shape[0])]
    
        fig = plt.figure()
        ax = fig.add_axes(axes)
        sns.barplot(y='Number of significant loci', x='StudySampleSize', hue='StudySampleSize', data=df, palette=[cmap[0], cmap[-1], cmap[1], cmap[-1], cmap[2]], legend=False)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_xlabel('')
        if show_numbers:
            for i, row in df.iterrows():
                N = row['Number of significant loci']
                ax.text(i, N, N, ha='center', va='bottom')
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], ylim[1] * 1.1)

        y1 = -0.4
        y2 = -0.5
        ax.plot([0, 0, 1, 1], [y1, y2, y2, y1], transform=ax.get_xaxis_transform(), lw=1, color='k', clip_on=False)
        ax.plot([2, 2, 3, 3], [y1, y2, y2, y1], transform=ax.get_xaxis_transform(), lw=1, color='k', clip_on=False)
        ax.plot([4, 4], [y1, y2], transform=ax.get_xaxis_transform(), lw=1, color='k', clip_on=False)
        ax.text(0.5, y2-0.03, 'caQTL', ha='center', va='top', transform=ax.get_xaxis_transform())
        ax.text(2.5, y2-0.03, 'eQTL', ha='center', va='top', transform=ax.get_xaxis_transform())
        ax.text(4, y2-0.03, 'pQTL', ha='center', va='top', transform=ax.get_xaxis_transform())

        #plt.tight_layout()
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

    def plot_number_raw_peaks(self, in_file='caQTL_number_raw_peaks.txt', out_file='caQTL_number_raw_peaks_boxplot.pdf', cmap='Blues', ylabel='Number of raw peaks (million)', add_stripplot=False, figsize=(4, 4)):
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if ylabel.find('million') != -1:
            df['number_of_peaks'] = df['number_of_peaks']/1e6
        sns.boxplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_type', palette=cmap, legend=False)
        if add_stripplot:
            sns.stripplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, color='C0', size=4)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_number_merged_peaks(self, in_files, out_file='caQTL_number_merged_peaks.txt', params={'consensus':'consensus peaks', 'summitExtended':'summit extended peaks'}):
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

    def plot_number_merged_peaks(self, in_file='caQTL_number_merged_peaks.txt', out_file='caQTL_number_merged_peaks_barplot.pdf', cmap='Blues', ylabel='Number of merged peaks (million)', figsize=(4, 4)):
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if ylabel.find('million') != -1:
            df['number_of_peaks'] = df['number_of_peaks']/1e6
        sns.barplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_method', palette=cmap)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        ax.legend(title=None)
        plt.tight_layout()
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

    def get_sig_variants(self, in_file='eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz'):
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

    def get_non_sig_variants(self, in_file='eQTL_nominal-1.0_w1M_PC25_extraInfo.txt.gz', params={'p_col': 'nom_pval', 'p_threshold': 0.05}):
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

    def split_Ensembl_regulatory_annotation(self, in_file='Ensembl_BioMart_RegulatoryAnnotation.txt.gz'):
        df = pd.read_table(in_file, header=0, sep='\t')
        for gi, g in df.groupby('Feature type'):
            out_file = in_file.replace('.txt.gz', f'_{gi}.bed')
            g_sub = g.loc[:, ['Chromosome/scaffold name', 'Start (bp)', 'End (bp)']]
            g_sub.columns = ['ch', 'start', 'end']
            g_sub['ch'] = [f'chr{x}' for x in g_sub['ch']]
            g_sub.to_csv(out_file, header=False, index=False, sep='\t')

    def get_nominal_sig_associations(self, in_files=['caQTL_nominal-1.0_w1k_qvalue_extraInfo_sig.txt.gz',
                                   'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz', 'pQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz'],
                                   out_file='QTL_nomnial_sig_associations.txt'):
        L = []
        for f in in_files:
            qtl_type = f.split('_')[0]
            with gzip.open(f, 'rt') as fin:
                head = fin.readline().strip().split('\t')
                phe_idx = head.index('phe_id')
                var_idx = head.index('var_id')
                p_idx = head.index('nom_pval')
                beta_idx = head.index('slope')
                for line in fin:
                    items = line.strip().split('\t')
                    if qtl_type in ['eQTL', 'pQTL']:
                        genes = [items[phe_idx].split('_')[-1]]
                    else:
                        genes = items[phe_idx].split('_')[-1].split(',')
                    for gene in genes:
                        L.append([items[phe_idx], items[var_idx], gene, items[beta_idx], items[p_idx], qtl_type])
        df = pd.DataFrame(L, columns=['phe_id', 'var_id', 'gene', 'beta', 'pval', 'qtl'])
        df.to_csv(out_file, index=False, sep='\t')

    def get_recurrent_associatoins(self, in_file='QTL_nominal_sig_associations.txt', N=3):
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

    def plot_heatmap_of_recurrent_associatoins(self, in_file='QTL_nominal_sig_associations_recurrent3.txt', cmap='coolwarm', figsize=(4, 8), fontsize=8, customize_cbar=True):
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

        row_indices = g.dendrogram_row.reordered_ind
        for i, idx in enumerate(row_indices):
            label = df['variant'].iloc[idx]
            ax.text(0-0.02, i + 0.5, label, ha='right', va='center', fontsize=fontsize)

        if customize_cbar:
            hm_pos = g.ax_heatmap.get_position()
            cax_left = (hm_pos.x0 + hm_pos.x1)/2 - 0.1
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
