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

    def bar_plot_significant_loci(self, in_file, axes=[0.3, 0.4, 0.6, 0.5], cmap='Dark2', show_numbers=True, figsize=(4, 4)):
        out_file = in_file.split('.txt')[0] + '.pdf'
        cmap = sns.color_palette(cmap)

        df = pd.read_table(in_file, header=0, sep='\t')
        df['StudySampleSize'] = [f"{df['Study'].iloc[n]}\n(N={df['SampleSize'].iloc[n]})" for n in range(df.shape[0])]
    
        fig = plt.figure(figsize=figsize)
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

        y1 = -0.55
        y2 = -0.6
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

    def plot_upset_qtl(self, in_file='QTL_upset_plot_table.txt', cmap='deep'):
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

    def plot_upset_donors(self, in_files=[], out_file='donors_upset.pdf', color='C0'):
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

    def plot_number_merged_peaks(self, in_file='caQTL_number_merged_peaks.txt', cmap='Blues', ylabel='Number of merged peaks (million)', figsize=(4, 4)):
        out_file = in_file.split('.txt')[0] + '_barplot.pdf'
        df = pd.read_table(in_file, header=0, sep='\t')
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        if ylabel.find('million') != -1:
            df['number_of_peaks'] = df['number_of_peaks']/1e6
        if df['peak_method'].nunique() > 1:
            sns.barplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_method', palette=cmap)
            ax.legend(title=None, loc='upper left', prop={'size': 7})
        else:
            sns.barplot(x='peak_type', y='number_of_peaks', data=df, ax=ax, hue='peak_type', palette=cmap, legend=False)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def plot_correlation_number_peaks_and_reads(self, in_file_reads='caQTL_number_mapped_reads.txt', in_file_peaks='caQTL_number_raw_peaks_qvalue.txt', out_file='correlation_number_peaks_and_reads.pdf', cmap='Blues', xlabel='Number of mapped reads (million)', ylabel='Number of raw peaks (million)', figsize=(4, 4), base=1e6):
        df_reads = pd.read_table(in_file_reads, header=0, sep='\t')
        df_peaks = pd.read_table(in_file_peaks, header=0, sep='\t')
        df = pd.merge(df_reads, df_peaks, on='sample')
        df['number_of_peaks'] = df['number_of_peaks']/base
        df['number_of_reads'] = df['number_of_reads']/base
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        sns.regplot(x='number_of_reads', y='number_of_peaks', data=df, ax=ax, color='C0')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_file)

    def get_length_distribution_merged_peaks(self, in_files=[], out_file='caQTL_length_distribution_merged_peaks.txt', params={'consensus':'consensus peaks', 'summitExtended':'summit extended peaks'}):
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

    def plot_length_distribution_merged_peaks(self, in_file='caQTL_length_distribution_merged_peaks.txt', cmap='Set2', xlabel='Length of merged peaks (bp)', figsize=(4, 4), peak_types=['consensus peaks', 'summit extended peaks'], show_summit=True):
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

    def get_variant_annotation(self, in_file='pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_variants.txt', vep_file='pQTL_genotyping_sampleRenamed_rsID_variantFiltered_vep.vcf.gz', canonical_transcript_only=True):
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

    def count_variant_consequence(self, in_files=['pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_annotated.txt'], out_file='QTL_variants_consequence_count.txt', extra_class=['splice', 'UTR', 'inframe', 'frameshift', 'stop', 'start']):
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

    def split_Ensembl_regulatory_annotation(self, in_file='Ensembl_BioMart_RegulatoryAnnotation.txt.gz'):
        # download the regulatory annotation from Ensembl BioMart manually, then split the file into different feature types for downstream analysis
        df = pd.read_table(in_file, header=0, sep='\t', low_memory=False)
        for gi, g in df.groupby('Feature type'):
            out_file = in_file.replace('.txt.gz', f'_{gi}.bed')
            g_sub = g.loc[:, ['Chromosome/scaffold name', 'Start (bp)', 'End (bp)']]
            g_sub.columns = ['ch', 'start', 'end']
            g_sub['ch'] = [f'chr{x}' for x in g_sub['ch']]
            g_sub.to_csv(out_file, header=False, index=False, sep='\t')

    def count_variant_regulatory(self, in_files=['pQTL_nominal-1.0_w1M_PC25_extraInfo_sig_variants.txt'], bed_files=['Ensembl_BioMart_RegulatoryAnnotation_Enhancer.bed'], out_file='QTL_variants_regulatory_count.txt'):
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

    def test_enrichment_using_fisher_exact(self, in_file='QTL_variants_regulatory_count.txt'):
        out_file = in_file.replace('.txt', '_enrichment.txt')
        L = []
        df = pd.read_table(in_file, header=0, sep='\t')
        df['qtl_type'] = [x.split('_')[0] for x in df['file']]
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

    def bar_plot_enrichment(self, in_file='QTL_variants_regulatory_count_enrichment.txt', subset_renaming_file='subset_renaming.txt', xlim=[0, 10], cmap='Dark2', title='Enrichment of QTL significant variants'):
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

    def get_nominal_sig_associations(self, in_files=['caQTL_nominal-1.0_w1k_qvalue_extraInfo_sig.txt.gz',
                                   'eQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz', 'pQTL_nominal-1.0_w1M_PC25_extraInfo_sig.txt.gz'],
                                   out_file='QTL_nomnial_sig_associations.txt', qtl_types={}):
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

    def plot_heatmap_of_recurrent_associatoins(self, in_file='QTL_nominal_sig_associations_recurrent3.txt', cmap='coolwarm', figsize=(4, 8), fontsize=8, customize_cbar=True, genes_highlight=['PTGFRN', 'STARD10', 'PEPD']):
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

    def bar_plot_sig_count_by_params(self, in_file, contrast=['consensus_peaks', 'summit_extended_peaks'], out_suffix='summit_vs_consensus', figsize=(4, 4), cmap='Set2', ylim=[0, 8000], xticklabels=[], fontsize=8):
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

    def prs_overlap_with_qtl(self, in_file, in_file2, flank=1e6, non_sig=False, params={'p_col': 'nom_pval', 'p_threshold': 0.05}):
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
            if non_sig:
                out_file = in_file.split('.csv')[0] + '_' + in_file2.split('.txt')[0] + '_non_sig_overlap.txt'
                df = df[df[params['p_col']] > params['p_threshold']]
            df.to_csv(out_file, index=False, sep='\t')
