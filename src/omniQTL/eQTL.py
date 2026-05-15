from .utils import *
from .qtl import QTL
from .qc import SeqQC

class EQTL(QTL, SeqQC):
    def __init__(self, QTLtools_env='QTLtools'):
        super().__init__()
        self.QTLtools_env = QTLtools_env

    def genome_indexing(self, fa_file='GRCh38.fa', gtf_file='GRCh38.115.gtf', genome_dir='GRCh38_STAR', out_file='indexing_genome.sh', n_threads=4, sjdbOverhang=100):
        cmd = f'STAR --runThreadN {n_threads} --runMode genomeGenerate --genomeDir {genome_dir} --genomeFastaFiles {fa_file} --sjdbGTFfile {gtf_file} --sjdbOverhang {sjdbOverhang}'
        with open(out_file, 'w') as out:
            out.write(cmd + '\n')

    def ranseq_mapping(self, out_file='mapping.sh', fq_dir='.', n_threads=4, flag='NH HI AS nM XS', genome_dir='GRCh38_STAR'):
        self.samples = []
        fqs = sorted([x for x in os.listdir(fq_dir) if x.endswith('.fq.gz') or x.endswith('.fastq.gz')])
        D = {}
        for fq in fqs:
            sample = fq.split('_1')[0].split('_2')[0].split('.fq.gz')[0].split('.fastq.gz')[0]
            D.setdefault(sample, {})
            D[sample].setdefault('fq', [])
            D[sample].setdefault('fq1', [])
            D[sample].setdefault('fq2', [])
            if fq_dir != '.':
                fq = f'{fq_dir}/{fq}'
            elif fq.find('_1') != -1:
                D[sample]['fq1'].append(fq)
            elif fq.find('_2') != -1:
                D[sample]['fq2'].append(fq)
            else:
                D[sample]['fq'].append(fq)

        with open(out_file, 'w') as outfile:
            for sample in D:
                if len(D[sample]['fq1']) > 0 and len(D[sample]['fq2']) > 0 and len(D[sample]['fq1']) == len(D[sample]['fq2']):
                    fs = ','.join(D[sample]['fq1']) + ' ' + ','.join(D[sample]['fq2'])
                    self.paired_end = True
                elif len(D[sample]['fq']) > 0:
                    fs = ','.join(D[sample]['fq'])
                    self.paired_end = False
                else:
                    raise ValueError(f'check fastq file names of {sample}')
                cmd = f'STAR --runThreadN {n_threads} --runMode alignReads --genomeDir {genome_dir} --readFilesIn {fs} --readFilesCommand zcat --outFileNamePrefix {sample}_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outSAMattributes {flag}; samtools index {sample}_Aligned.sortedByCoord.out.bam'
                outfile.write(cmd + '\n')
                self.samples.append(sample)

    def counting_genes(self, out_file='counting_genes.sh', gtf_file='GRCh38.115.gtf', strand=0, min_quality=30, n_threads=4):
        with open(out_file, 'w') as outfile:
            for sample in self.samples:
                bam = f'{sample}_Aligned.sortedByCoord.out.bam'
                if self.paired_end:
                    cmd = f'featureCounts -p -T {n_threads} -s {strand} -Q {min_quality} -a {gtf_file} -o {sample}_geneCounts.txt {bam}'
                else:
                    cmd = f'featureCounts -T {n_threads} -s {strand} -Q {min_quality} -a {gtf_file} -o {sample}_geneCounts.txt {bam}'
                outfile.write(cmd + '\n')

    def counting_exons(self, out_file='counting_exons.sh', gtf_file='GRCh38.115.gtf', strand=0, min_quality=30, n_threads=4):
        with open(out_file, 'w') as outfile:
            for sample in self.samples:
                bam = f'{sample}_Aligned.sortedByCoord.out.bam'
                if self.paired_end:
                    cmd = f'featureCounts -J -f -t exon -O -p -T {n_threads} -s {strand} -Q {min_quality} -a {gtf_file} -o {sample}_exonCounts.txt {bam}'
                else:
                    cmd = f'featureCounts -J -f -t exon -O -T {n_threads} -s {strand} -Q {min_quality} -a {gtf_file} -o {sample}_exonCounts.txt {bam}'

    def get_exonIRER(self, exon_counts_dir='IRER', exon_table='Homo_sapiens.GRCh38.115.Exons', chrom='chr', suffix='.bam', knownJuncsOnly=True):
        exonCountsTables = [f'{exon_counts_dir}/{x}' for x in os.listdir(exon_counts_dir) if x.endswith('_exonCounts.txt')]
        for exonCountsTable in sorted(exonCountsTables):
            jcountsTable = exonCountsTable.replace('_exonCounts.txt', '_jcounts.txt')
            print(f'processing {exonCountsTable} and {jcountsTable}', flush=True)
            IR = {}
            ER = {}
            JC = {}
            EL = []
            knownJuncs = {}
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
                        if fields[0] =='Geneid':
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
                    JCL = np.array([0]*len(sampleExon))
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
                f.write('ExonID\tChr\tStart\tEnd\tGeneID\tGeneName' + '\t' + '\t'.join(['%s_IR\t%s_ER'%(x, x) for x in sampleExon]) + '\n')
                for E in EL:
                    k = '\t'.join([E[1], E[2], E[3]])
                    if k in IR and E[0] in ER:
                        f.write('\t'.join(E) + '\t' + '\t'.join([str(IR[k][n] + '\t' + str(ER[E[0]][n])) for n in range(0, len(IR[k]))]) + '\n')
                    elif k in IR:
                        f.write('\t'.join(E) + '\t' + '\t'.join([str(IR[k][n] + '\t' + '-1') for n in range(0, len(IR[k]))]) + '\n')
                    elif E[0] in ER:
                        f.write('\t'.join(E) + '\t' + '\t'.join(['-1' + '\t' + str(ER[E[0]][n]) for n in range(0, len(ER[E[0]]))]) + '\n')
                    else:
                        f.write('\t'.join(E) + '\t' + '\t'.join(['-1' + '\t' + '-1' for x in sampleExon]) + '\n')

    def merge_counts_tables(self, counts_tables='counts_tables.txt', out_file='eQTL_geneCounts.txt'):
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
                    df3.columns = ['ExonID', 'Chr', 'Start', 'End', 'GeneID', 'GeneName', f'{sample}_IR', f'{sample}_ER']
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

    def get_exonPSI(self, in_file='sQTL_exonIRER.txt'):
        out_file = in_file.split('.txt')[0] + '_PSI.txt'
        with open(in_file) as f, open(out_file, 'w') as fout:
            head = f.readline().strip().split('\t')
            H = head
            for n in range(6, len(head), 2):
                H.append(head[n].split('_IR')[0] + '_PSI')
            fout.write('\t'.join(H) + '\n')

            for line in f:
                line = line.strip()
                fields = line.split('\t')
                L = fields
                for n in range(6, len(fields), 2):
                    IR = float(fields[n])
                    ER = float(fields[n + 1])
                    if IR + ER == 0:
                        PSI = (IR + 1)/(IR + ER + 1)
                    else:
                        PSI = IR/(IR + ER)
                    L.append(f'{PSI:.4f}')
                fout.write('\t'.join(L) + '\n')

    def annotate_gene_name(self, gene_table='GRCh38.115_GenePosType.txt', in_file='eQTL_geneCounts.txt'):
        if not os.path.exists(gene_table):
            raise ValueError(f'{gene_table} is not found, run gtf_to_GenePosType in utils.py on the gtf file first')

        out_file = in_file.replace('.txt', '_geneName.txt')
        df = pd.read_table(gene_table, header=None)
        D = dict(zip(df[0], df[1]))
        with open(in_file, 'r') as f, open(out_file, 'w') as fout:
            head = f.readline().strip().split('\t')
            fout.write('\t'.join(head[0:1] + ['GeneName'] + head[1:]) + '\n')
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                gene_id = fields[0]
                gene_name = D.get(gene_id, gene_id)
                fout.write('\t'.join([gene_id, gene_name] + fields[1:]) + '\n')

    def annotate_exon_name(self, exon_table='GRCh38.115_Exons.txt', in_file='eQTL_exonCounts.txt'):
        if not os.path.exists(exon_table):
            raise ValueError(f'{exon_table} is not found')

        out_file = in_file.replace('.txt', '_geneName.txt')
        D = {}
        G = {}
        with open(exon_table) as f:
            for line in f:
                line = line.strip()
                fields = line.split('\t')
                exonID = '_'.join(fields[1:4])
                geneID = fields[4]
                geneName = fields[5]
                if exonID not in D:
                    D[exonID] = fields[1:4] + fields[0:1] + fields[4:6]
                G[geneID] = geneName

        with open(in_file, 'r') as f, open(out_file, 'w') as fout:
            head = f.readline().strip().split('\t')
            fout.write('ExonID\tGeneName\t'+'\t'.join(head[6:])+'\n')
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

    def counts_to_tpm(self, counts_table='eQTL_geneCounts_geneName.txt', counts_sample='sample_geneCounts.txt', sample_start_idx=2, norm_base=1e6):
        '''
        the same as edger, TPM <- t(t(RPKM) / colSums(RPKM)) * 1e6
        '''
        Length = []
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
        print(f'Total Reads (million):\n{matTotalRaw/norm_base}')
        matLength = (mat.T/Length).T
        matTotal = matLength.sum(axis=0)
        M = matLength/matTotal*norm_base
        df = pd.concat([mat2, M], axis=1)
        df.to_csv(out_file, header=True, index=False, sep='\t', float_format='%.4f')
