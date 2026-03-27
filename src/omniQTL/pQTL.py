from .utils import *
from .qtl import QTL

class PQTL(QTL):
    def __init__(self, QTLtools_env='QTLtools'):
        super().__init__()
        self.QTLtools_env = QTLtools_env

    def impute_missing_values(self):
        pass

    def filter_phenotype_features(self, in_file='pQTL_proteomics_imputed.txt', params={'density':np.log2(0.2), 'sample_percent':0.2}):
        if os.path.exists(in_file):
            df = pd.read_table(in_file, header=0, sep='\t')
        else:
            raise FileNotFoundError(f'{in_file} not found.')

        wh = []
        for n in range(df.shape[0]):
            L = df.iloc[n, 2:] >= params['density']
            flag = False
            if sum(L)/len(L) >= params['sample_percent']:
                flag = True
            wh.append(flag)
        out_file = in_file.replace('.txt', f'_proteinFiltered.txt')
        df = df.loc[wh, ]
        df.to_csv(out_file, header=True, index=False, sep='\t')
