"""Proteomics (pQTL) phenotype preparation."""

from typing import Any

from .utils import *
from .qtl import QTL


class PQTL(QTL):
    """Proteomics-specific phenotype preparation steps for pQTL mapping."""

    def __init__(self, QTLtools_env: str = 'QTLtools') -> None:
        """Initialize the pQTL pipeline.

        Args:
            QTLtools_env: Conda environment containing the ``QTLtools``
                executable.
        """
        super().__init__()
        self.QTLtools_env = QTLtools_env

    def impute_missing_values(self) -> None:
        """Impute missing proteomics values.

        Not yet implemented.
        """
        pass

    def filter_phenotype_features(
        self,
        in_file: str = 'pQTL_proteomics_imputed.txt',
        params: dict[str, Any] = {'density': np.log2(0.2), 'sample_percent': 0.2},
    ) -> None:
        """Filter out proteins with insufficient detection across samples.

        A protein is kept if the fraction of samples with a value at or
        above `params['density']` is at least `params['sample_percent']`.

        Args:
            in_file: Path to a tab-delimited proteomics table; the first
                two columns are feature identifiers and the remaining
                columns are per-sample intensity values.
            params: Dict with keys ``density`` (the minimum per-sample
                intensity threshold, e.g. log2-transformed) and
                ``sample_percent`` (the minimum fraction of samples that
                must meet the threshold).

        Raises:
            FileNotFoundError: If `in_file` does not exist.
        """
        if os.path.exists(in_file):
            df = pd.read_table(in_file, header=0, sep='\t')
        else:
            raise FileNotFoundError(f'{in_file} not found.')

        wh = []
        for n in range(df.shape[0]):
            L = df.iloc[n, 2:] >= params['density']
            flag = False
            if sum(L) / len(L) >= params['sample_percent']:
                flag = True
            wh.append(flag)
        out_file = in_file.replace('.txt', '_proteinFiltered.txt')
        df = df.loc[wh, ]
        df.to_csv(out_file, header=True, index=False, sep='\t')
