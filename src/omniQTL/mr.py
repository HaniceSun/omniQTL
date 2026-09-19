"""Two-sample Mendelian randomization helpers.

Prepares exposure (QTL) and outcome (GWAS) summary statistics files for
Mendelian randomization and dispatches the analysis to the bundled
``TwoSampleMR`` R script.
"""

from typing import Any

from .utils import *


class MR:
    """Prepare and run two-sample Mendelian randomization analyses."""

    def __init__(self) -> None:
        """Initialize the MR helper (currently stateless)."""
        pass

    def subset_summary_stats(
        self,
        fetures: list[str] = ['ENSG00000134247_PTGFRN', 'ENSG00000214530_STARD10'],
        summary_stats: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo.txt.gz',
    ) -> None:
        """Extract summary statistics rows for a set of phenotype features.

        Args:
            fetures: Phenotype IDs (first column of `summary_stats`) to
                keep.
            summary_stats: Path to a gzipped, tab-delimited summary
                statistics file. The output is written to the same path
                with ``.txt`` replaced by ``_subset.txt``.
        """
        out_file = summary_stats.split('.txt')[0] + '_subset.txt'
        with gzip.open(summary_stats, 'rt') as f_in, open(out_file, 'w') as f_out:
            header = f_in.readline()
            f_out.write(header)
            for line in f_in:
                if line.split('\t')[0] in fetures:
                    f_out.write(line)

    def get_exposure(
        self,
        feature: str = 'ENSG00000134247_PTGFRN',
        summary_stats: str = 'eQTL_nominal-1.0_w1M_PC25_extraInfo_subset.txt',
        clumping: bool = True,
        clumping_bfile: str = 'eQTL_genotyping_sampleRenamed_rsID_variantFiltered',
        clumping_params: dict[str, Any] = {'r2': 0.001, 'kb': 10000, 'p1': 1},
        phenotype_name: str | None = None,
        cols: dict[str, str] = {
            'var_id': 'SNP',
            'slope': 'beta',
            'slope_se': 'se',
            'nom_pval': 'P',
            'effective_allele': 'effect_allele',
            'non_effective_allele': 'other_allele',
            'var_chr': 'chr',
            'var_from': 'pos',
        },
    ) -> None:
        """Build an MR exposure file for one feature, optionally LD-clumped.

        Subsets `summary_stats` to `feature`, renames columns to the
        ``TwoSampleMR``-expected names via `cols`, and, if `clumping` is
        True, LD-clumps the variants with PLINK using `clumping_bfile` and
        `clumping_params` before writing the clumped exposure file.

        Args:
            feature: Phenotype ID (first column of `summary_stats`) to
                extract.
            summary_stats: Path to a tab-delimited summary statistics file.
            clumping: If True, LD-clump the exposure variants with PLINK.
            clumping_bfile: PLINK binary fileset prefix used as the LD
                reference panel for clumping.
            clumping_params: Dict with keys ``r2``, ``kb``, and ``p1``
                passed to ``plink --clump-r2``, ``--clump-kb``, and
                ``--clump-p1`` respectively.
            phenotype_name: Value written to the ``Phenotype`` column of the
                clumped output; defaults to `feature` when not given.
            cols: Mapping from `summary_stats` column names to the column
                names expected by ``TwoSampleMR`` (e.g. ``SNP``, ``beta``).
        """
        df = pd.read_csv(summary_stats, sep='\t')
        df = df[df.iloc[:, 0] == feature]
        df_subset = df[cols.keys()].dropna()
        df_subset = df_subset.rename(columns=cols)
        out_file = summary_stats.split('.txt')[0] + f'_{feature}_exposure.txt'
        out_file_prefix = out_file.split('.txt')[0]
        out_file_clumped = out_file_prefix + '_clumped.txt'
        df_subset.to_csv(out_file, sep='\t', index=False)
        if clumping:
            cmd = f"plink --bfile {clumping_bfile} --clump {out_file} --clump-p1 {clumping_params['p1']} --clump-r2 {clumping_params['r2']} --clump-kb {clumping_params['kb']} --out {out_file_prefix}"
            print(cmd, flush=True)
            subprocess.run(cmd, shell=True)

            df_clumped = pd.read_csv(out_file_prefix + '.clumped', sep=r'\s+')
            df_subset_clumped = df_subset[df_subset['SNP'].isin(df_clumped['SNP'])]
            if phenotype_name is not None:
                df_subset_clumped['Phenotype'] = phenotype_name
            else:
                df_subset_clumped['Phenotype'] = feature

            df_subset_clumped.to_csv(out_file_clumped, sep='\t', index=False)

    def get_outcome(
        self,
        summary_stats_exposure: str = (
            'eQTL_nominal-1.0_w1M_PC25_extraInfo_subset_ENSG00000134247_PTGFRN_exposure_clumped.txt'
        ),
        summary_stats_outcome: str = 'T2D_GGI_EUR_sumstat_harmoniser.h.tsv.gz',
        phenotype_name: str | None = None,
        cols: dict[str, str] = {
            'rsid': 'SNP',
            'beta': 'beta',
            'standard_error': 'se',
            'p_value': 'P',
            'effect_allele': 'effect_allele',
            'other_allele': 'other_allele',
            'chromosome': 'chr',
            'base_pair_location': 'pos',
        },
        flank: int = 100,
        using_proxy: bool = False,
    ) -> None:
        """Build an MR outcome file matching the exposure SNPs.

        Tabix-queries `summary_stats_outcome` for the genomic range spanned
        by the exposure SNPs (plus `flank`), subsets to variants shared
        with the exposure, renames columns via `cols`, and writes the
        result next to `summary_stats_exposure`.

        Args:
            summary_stats_exposure: Path to the exposure file produced by
                :meth:`get_exposure`, with ``chr``, ``pos``, and ``SNP``
                columns.
            summary_stats_outcome: Path to a tabix-indexed, gzipped,
                tab-delimited GWAS summary statistics file for the outcome.
            phenotype_name: Value written to the ``Phenotype`` column of the
                output; defaults to a name derived from
                `summary_stats_outcome` when not given.
            cols: Mapping from `summary_stats_outcome` column names to the
                column names expected by ``TwoSampleMR``.
            flank: Number of base pairs to extend the query window beyond
                the exposure SNPs' min/max positions.
            using_proxy: If True, look up proxy SNPs for exposure SNPs
                missing from the outcome data. Not yet implemented.
        """
        df_exposure = pd.read_csv(summary_stats_exposure, sep='\t')
        header = pd.read_table(summary_stats_outcome, sep='\t', nrows=0).columns.tolist()
        chrom = df_exposure['chr'].iloc[0].strip('chr')
        pos_min = df_exposure['pos'].min()
        pos_max = df_exposure['pos'].max()
        tb = tabix.open(summary_stats_outcome)
        records = tb.query(chrom, pos_min - flank, pos_max + flank)
        if records is not None:
            if using_proxy:
                # to be implemented: find proxy SNPs for the exposure SNPs and subset the outcome summary stats accordingly
                pass
            else:
                df_outcome = pd.DataFrame(records, columns=header)
                df_outcome_subset = df_outcome[df_outcome['rsid'].isin(df_exposure['SNP'])]
                df_outcome_subset = df_outcome_subset[cols.keys()].dropna()
                df_outcome_subset = df_outcome_subset.rename(columns=cols)
                if phenotype_name is not None:
                    df_outcome_subset['Phenotype'] = phenotype_name
                else:
                    df_outcome_subset['Phenotype'] = summary_stats_outcome.split('_sumstat')[0]

                out_file = summary_stats_exposure.split('.txt')[0] + '_' + summary_stats_outcome.split('.h.')[0] + '_outcome.txt'
                df_outcome_subset.to_csv(out_file, sep='\t', index=False)
        else:
            print(f"No records found for chromosome {chrom} in the range {pos_min-flank}-{pos_max+flank}")

    def run_mr(
        self,
        in_file_exposure: str = (
            'eQTL_nominal-1.0_w1M_PC25_extraInfo_subset_ENSG00000134247_PTGFRN_exposure_clumped.txt'
        ),
        in_file_outcome: str = (
            'eQTL_nominal-1.0_w1M_PC25_extraInfo_subset_ENSG00000134247_PTGFRN_exposure_clumped_'
            'T2D_GGI_EUR_sumstat_harmoniser_outcome.txt'
        ),
        out_file: str = 'ENSG00000134247_PTGFRN_T2D_GGI_EUR_mr_results.txt',
        harmonise_action: int = 1,
        conda_env: str | None = 'QTLtools',
    ) -> None:
        """Run the bundled ``TwoSampleMR`` R script on an exposure/outcome pair.

        Args:
            in_file_exposure: Path to the exposure file (as produced by
                :meth:`get_exposure`).
            in_file_outcome: Path to the outcome file (as produced by
                :meth:`get_outcome`).
            out_file: Path to write the MR results table.
            harmonise_action: Action code passed through to
                ``TwoSampleMR::harmonise_data`` (1, 2, or 3).
            conda_env: Conda environment containing the R/``TwoSampleMR``
                installation, or None to call ``Rscript`` directly.
        """
        mr_script = BASE / 'scripts/mr_TwoSampleMR.R'
        cmd = f'Rscript {mr_script} {in_file_exposure} {in_file_outcome} {out_file} {harmonise_action}'
        if conda_env is not None:
            cmd = f'conda run -n {conda_env} {cmd}'
        print(cmd, flush=True)
        subprocess.run(cmd, shell=True)
