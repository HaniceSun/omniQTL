"""Prepare and run GWAS/QTL colocalization analysis.

This module prepares GWAS summary statistics for the EBI GWAS-SSF harmoniser, runs the
harmoniser, formats matched QTL/GWAS summary statistics windows for colocalization
(coloc/SuSiE), and merges/filters the resulting per-locus coloc results.
"""
from typing import Any

from .utils import *


class Coloc:
    """Prepare inputs for and orchestrate GWAS/QTL colocalization analysis.

    This class wraps the workflow used to: harmonise GWAS summary statistics with the
    EBI GWAS-SSF harmoniser; format matched QTL and GWAS summary statistics windows into
    the input format expected by the bundled R ``coloc``/SuSiE script; generate the shell
    commands needed to run that script across many loci; and merge/filter the resulting
    per-locus colocalization results.

    Attributes:
        coloc_R_script: Path to the ``coloc_SuSiE.R`` script bundled with the package,
            used to run colocalization analysis for each locus.
    """

    def __init__(self) -> None:
        """Locate the bundled coloc/SuSiE R script.

        Raises:
            FileNotFoundError: If the ``coloc_SuSiE.R`` script cannot be found under
                the package's ``scripts`` directory.
        """
        self.coloc_R_script = BASE / 'scripts/coloc_SuSiE.R'
        if not os.path.exists(self.coloc_R_script):
            raise FileNotFoundError(f"Coloc R script not found at {self.coloc_R_script}.")

    def download_gwas_harmoniser_reference(
        self,
        out_dir: str = 'harmoniser_reference',
        url: str = 'https://ftp.ebi.ac.uk/pub/databases/gwas/harmonisation_resources',
    ) -> None:
        """Download the EBI GWAS-SSF harmoniser reference files.

        Downloads per-chromosome ``.parquet``, ``.vcf.gz`` and ``.vcf.gz.tbi`` reference
        files, plus ``rsID.sql`` and ``md5sums.txt``, from the EBI GWAS harmonisation
        resources site into ``out_dir``. If ``out_dir`` already exists, the download is
        skipped entirely.

        Args:
            out_dir: Directory to download the reference files into. Created if it does
                not already exist.
            url: Base URL of the EBI GWAS harmonisation resources directory.
        """
        if os.path.exists(out_dir):
            print(f"Directory {out_dir} already exists. Skipping download.")
        else:
            cwd = os.getcwd()
            os.makedirs(out_dir, exist_ok=True)
            os.chdir(out_dir)
            chrs = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY', 'chrMT']
            for ch in chrs:
                for suffix in ['.parquet', '.vcf.gz', '.vcf.gz.tbi']:
                    print(f'Downloading {ch}{suffix}...')
                    cmd = f'wget {url}/homo_sapiens-{ch}{suffix}'
                    subprocess.run(cmd, shell=True)
            for x in ['rsID.sql', 'md5sums.txt']:
                print(f'Downloading {x}...')
                cmd = f'wget {url}/{x}'
                subprocess.run(cmd, shell=True)
            os.chdir(cwd)

    def prepare_gwas_harmoniser_input(
        self,
        in_file: str,
        params: dict[str, Any] = {
            'chromosome': 'Chromsome',
            'base_pair_location': 'Position',
            'effect_allele': 'EffectAllele',
            'other_allele': 'NonEffectAllele',
            'p_value': 'Pval',
            'beta': 'Beta',
            'standard_error': 'SE',
        },
        genome_assembly: str = '37',
        coordinate_system: str = '1-based',
        file_type: str = 'GWAS-SSF v0.1',
        is_harmonised: str = 'false',
        is_sorted: str = 'false',
    ) -> None:
        """Reformat a summary statistics file into GWAS-SSF harmoniser input.

        Reads ``in_file``, builds the GWAS-SSF standard columns (as mapped by ``params``,
        filling missing source columns with ``'NA'``) while also keeping any remaining
        original columns, adds a ``chr_pos_GRCh{genome_assembly}`` identifier column, and
        writes the result to a new ``*_harmoniser.tsv`` file alongside a matching
        ``-meta.yaml`` metadata file describing the dataset for the harmoniser.

        Args:
            in_file: Path to the whitespace-delimited input summary statistics file.
            params: Mapping from GWAS-SSF standard column name to the corresponding
                column name in ``in_file``. Columns missing from ``in_file`` are filled
                with the string ``'NA'``.
            genome_assembly: Genome assembly of the input coordinates, used both in the
                derived ``chr_pos_GRCh{genome_assembly}`` column name and written to the
                metadata file.
            coordinate_system: Coordinate system of the input file, written to the
                metadata file.
            file_type: File type identifier written to the metadata file.
            is_harmonised: Whether the input file is already harmonised, written to the
                metadata file.
            is_sorted: Whether the input file is already sorted, written to the metadata
                file.
        """
        df = pd.read_table(in_file, header=0, sep=r'\s+', low_memory=False)
        cols = df.columns
        df_out = pd.DataFrame()
        for key, value in params.items():
            if value in cols:
                df_out[key] = df[value]
            else:
                df_out[key] = 'NA'
        for col in cols:
            if col not in params.values():
                df_out[col] = df[col]
        df_out[f'chr_pos_GRCh{genome_assembly}'] = df_out[['chromosome', 'base_pair_location']].astype(
            str
        ).agg('_'.join, axis=1)
        output_file = in_file.split('.tsv')[0].split('.txt')[0] + '_harmoniser.tsv'
        df_out.to_csv(output_file, sep='\t', index=False)

        output_config = output_file + '-meta.yaml'
        today = datetime.datetime.today().strftime('%Y-%m-%d')
        md5 = subprocess.check_output(f'md5sum {in_file}', shell=True, text=True).strip().split()[0]
        with open(output_config, 'w') as f:
            f.write(f'date_metadata_last_modified: {today}\n')
            f.write(f'genome_assembly: {genome_assembly}\n')
            f.write(f'coordinate_system: {coordinate_system}\n')
            f.write(f'data_file_name: {output_file}\n')
            f.write(f'file_type: {file_type}\n')
            f.write(f'data_file_md5sum: {md5}\n')
            f.write(f'is_harmonised: {is_harmonised}\n')
            f.write(f'is_sorted: {is_sorted}\n')

    def harmonise_gwas_sumstats(
        self,
        in_file: str,
        reference_dir: str = 'harmonisation_reference',
        from_build: str = '37',
        to_build: str = '38',
        cooridnate: str = '1-based',
        chroms: list[str] | None = None,
        version: str = 'v1.1.11',
        profile: str = 'standard,singularity',
        resume: bool = True,
    ) -> None:
        """Run the EBI GWAS-SSF harmoniser Nextflow pipeline on a summary statistics file.

        Builds and runs the ``EBISPOT/gwas-sumstats-harmoniser`` Nextflow command to lift
        over and harmonise ``in_file`` from ``from_build`` to ``to_build``.

        Args:
            in_file: Path to the tab-delimited GWAS-SSF format summary statistics file
                to harmonise.
            reference_dir: Directory containing the harmoniser reference files (see
                `download_gwas_harmoniser_reference`).
            from_build: Genome build of the input coordinates.
            to_build: Genome build to lift over the coordinates to.
            cooridnate: Coordinate system of the input file (e.g. ``'1-based'``).
            chroms: Chromosomes to process. If None, inferred from the unique values of
                the ``chromosome`` column in ``in_file``.
            version: Version of the ``EBISPOT/gwas-sumstats-harmoniser`` pipeline to run.
            profile: Nextflow configuration profile(s) to use.
            resume: Whether to pass the ``-resume`` flag to Nextflow to resume a
                previous run.
        """
        if chroms is None:
            df = pd.read_table(in_file, header=0, sep='\t', low_memory=False)
            chroms = list(sorted(df['chromosome'].astype(str).unique()))
        cmd = (
            f'nextflow run  EBISPOT/gwas-sumstats-harmoniser -profile {profile} -r {version} --harm '
            f'--file {in_file} --ref {reference_dir} --to_build {to_build} --from_build {from_build} '
            f'--coordinate {cooridnate} --chromlist {",".join(chroms)}'
        )
        if resume:
            cmd += ' -resume'
        print(
            'Ruuning the harmoniser requires >28Gb memory. Please ensure you have sufficient resources '
            'before running this command.'
        )
        print(cmd)
        subprocess.run(cmd, shell=True)

    def prepare_coloc_input(
        self,
        sumstats1: str = 'sumstats_from_QTLtools.txt',
        sumstats2: str = 'sumstats_gwas_harmonised.txt',
        sumstats1_type: str = 'qtl',
        sumstats2_type: str = 'gwas',
        sumstats1_sample_size: int = 100,
        sumstats2_sample_size: int = 1000000,
        sumstats1_study_type: str = 'quant',
        sumstats2_study_type: str = 'cc',
        sumstats1_sig_file: str | None = 'sumstats_from_QTLtools_permute_sig.txt',
        pos_flank: float = 1e6,
        out_dir: str | None = None,
        bfile_for_ld: str | None = None,
        external_ld: str | None = None,
        sumstats_suffixes: list[str] = ['_ss1', '_ss2'],
        same_gene_only: bool = False,
        params1: dict[str, Any] = {
            'var_key': ['var_id'],
            'feature_id': 'phe_id',
            'chrom_col': 'var_chr',
            'pos_col': 'var_from',
            'beta_col': 'slope',
            'se_col': 'slope_se',
            'maf_col': 'MAF',
        },
        params2: dict[str, Any] = {
            'var_key': ['rsid'],
            'chrom_col': 'chromosome',
            'pos_col': 'base_pair_location',
            'beta_col': 'beta',
            'se_col': 'standard_error',
            'maf_col': 'MAF',
        },
    ) -> None:
        """Build per-locus coloc/SuSiE input files from two summary statistics files.

        For each phenotype feature in ``sumstats1`` (optionally restricted to
        significant features listed in ``sumstats1_sig_file``), tabix-queries
        ``sumstats2`` in a ``pos_flank``-sized window around that feature (retrying with
        a toggled ``'chr'`` prefix on the chromosome name if the first query fails), then
        matches variants between the two datasets on a constructed variant key. If
        ``sumstats2_type`` is ``'qtl'``, matched variants are further split by
        ``sumstats2``'s own feature/phenotype grouping (optionally restricted to pairs
        that map to the same gene via ``same_gene_only``); if ``'gwas'``, all matched
        variants for the locus are treated as a single group. Each matched
        dataset1/dataset2 pair is reshaped into the ``coloc``-package input format
        (``snp``, ``pos``, ``beta``, ``varbeta``, ``type``, ``N``, ``MAF`` for ``quant``
        study types, and ``phe_id``) and written to a pair of tab-delimited files in
        ``out_dir``, along with a plain list of matched SNP IDs. If ``bfile_for_ld`` is
        given, PLINK is additionally run to compute a square LD (``r``) matrix restricted
        to those SNPs.

        Args:
            sumstats1: Path to the tab-delimited QTL (or other) summary statistics file
                (dataset 1).
            sumstats2: Path to the tabix-indexed, tab-delimited summary statistics file
                to query for matching variants (dataset 2).
            sumstats1_type: Type label for dataset 1 (e.g. ``'qtl'``); informational only.
            sumstats2_type: Type of dataset 2; controls whether matched variants are
                further split by feature (``'qtl'``) or kept as one group (``'gwas'``).
            sumstats1_sample_size: Sample size recorded as the ``N`` column for dataset 1
                in the coloc input.
            sumstats2_sample_size: Sample size recorded as the ``N`` column for dataset 2
                in the coloc input.
            sumstats1_study_type: ``coloc`` study type for dataset 1 (e.g. ``'quant'`` or
                ``'cc'``); when ``'quant'`` a ``MAF`` column is also included.
            sumstats2_study_type: ``coloc`` study type for dataset 2; when ``'quant'`` a
                ``MAF`` column is also included.
            sumstats1_sig_file: Optional whitespace-delimited file (no header) whose
                first column lists significant feature IDs to restrict ``sumstats1`` to.
                Ignored if None or if the file does not exist.
            pos_flank: Number of base pairs to extend on either side of each feature's
                min/max variant position when querying ``sumstats2``.
            out_dir: Directory to write the per-locus coloc input files to. If None,
                derived from ``sumstats1`` and ``sumstats2`` file names (with a
                ``_sameGene`` suffix when ``same_gene_only`` is True). Created if it does
                not already exist.
            bfile_for_ld: Prefix of a PLINK binary fileset (``.bed``/``.bim``/``.fam``)
                used to compute an LD matrix for the matched SNPs of each locus. Mutually
                exclusive with ``external_ld``.
            external_ld: Path to a pre-computed external LD file (not yet implemented).
                Mutually exclusive with ``bfile_for_ld``.
            sumstats_suffixes: Two-element list of suffixes appended to dataset 1 and
                dataset 2 column names respectively, to disambiguate them after merging.
            same_gene_only: If True, restrict matched variant pairs to those where
                dataset 1's feature and dataset 2's feature map to the same gene.
            params1: Column-mapping configuration for dataset 1. Expected keys include
                ``var_key`` (list of columns to concatenate into a variant key),
                ``feature_id``, ``chrom_col``, ``pos_col``, ``beta_col``, ``se_col`` and
                ``maf_col``.
            params2: Column-mapping configuration for dataset 2, analogous to
                ``params1``. When ``sumstats2_type`` is ``'qtl'`` it must also supply a
                ``feature_id`` key.

        Raises:
            ValueError: If both ``bfile_for_ld`` and ``external_ld`` are given.
            FileNotFoundError: If ``bfile_for_ld`` (its ``.bed`` file) or ``external_ld``
                is given but does not exist.
        """
        if out_dir is None:
            out_dir = (
                sumstats1.split('.txt')[0].split('.tsv')[0]
                + '_'
                + sumstats2.split('.txt')[0].split('.tsv')[0]
                + '_coloc'
            )
            if same_gene_only:
                out_dir += '_sameGene'
        os.makedirs(out_dir, exist_ok=True)
        if bfile_for_ld is not None and external_ld is not None:
            raise ValueError("Please provide either bfile_for_ld or external_ld, not both.")
        if bfile_for_ld is not None:
            if not os.path.exists(bfile_for_ld + '.bed'):
                raise FileNotFoundError(f"Genotype data not found at {bfile_for_ld}.")
        elif external_ld is not None:
            if not os.path.exists(external_ld):
                raise FileNotFoundError(f"External LD file not found at {external_ld}.")

        df1 = pd.read_table(sumstats1, sep='\t', header=0)
        if sumstats1_sig_file is not None and os.path.exists(sumstats1_sig_file):
            df_sig = pd.read_table(sumstats1_sig_file, header=None, sep=r'\s+')
            wh = df1.iloc[:, 0].isin(df_sig.iloc[:, 0])
            df1 = df1[wh].copy()

        df1['var_key'] = df1[params1['var_key']].apply(lambda x: '_'.join(x.astype(str)).strip(), axis=1)
        df1.columns = [x + sumstats_suffixes[0] for x in df1.columns]

        tb = tabix.open(sumstats2)
        df2_header = pd.read_table(sumstats2, sep='\t', header=0, nrows=0).columns

        for feature, df1_sub in df1.groupby(params1['feature_id'] + sumstats_suffixes[0]):
            try:
                chrom = str(df1_sub[params1['chrom_col'] + sumstats_suffixes[0]].iloc[0])
                start = int(df1_sub[params1['pos_col'] + sumstats_suffixes[0]].min() - pos_flank)
                end = int(df1_sub[params1['pos_col'] + sumstats_suffixes[0]].max() + pos_flank)
                print(['processing:', feature, chrom, start, end], flush=True)
            except:
                print(f'Error processing feature {feature}. Skipping.')
                continue

            try:
                res = tb.query(chrom, start, end)
                df2_sub = pd.DataFrame(res)
            except:
                if chrom.startswith('chr'):
                    chrom = chrom.split('chr')[-1]
                else:
                    chrom = 'chr' + chrom

                try:
                    res = tb.query(chrom, start, end)
                    df2_sub = pd.DataFrame(res)
                except:
                    print(f'Error querying sumstats2 for feature {feature}. Skipping.')
                    continue
            if df2_sub.shape[0] == 0:
                print(f'No variants found in sumstats2 for feature {feature}. Skipping.')
                continue
            df2_sub.columns = df2_header
            df2_sub['var_key'] = df2_sub[params2['var_key']].apply(
                lambda x: '_'.join(x.astype(str)).strip(), axis=1
            )
            df2_sub.columns = [x + sumstats_suffixes[1] for x in df2_sub.columns]
            if sumstats2_type == 'gwas':
                df2_subs = [['none', df2_sub]]
            elif sumstats2_type == 'qtl':
                df2_subs = []
                for gi, g in df2_sub.groupby(params2['feature_id'] + sumstats_suffixes[1]):
                    df2_subs.append([gi, g])

            for feature2, df2_sub in df2_subs:
                if same_gene_only:
                    gene1 = set([feature.split('_')[1]])
                    gene2 = set(feature2.split('_')[-1].split(','))
                    if not gene1.intersection(gene2):
                        print(f'Skipping feature pair {feature} and {feature2} due to same_gene_only=True.')
                        continue
                df_merged = pd.merge(
                    df1_sub,
                    df2_sub,
                    left_on='var_key' + sumstats_suffixes[0],
                    right_on='var_key' + sumstats_suffixes[1],
                ).sort_values(by=params1['pos_col'] + sumstats_suffixes[0])
                if df_merged.shape[0]:
                    df1x = pd.DataFrame()
                    df2x = pd.DataFrame()
                    df1x['snp'] = df_merged['var_key' + sumstats_suffixes[0]]
                    df1x['pos'] = df_merged[params1['pos_col'] + sumstats_suffixes[0]].astype(int)
                    df1x['beta'] = df_merged[params1['beta_col'] + sumstats_suffixes[0]].astype(float)
                    df1x['varbeta'] = df_merged[params1['se_col'] + sumstats_suffixes[0]].astype(float) ** 2
                    df1x['type'] = sumstats1_study_type
                    df1x['N'] = sumstats1_sample_size
                    if sumstats1_study_type == 'quant':
                        df1x['MAF'] = df_merged[params1['maf_col'] + sumstats_suffixes[0]].astype(float)
                    df1x['phe_id'] = feature

                    df2x['snp'] = df_merged['var_key' + sumstats_suffixes[1]]
                    df2x['pos'] = df_merged[params2['pos_col'] + sumstats_suffixes[1]].astype(int)
                    df2x['beta'] = df_merged[params2['beta_col'] + sumstats_suffixes[1]].astype(float)
                    df2x['varbeta'] = df_merged[params2['se_col'] + sumstats_suffixes[1]].astype(float) ** 2
                    df2x['type'] = sumstats2_study_type
                    df2x['N'] = sumstats2_sample_size
                    if sumstats2_study_type == 'quant':
                        df2x['MAF'] = df_merged[params2['maf_col'] + sumstats_suffixes[1]].astype(float)
                    df2x['phe_id'] = feature2

                    output_file1 = os.path.join(out_dir, f'{feature}-{feature2}{sumstats_suffixes[0]}.txt')
                    output_file2 = os.path.join(out_dir, f'{feature}-{feature2}{sumstats_suffixes[1]}.txt')
                    output_snps = os.path.join(out_dir, f'{feature}-{feature2}_snps.txt')
                    output_ld = os.path.join(out_dir, f'{feature}-{feature2}')

                    df1x.to_csv(output_file1, sep='\t', index=False)
                    df2x.to_csv(output_file2, sep='\t', index=False)
                    df1x['snp'].to_csv(output_snps, sep='\t', index=False, header=False)

                    if bfile_for_ld is not None:
                        cmd = f'plink --bfile {bfile_for_ld} --r square --extract {output_snps} --out {output_ld}'
                        print(cmd)
                        subprocess.run(cmd, shell=True, check=True)
                        try:
                            os.remove(output_ld + '.log')
                            os.remove(output_ld + '.nosex')
                        except:
                            pass
                    elif external_ld is not None:
                        print('Using external LD file, to be implemented.')

    def get_coloc_script(
        self,
        in_dir: str,
        out_script: str = 'run_coloc.sh',
        R_env: str | None = 'QTLtools',
    ) -> None:
        """Generate a shell script of coloc/SuSiE commands for all loci in a directory.

        Scans ``in_dir`` for per-locus LD files (``*.ld``) produced by
        `prepare_coloc_input`, and for each locus that does not already have a
        corresponding ``*_coloc.txt`` output, writes an ``Rscript`` command (optionally
        wrapped in ``conda run -n {R_env}``) invoking `coloc_R_script` on that locus's
        dataset 1, dataset 2 and LD files.

        Args:
            in_dir: Directory containing the per-locus ``*_ss1.txt``, ``*_ss2.txt`` and
                ``*.ld`` files, and where ``*_coloc.txt`` outputs are expected.
            out_script: Path to the shell script file to write the coloc commands to.
            R_env: Name of a conda environment to run the R script in via
                ``conda run -n {R_env}``. If None, the ``Rscript`` command is used
                directly without a conda wrapper.
        """
        in_files = os.listdir(in_dir)
        ld_files = sorted([f for f in os.listdir(in_dir) if f.endswith('.ld')])
        with open(out_script, 'w') as f:
            for ld in ld_files:
                ss1 = ld.split('.ld')[0] + '_ss1.txt'
                ss2 = ld.split('.ld')[0] + '_ss2.txt'
                out_file = ld.split('.ld')[0] + '_coloc.txt'
                if out_file not in in_files:
                    cmd = f'Rscript {self.coloc_R_script} {in_dir}/{ss1} {in_dir}/{ss2} {in_dir}/{ld}'
                    if R_env is not None:
                        cmd = f'conda run -n {R_env} ' + cmd
                    f.write(cmd + '\n')

    def merge_coloc_results(
        self,
        in_dir: str,
        parmas: dict[str, float] = {'PP.H4.abf': 0.8},
    ) -> None:
        """Merge per-locus coloc results and filter by posterior probability.

        Reads and concatenates every ``*_coloc.txt`` file in ``in_dir`` (tagging each
        row with the feature name derived from its file name), writes the concatenated
        table to ``{in_dir}_results.txt``, and writes the subset of rows whose
        ``PP.H4.abf`` column meets or exceeds the threshold in ``parmas`` to
        ``{in_dir}_results_H4.txt``.

        Args:
            in_dir: Directory containing the per-locus ``*_coloc.txt`` result files.
            parmas: Mapping containing the ``'PP.H4.abf'`` posterior probability
                threshold used to filter the merged results.
        """
        out_file = in_dir + '_results.txt'
        out_file_sub = in_dir + '_results_H4.txt'
        coloc_files = sorted([f for f in os.listdir(in_dir) if f.endswith('_coloc.txt')])
        df_list = []
        for coloc in coloc_files:
            df = pd.read_table(os.path.join(in_dir, coloc), sep='\t')
            df['feature'] = coloc.split('_coloc')[0]
            df_list.append(df)
        if not df_list:
            print("No coloc result files found in the directory.")
            return
        df_merged = pd.concat(df_list, axis=0)
        k = 'PP.H4.abf'
        v = parmas['PP.H4.abf']
        df_sub = df_merged[df_merged[k] >= v]
        df_merged.to_csv(out_file, sep='\t', index=False)
        df_sub.to_csv(out_file_sub, sep='\t', index=False)
