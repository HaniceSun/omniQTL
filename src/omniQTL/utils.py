import os
import glob
import json
import numpy as np
import pandas as pd
import subprocess
import pyranges
from importlib import resources
import yaml
import datetime
import gzip
import tabix
import liftover
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.sans-serif'] = ['Arial'] + matplotlib.rcParams['font.sans-serif']
matplotlib.rcParams['savefig.dpi'] = 300
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import pylab as plt
import seaborn as sns
import scipy

plt.rcParams.update({
    'figure.figsize': [6.4, 4.8],
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.titlesize': 18,
    'figure.labelsize': 18,
    'legend.fontsize': 10,
})
BASE = resources.files(__package__.split(".")[0])

# modify plt.savefig to save both pdf, svg, and png
def plt_savefig(*args, **kwargs):
	formats = ['.png', '.svg', '.pdf']
	fig = plt.gcf()
	name, ext = os.path.splitext(args[0])
	if ext in formats:
		for x in formats:
			f = name + x
			args = tuple([f]) + args[1:]
			res = fig.savefig(*args, **kwargs)
	else:
		res = fig.savefig(*args, **kwargs)
	fig.canvas.draw_idle()
	return res
plt.savefig = plt_savefig


def get_dbsnp_vcf(vcf_url='https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz', assembly_report_url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt', out_file='dbSNP157_GRCh38.vcf.gz'):
    # vcf_url='https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz'
    # assembly_report_url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt'
    # out_file='dbSNP157_GRCh37.vcf.gz'

    assembly_report_file = assembly_report_url.split('/')[-1]
    mapping_file = assembly_report_file.split('.txt')[0] + '_mapping.txt'
    cmd = f'wget {assembly_report_url} -O {assembly_report_file}'
    subprocess.run(cmd, shell=True)
    df = pd.read_table(assembly_report_file, comment='#', header=None)
    df_sub = df.iloc[:, [6, -1]]
    df_sub.to_csv(mapping_file, sep='\t', index=False, header=False)

    vcf_file = vcf_url.split('/')[-1]
    cmd = f'wget {vcf_url} -O {vcf_file}'
    subprocess.run(cmd, shell=True)
    cmd = f'wget {vcf_url}.tbi -O {vcf_file}.tbi'
    subprocess.run(cmd, shell=True)

    cmd = f'bcftools annotate --rename-chrs {mapping_file} -Oz -o {out_file} {vcf_file}'
    subprocess.run(cmd, shell=True)
    try:
        cmd = f'tabix -p vcf {out_file}'
        subprocess.run(cmd, shell=True, check=True)
    except Exception as e:
        print(f'Error indexing VCF file: {e}, attempting to sort and index again.')
        out_sorted = out_file.replace('.vcf.gz', '.sorted.vcf.gz')
        try:
            cmd = f'bcftools sort {out_file} -Oz -o {out_sorted}; tabix -p vcf {out_sorted}'
            subprocess.run(cmd, shell=True, check=True)
            os.rename(out_sorted, out_file)
            os.rename(out_sorted + '.tbi', out_file + '.tbi')
        except Exception as e:
            print(f'Error sorting and indexing VCF file: {e}. Please check the VCF file for issues.')

def gtf_to_GenePosType(in_file):
    D = {}
    fin = open(in_file)
    fout = open(in_file.split('.gtf')[0] + '_GenePosType.txt', 'w')
    for line in fin:
        if line[0] != '#':
            fields = line.split('\t')
            if fields[2] == 'gene':
                info = fields[8].split(';')
                for item in info:
                    if item.find('gene_id') != -1:
                        gene_id = item.split('"')[1]
                        gene_name = gene_id
                    elif item.find('gene_name') != -1:
                        gene_name = item.split('"')[1]
                    elif item.find('gene_biotype') != -1:
                        gene_biotype = item.split('"')[1]
                D.setdefault(gene_id, [])
                D[gene_id].append('\t'.join([gene_name, fields[0], fields[3], fields[4], fields[6], gene_biotype]))
    fin.close()
    d = sorted(D.items(), key = lambda x : x[0])
    for item in d:
        if len(set(item[1])) != 1:
            print('Warning:gene_id to multiple gene_name')
        else:
            fout.write(item[0] + '\t' + item[1][0] + '\n')
    fout.close()

def annotate_vcf_with_vep(vcf_file, vep_env='vep115', vep_cache='vep_cache', species='homo_sapiens', assembly='GRCh38'):
    out_script = 'run_vep_' + vcf_file.split('.vcf')[0] + '.sh'
    out_file = vcf_file.split('.vcf')[0] + '_vep.vcf'
    with open(out_script, 'w') as fout:
        cmd = f'conda run -n {vep_env} vep --vcf -i {vcf_file} -o {out_file} --species {species} --assembly {assembly} --cache --dir_cache {vep_cache} --canonical --regulatory --show_ref_allele --force_overwrite'
        fout.write(cmd + '\n')
        cmd = f'gzip {out_file}'
        fout.write(cmd + '\n')
