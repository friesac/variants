import os
import pandas as pd
import json
from pathlib import Path
from datetime import date

import tarfile

# Goal of this code is to annotated the variants we identified from previously selected cohorts
# Steps:
# 1. Import list of EPI_ISL IDers of samples in each cohort
# 2. Grab the associated FASTA sequences from GISAID data
# 3. Call VCFs
# 4. Annotate with SnpEff
# 5. Filter with SnpSift

# so now the next step is to iterate through the list of chunks and create the multifasta file
p = Path.home() / "mskar" / "variants" / "data" / "2020-10-21_gisaid-fasta.tar.bz2"

tf = tarfile.open(p)
tf.extractall()

def write_out_single_fasta(filename, targetdir):
    if not os.path.isdir(targetdir):
        os.mkdir(targetdir)

    outFile = None
    for line in open(filename):
        if line[0] == '>':
            if outFile:
                outFile.close()
            line = line.replace(" ", "_").replace("|", "_").replace("/", "_")
            outFile = open(os.path.join(targetdir, line[1:].strip() + '.fa'), 'w')
        if outFile:
            outFile.write(line)
    if outFile:  # if filename is empty
        outFile.close()

def make_uber_vcf(target_dir, out):
    mylist = []
    for entry in os.listdir(target_dir):
        if entry.endswith('.vcf'):
            mylist.append(entry)
    thing = (' '.join(mylist))

    get_ipython().system('bcftools concat {thing} -o {out}')


vcfs = pd.read_csv('/Users/erinmcauley/Desktop/gisaid_processing/pid_list.csv')

# Step 2: Grab FASTA seqs of interest

with open('/Users/erinmcauley/2020-10-21-gisaidpull/2020-10-21-gisaid_pull.json') as f:
    df = pd.DataFrame(json.loads(line) for line in f)

fasta_setup = df[['covv_accession_id', 'sequence']]
fasta_setup['pid_mod'] = fasta_setup['covv_accession_id'].str.extract(pat='L\_(.*)$')
fasta_setup['pid_mod'] = fasta_setup['pid_mod'].astype('int64')

# left join because we want everything from the list of vcfs of interest
seqs_of_interest = pd.merge(vcfs, fasta_setup, left_on='pid', right_on='pid_mod', how='left')
seqs_of_interest_multifasta = seqs_of_interest[['covv_accession_id', 'sequence']]

# Step 2b - get individual FASTA sequences from multifastas
input_dir = '/Users/erinmcauley/Desktop/gisaid_processing/sat_singleton_test/'
target_dir = '/Users/erinmcauley/Desktop/gisaid_processing/split'

for f in os.listdir(input_dir):
    if (f.endswith('.fa')):
        write_out_single_fasta(os.path.join(input_dir, f), target_dir)

# Step 3 - Call individual VCFs
minimap_path = '/Users/erinmcauley/minimap2/minimap2'
k8_path = '~/minimap2/k8'
ref_path = '/Users/erinmcauley/covid19_variation_analysis/NC_045512.2.fa'
paf_path = '/Users/erinmcauley/minimap2/misc/paftools.js'
input_dir = '/Users/erinmcauley/Desktop/gisaid_processing/split/'
target_dir = '/Users/erinmcauley/Desktop/new_dir'
new_format = '.vcf'

Path.mkdir(target_dir)

for f in os.listdir(input_dir):
    if (f.endswith('.fa')):
        print("Now working on: ", f)
        base = os.path.splitext(os.path.basename(f))[0]
        into = os.path.join(input_dir, f)
        new_file = base + new_format
        out_of = os.path.join(target_dir, new_file)
        get_ipython().system(
            '{minimap_path} -c --cs {ref_path} {into} 2>/dev/null | sort -k6,6 -k8,8n | {k8_path} {paf_path} call -L20000 -f {ref_path} - > {out_of}    ')

# check how many files we produced
get_ipython().system('ls {target_dir} | grep -c "EPI"')

# Step 4 - make multi-vcf and Use SnpEff to annotate individual VCFs
target_dir = '/Users/erinmcauley/Desktop/new_dir'
out = '2020-11-30-SnpEffannos3612frompids.vcf'
anno_out = '2020-11-30-3612_SnpEff_annotated.vcf'
make_uber_vcf(target_dir, os.path.join(target_dir, out))
merged_vcf_full_path = os.path.join(target_dir, out)
print(merged_vcf_full_path)
get_ipython().system('java -Xmx8g -jar ~/snpEff/snpEff.jar NC_045512.2 {merged_vcf_full_path} > {anno_out}')

# Step 5 - filter annotations with SnpSift
eff_sift_output = str(date.today()) + 'test.vcf'
output_full_path = os.path.join(target_dir, eff_sift_output)

get_ipython().system(
    'cat {anno_out} | /Users/erinmcauley/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /Users/erinmcauley/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT QUAL DP AF SB "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].EFFECT" "EFF[*].GENE" "EFF[*].CODON" > {output_full_path}')
