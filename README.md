# Common Variants in SARS-CoV-2 are Associated with Severity

[medRxiv Preprint](https://www.medrxiv.org/content/10.1101/2020.12.01.20242149v1)

Running `main.py` will complete the following nine steps:

1. Unpack `2020-10-21_gisaid-fasta.tar.bz2`<sup>[1](#data_footnote)</sup>

```sh
tar -xjvf data/2020-10-21_gisaid-fasta.tar.bz2 -C fasta
```

2. Create VCF files from FASTA files<sup>[1](#data_footnote)</sup>

```sh
python 00_fasta.py
```

3. Create JSON files from `2020-10-21_gisaid-patient.json`<sup>[1](#data_footnote)</sup>

```sh
python 00_recode.py
```

4. Run VCF data processing scripts

```sh
python 01_long.py 02_wide.py
```

5. Join VCF and JSON data

```sh
python 03_join.py
```

6. Calculate variant frequency from VCF data

```sh
python 03_var-freq.py
```

7. Clean data in preparation for modeling

```sh
python 04_clean.py
```

8. Run logistic regression and create Figure 2

```sh
python 05_logit.py
```

9. Create Figures 1, 3, and S1-S3

```sh
python 06_plot_variants.py 06_fig-s1.py 06_fig-s2.py
```
<a name="data_footnote">1</a>: Datasets are not provided. Users must obtain data from GISAID after signing a GISAID Data Use Agreement.
