
# 1: Items to install

For the tool to work, copy and paste the following commands into the terminal depending on whether you are using macOS/Linux or Windows

## 1.1: MacOS and Linux

```bash
git clone https://github.com/Albeeertt/ss_utr.git
cd ss_utr
python3.10 -m venv env_ss_utr
source env_ss_utr/bin/activate
pip install .
```

## 1.2: Windows

```bash
git clone https://github.com/Albeeertt/ss_utr.git
cd ss_utr
python3.10 -m venv env_ss_utr
env_ss_utr\Scripts\activate
pip install .
```

# 2: Arguments

- **--gff**: Path to the GFF file.
- **--gtf**: Path to the GTF file.
- **--out**: Output path where the newly generated GFF3 file will be stored.
- **--all_genes**: Some genes in your annotation (from the GFF3 file provided as an argument) may already have UTRs annotated. If you include this argument when running the tool, UTRs will be calculated for all genes. If you omit it, only genes that don’t yet have annotated UTRs will be processed.
- **--stringtie**: ...
- **bams**: ...

# 3: Example

```bash
ss_utr --gff ../Athaliana_447_Araport11.gene_exons.gff3 --gtf ../Artha_AllRNASeq.STAR.TAIR10.gtf --all_genes --out prueba2.gff3
```
