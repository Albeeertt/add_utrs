# 1: What is ss_utr?

ss_utr is a tool for annotating the UTRs of each isoform of a gene from the transcripts obtained from **StringTie**.

The first step to understanding this tool is knowing what StringTie is. It is a bioinformatics program that allows the reconstruction and quantification of **transcripts**. From sequencing reads, this tool can reconstruct transcripts, detect isoforms, and quantify gene expression.

And why are these **transcripts** so important for obtaining the UTRs of each gene?

This importance stems from the fact that UTRs are defined by the boundaries of the transcribed RNA. In other words, UTRs are part of the transcript, and if we can obtain the first and last CDS of a gene, along with a given transcript, we can identify the UTRs by observing which parts of the transcript extend beyond these CDS.

Therefore, an annotation file (GFF3) and the StringTie output are necessary. **Both are generated from the same assembly.**

## 1.1: Criteria for selecting the UTRs of a gene.

Given a gene isoform, it is common to find multiple transcripts that overlap with it, either because they correspond to alternative isoforms of the same gene or to transcripts from nearby genes.
Therefore, the best-matching transcript for a given gene isoform is defined according to the following criteria:
1. The transcript must overlap with at least one third of the nucleotides of the isoform.
2. The CDS of the isoform is used as an error metric. Each CDS region must be covered by a corresponding exon in the candidate transcript; otherwise, the number of nucleotides in the non-overlapping CDS regions is counted as error. This allows filtering out transcripts that do not correspond to the isoform under analysis.
3. If multiple transcripts yield the same error, the total length of their UTR regions is used as a tiebreaker. Transcripts with longer UTRs are prioritized.

# 2: Items to install

For the tool to work, copy and paste the following commands into the terminal depending on whether you are using macOS/Linux or Windows

## 2.1: MacOS and Linux

```bash
git clone https://github.com/Albeeertt/ss_utr.git
cd ss_utr
python3.10 -m venv env_ss_utr
source env_ss_utr/bin/activate
pip install .
```

## 2.2: Windows

```bash
git clone https://github.com/Albeeertt/ss_utr.git
cd ss_utr
python3.10 -m venv env_ss_utr
env_ss_utr\Scripts\activate
pip install .
```

# 3: Arguments

- **--gff**: Path to the GFF file.
- **--gtf**: Path to the GTF file.
- **--out**: Output path where the newly generated GFF3 file will be stored.
- **--all_genes**: Some genes in your annotation (from the GFF3 file provided as an argument) may already have UTRs annotated. If you include this argument when running the tool, UTRs will be calculated for all genes. If you omit it, only genes that don’t yet have annotated UTRs will be processed.
- **--stringtie**: ...
- **bams**: ...

# 4: Example

```bash
ss_utr --gff ../Athaliana_447_Araport11.gene_exons.gff3 --gtf ../Artha_AllRNASeq.STAR.TAIR10.gtf --all_genes --out prueba2.gff3
```
