
import argparse
from collections import defaultdict


class Info:

    def __init__(self):
        self.n_transcript_p_chr = defaultdict(int)
        self.n_isoformTranscript_p_chr = defaultdict(int)
        self.n_genes_p_chr = defaultdict(int)
        self.n_isoformGenes_p_chr = defaultdict(int)
        self.n_genes_noUTR_p_chr = defaultdict(int)
        self.n_isoformGenes_noUTR_p_chr = defaultdict(int)

        self.genes_struct = None
        

    def print_arguments(self, args: argparse.Namespace):

        print("\t Applied arguments:")
        print("#####################################")
        print("#####################################")
        print("#####################################")

        for i, (name, value) in enumerate(vars(args).items(), start=1):
            print(f"\t{i}. {name}: {value}")

        print("#####################################")
        print("#####################################")
        print("#####################################\n\n\n")

    def info_count_transcripts(self, dict_transcripts):
        for key_chromosome in dict_transcripts.keys():
            for key_transcript in dict_transcripts[key_chromosome].keys():
                self.n_transcript_p_chr[key_chromosome] += 1
                for _ in dict_transcripts[key_chromosome][key_transcript].keys():
                    self.n_isoformTranscript_p_chr[key_chromosome] += 1

    def info_count_genes(self, dict_isoforms):
        self.genes_struct = dict_isoforms
        for key_chromosome in dict_isoforms.keys():
            for key_gene in dict_isoforms[key_chromosome].keys():
                self.n_genes_p_chr[key_chromosome] += 1
                for _ in dict_isoforms[key_chromosome][key_gene].keys():
                    self.n_isoformGenes_p_chr[key_chromosome] += 1

    def info_noUTRadd(self, record):
        self.n_genes_noUTR_p_chr[record['chr']] += 1
        self.n_isoformGenes_noUTR_p_chr[record['chr']] += len(self.genes_struct[record['chr']][record['ID']])
        
    def print_result(self):

        print("\t Results:")
        print("#####################################")
        print("#####################################")
        print("#####################################")

        print("Transcript: ")
        print(self.n_transcript_p_chr)

        print("Isoform transcript: ")
        print(self.n_isoformTranscript_p_chr)

        print("Genes: ")
        print(self.n_genes_p_chr)

        print("Isoform genes: ")
        print(self.n_isoformGenes_p_chr)

        print("----------")

        print("no UTR: " )
        print(self.n_genes_noUTR_p_chr)

        print("no isoform UTR: ")
        print(self.n_isoformGenes_noUTR_p_chr)


        print("#####################################")
        print("#####################################")
        print("#####################################\n\n\n")