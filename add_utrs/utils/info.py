
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

        print("\n\n\n")
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

    def merge_element(self):
        return {'n_transcript_p_chr': self.n_transcript_p_chr, 'n_isoformTranscript_p_chr': self.n_isoformTranscript_p_chr, 'n_genes_p_chr': self.n_genes_p_chr, 'n_isoformGenes_p_chr': self.n_isoformGenes_p_chr, 'n_genes_noUTR_p_chr': self.n_genes_noUTR_p_chr, 'n_isoformGenes_noUTR_p_chr': self.n_isoformGenes_noUTR_p_chr}

    def print_result(self, dict_elements_chr = None, merge: bool = False):

        if merge:
            dict_elements_chr = self.merge_element()

        print(f"\t Results:")
        print("#####################################")
        print("#####################################")
        print("#####################################")

        for key_chr in dict_elements_chr['n_transcript_p_chr']:
            print(f"There are {dict_elements_chr['n_transcript_p_chr'][key_chr]} transcripts for the {key_chr} chromosome.")
            print(f"There are {dict_elements_chr['n_isoformTranscript_p_chr'][key_chr]} isoforms of transcripts for the {key_chr} chromosome.")

        for key_chr in dict_elements_chr['n_genes_p_chr']:
            print(f"There are {dict_elements_chr['n_genes_p_chr'][key_chr]} genes for the {key_chr} chromosome.")
            print(f"There are {dict_elements_chr['n_isoformGenes_p_chr'][key_chr]} gene isoforms for the {key_chr} chromosome.")

        for key_chr in dict_elements_chr['n_genes_noUTR_p_chr']:
            print(f"The number of genes to which UTRs have not been added is {dict_elements_chr['n_genes_noUTR_p_chr'][key_chr]} on the {key_chr} chromosome.")
            print(f"The number of gene isoforms to which UTRs have not been added is {dict_elements_chr['n_isoformGenes_noUTR_p_chr'][key_chr]} on the {key_chr} chromosome.")


        print("#####################################")
        print("#####################################")
        print("#####################################\n\n\n")