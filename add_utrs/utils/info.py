
import argparse
from collections import defaultdict


class Info:

    def __init__(self):
        self.n_transcript_p_chr = defaultdict(int)

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

    def info_extract_info_gtf(self, record):
        pass
        