
import argparse

def print_arguments(args: argparse.Namespace):

    print("\t\t Applied arguments:")
    print("#####################################")
    print("#####################################")
    print("#####################################")



    for i, (name, value) in enumerate(vars(args).items(), start=1):
        print(f"\t{i}. {name}: {value}")

    print("#####################################")
    print("#####################################")
    print("#####################################")