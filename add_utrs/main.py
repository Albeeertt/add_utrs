
import argparse
import subprocess
import os
import pandas as pd
import time
import resource
from operator import itemgetter


from add_utrs.core.handleFile import HandleGFF, HandleGTF
from add_utrs.core.compare import Compare
from add_utrs.utils.postProcess import ProcessTranscript
from add_utrs.main_parallelize import parallelize_main_part

def obtain_arguments():
    '''
    - Defines the program's arguments.
    '''

    parser = argparse.ArgumentParser()

    parser.add_argument('--gff', type=str, required=True, help="Path to GFF file.")
    parser.add_argument('--gtf', type=str, help="Path to stringtie output.")
    parser.add_argument('--stringtie', action="store_true", help="Execute stringtie.")
    parser.add_argument('--bams', type=str, help="Path to bams dir.")
    parser.add_argument('--out', type=str, required=True, help='Path to output directory')
    parser.add_argument('--all_genes', action="store_true", help="Some genes in your annotation (from the GFF3 file provided as an argument) may already have UTRs annotated. If you include this argument when running the tool, UTRs will be calculated for all genes. If you omit it, only genes that don’t yet have annotated UTRs will be processed.")
    parser.add_argument('--length_overlap', type=float, default=0.8, help="")
    parser.add_argument('--length_utrs', type=float, default=0.5, help="")
    parser.add_argument('--overlap_genes', action='store_true', help="")
    parser.add_argument('--n_cpus', type=int, default=1, help="")
    parser.add_argument('--mem', type=int, default=500, help="")
    return parser.parse_args()

def execute_main_program():
    '''
    - This function contains the program flow, where all calls to the different functions are made. The main steps performed are:
    1. Read the arguments.
    2. Execute stringtie if the argument are present.
    3. Extract information from the GTF. (output of the stringtie)
    4. Extract information from the GFF.
    5. Obtain the UTRs.
    6. Modify the values of the original GFF.
    7. Add the new UTRs.
    8. Write the new GFF.

    The main dependencies are: HandleGFF, HandleGTF, Compare and ProcessTranscript.
    '''


    OUTPUT_GFF3: str = 'output.gff3'
    OUTPUT_OVERLAP: str = 'overlap.json'

    args = obtain_arguments()

    max_heap_size = args.mem * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (max_heap_size, max_heap_size))

    print("Running program...")
    try:
        gff: str = args.gff
        gtf: str = args.gtf
        bams: str = args.bams

        if not os.path.exists(args.out):
            os.mkdir(args.out)

        route_gff3: str = args.out+OUTPUT_GFF3 if args.out.endswith('/') else args.out+'/'+OUTPUT_GFF3
        route_overlap: str = args.out+OUTPUT_OVERLAP if args.out.endswith('/') else args.out+'/'+OUTPUT_OVERLAP

        if args.stringtie:
            string_bams = ""
            for file in os.listdir(bams):
                if file.endswith('.bam'):
                    string_bams += file + " "
            subprocess.run([f'stringtie {string_bams} -o ss_utr.gtf'], shell=True)
            gtf = './ss_utr.gtf'

        instance_handle_gff = HandleGFF()
        instance_handle_gtf = HandleGTF()
        instance_compare = Compare(args.length_overlap, args.length_utrs, args.overlap_genes)

        df_gff: pd.DataFrame = instance_handle_gff.obtain_gff(gff)
        df_gtf: pd.DataFrame = instance_handle_gtf.obtain_gtf(gtf)

        if args.n_cpus > 1:
            records_gene_mRNA, utrs, list_idx_gene, list_value_idx_gene, list_idx_mRNA, list_value_idx_mRNA, list_idx_five, list_value_idx_five, list_idx_three, list_value_idx_three, n_gen_without_utrs = parallelize_main_part(instance_handle_gff, instance_handle_gtf, instance_compare, df_gff, df_gtf, args)
        else:
            records_transcript, structure_transcript = instance_handle_gtf.extract_info_gtf(df_gtf)
            records_gene_mRNA, structure_gene, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five = instance_handle_gff.obtain_gene_w_mRNA(df_gff, args.all_genes)
            dict_limits_genes = instance_handle_gff.extract_all_limits_gene(records_gene_mRNA)
            utrs, list_idx_gene, list_value_idx_gene, list_idx_mRNA, list_value_idx_mRNA, list_idx_five, list_value_idx_five, list_idx_three, list_value_idx_three, n_gen_without_utrs = instance_compare.compare_gff_gtf(records_gene_mRNA, records_transcript, structure_transcript, structure_gene, dict_limits_genes, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five)
        
        df_gff = instance_handle_gff.change_value(df_gff, list_idx_gene, [v[0] for v in list_value_idx_gene], 'start', 0)
        df_gff = instance_handle_gff.change_value(df_gff, list_idx_gene, [v[1] for v in list_value_idx_gene], 'end', 0)

        df_gff = instance_handle_gff.change_value(df_gff, list_idx_mRNA, [v[0] for v in list_value_idx_mRNA], 'start', 0)
        df_gff = instance_handle_gff.change_value(df_gff, list_idx_mRNA, [v[1] for v in list_value_idx_mRNA], 'end', 0)

        df_gff = instance_handle_gff.change_value(df_gff, list_idx_three, list_value_idx_three, 'end', 0)
        df_gff = instance_handle_gff.change_value(df_gff, list_idx_five, list_value_idx_five, 'start', 0)

        if args.n_cpus > 1:
            utrs = sorted(utrs, key = itemgetter('chr', 'start')) # innecesario, el orden sigue estando mal.
        df_gff, n_five, n_three = instance_handle_gff.add_utrs(df_gff, utrs, clean_columns=True)

        instance_handle_gff.write_gff(df_gff, route_gff3)

        print("---")
        print("Number of valid genes: ", len(records_gene_mRNA))
        print("Number of genes without added UTRs: ", n_gen_without_utrs)
        print("Number of 5'UTR added: ", n_five)
        print("Number of 3'UTR added: ", n_three)

        print(".......................................................................................")
    except MemoryError:
        print("Memory limit reached.")

    # instance_ProcessTranscript = ProcessTranscript(instance_compare.get_overlap_transcript_over_all_genes())
    # instance_ProcessTranscript.valid_genes(write_file=True,route=route_overlap)