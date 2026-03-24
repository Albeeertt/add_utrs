from typing import List, Dict, Tuple

import argparse
import subprocess
import os
import pandas as pd
import time

from ss_utr.core.handleFile import HandleGFF, HandleGTF
from ss_utr.core.compare import Compare

def obtener_argumentos():

    parser = argparse.ArgumentParser()

    parser.add_argument('--gff', type=str, required=True, help="Path to GFF file.")
    parser.add_argument('--gtf', type=str, required=False, help="Path to stringtie output.")
    parser.add_argument('--stringtie', action="store_true", help="Execute stringtie.")
    parser.add_argument('--bams', type=str, required=False, help="Path to bams dir.")
    parser.add_argument('--out', type=str, required=True, help='Path for output csv.')
    parser.add_argument('--all_genes', action="store_true", help="")
    return parser.parse_args()

def ejecutar():
    print("Reading arguments...")
    args = obtener_argumentos()

    gff: str = args.gff
    gtf: str = args.gtf
    bams: str = args.bams

    if args.stringtie:
        string_bams = ""
        for file in os.listdir(bams):
            if file.endswith('.bam'):
                string_bams += file + " "
        subprocess.run([f'stringtie {string_bams} -o ss_utr.gtf'], shell=True)
        gtf = './ss_utr.gtf'

    # time_now = time.time()
    instance_handle_gff = HandleGFF()
    instance_handle_gtf = HandleGTF()
    instance_compare = Compare()

    print("Reading GFF file...")
    df_gff: pd.DataFrame = instance_handle_gff.obtain_gff(gff)
    print("Reading GTF file...")
    df_gtf: pd.DataFrame = instance_handle_gtf.obtain_gtf(gtf)

    df_gff_sorted: pd.DataFrame = df_gff.sort_values(by=['chr', 'start'])
    df_gtf: pd.DataFrame = df_gtf.sort_values(by=['chr', 'start'])

    # print("first step: ", time.time()-time_now)
    # time_now = time.time()
    print("Extracting information from the GTF...")
    records_transcript, structure_transcript = instance_handle_gtf.extract_info_gtf(df_gtf)
    # print("second step: ", time.time()-time_now)
    # time_now = time.time()
    print("Extracting information from the GFF...")
    records_gene_mRNA, structure_gene, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five = instance_handle_gff.obtain_gene_w_mRNA(df_gff_sorted, args.all_genes)
    # print("third step: ",time.time()-time_now)
    # time_now = time.time()
    print("Obtaining UTRs...")
    utrs, list_idx_gene, list_value_idx_gene, list_idx_mRNA, list_value_idx_mRNA, list_idx_five, list_value_idx_five, list_idx_three, list_value_idx_three, n_gen_without_utrs = instance_compare.compare_gff_gtf(records_gene_mRNA, records_transcript, structure_transcript, structure_gene, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five)
    # print("fourth step: ", time.time()-time_now)
    # time_now = time.time()
    print("Changing the sample values...")
    df_gff = instance_handle_gff.change_value(df_gff, list_idx_gene, [v[0] for v in list_value_idx_gene], 'start', 0)
    df_gff = instance_handle_gff.change_value(df_gff, list_idx_gene, [v[1] for v in list_value_idx_gene], 'end', 0)

    df_gff = instance_handle_gff.change_value(df_gff, list_idx_mRNA, [v[0] for v in list_value_idx_mRNA], 'start', 0)
    df_gff = instance_handle_gff.change_value(df_gff, list_idx_mRNA, [v[1] for v in list_value_idx_mRNA], 'end', 0)

    df_gff = instance_handle_gff.change_value(df_gff, list_idx_three, list_value_idx_three, 'end', 0)
    df_gff = instance_handle_gff.change_value(df_gff, list_idx_five, list_value_idx_five, 'start', 0)
    # print("five step: ", time.time()-time_now)
    # time_now = time.time()

    print("Adding the new UTRs...")
    list_df_gff = df_gff.to_dict(orient='records')
    n_five: int = 0
    n_three: int = 0
    for idx, utr in enumerate(utrs):
        if utr['type'] == 'five_prime_UTR':
            n_five += 1
        elif utr['type'] == 'three_prime_UTR':
            n_three += 1
        new_idx = utr['old_idx']
        del utr['ID']
        del utr['Parent']
        if utr.get('five', -1) != -1:
            del utr['five']
        if utr.get('three', -1) != -1:
            del utr['three']
        list_df_gff.insert(new_idx+idx, utr)


    df_gff = pd.DataFrame(list_df_gff)
    del df_gff['old_idx']
    # print("six step: ", time.time()-time_now)
    # time_now = time.time()

    print("Writing the new GFF3...")
    instance_handle_gff.write_gff(df_gff, args.out)
    # print("seven step: ", time.time()- time_now)

    print("---")
    print("Número de genes válidos: ", len(records_gene_mRNA))
    print("Número de genes sin UTRs añadidos: ", n_gen_without_utrs)
    print("Número de 5'UTR añadidos: ", n_five)
    print("Número de 3'UTR añadidos: ", n_three)

