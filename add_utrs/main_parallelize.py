
from typing import List, Tuple

import pandas as pd
from functools import partial

from add_utrs.utils.split import split_into_chunks
from add_utrs.utils.stage import Stage
from add_utrs.utils.scheduler import Scheduler

def parallelize_main_part(instance_handle_gff, instance_handle_gtf, instance_compare, df_gff, df_gtf, args):
    list_chunks: List[Tuple[pd.DataFrame, pd.DataFrame]] = split_into_chunks(df_gtf, df_gff)
    stage1 = Stage(
        func=instance_handle_gtf.extract_info_gtf,
        inputs=['df_gtf'],
        outputs=['records_transcript', 'structure_transcript']
    )
    stage2 = Stage(
        func=partial(instance_handle_gff.obtain_gene_w_mRNA, all_genes=args.all_genes),
        inputs=['df_gff'],
        outputs=['records_gene_mRNA', 'structure_gene', 'dict_idx_gen', 'dict_idx_mRNA', 'dict_idx_exon_three', 'dict_idx_exon_five']
    )
    stage3 = Stage(
        func=instance_handle_gff.extract_all_limits_gene,
        inputs=['records_gene_mRNA'],
        outputs=['dict_limits_genes']
    )
    stage4 = Stage(
        func=instance_compare.compare_gff_gtf,
        inputs=['records_gene_mRNA', 'records_transcript', 'structure_transcript', 'structure_gene', 'dict_limits_genes', 'dict_idx_gen', 'dict_idx_mRNA', 'dict_idx_exon_three', 'dict_idx_exon_five'],
        outputs=['utrs', 'list_idx_gene', 'list_value_idx_gene', 'list_idx_mRNA', 'list_value_idx_mRNA', 'list_idx_five', 'list_value_idx_five', 'list_idx_three', 'list_value_idx_three', 'n_gen_without_utrs']
    )
    result_chr = []
    for i in range(0, len(list_chunks), args.n_cpus):
        chr_to_process = list_chunks[i:i+args.n_cpus]
        scheduler = Scheduler([stage1, stage2, stage3, stage4], chr_to_process)
        result_chunks = scheduler.run()
        result_chr.extend(result_chunks)
    utrs = []
    list_idx_gene = []
    list_value_idx_gene = []
    list_idx_mRNA = []
    list_value_idx_mRNA = []
    list_idx_five = []
    list_value_idx_five = []
    list_idx_three = []
    list_value_idx_three = []
    n_gen_without_utrs = 0
    records_gene_mRNA = []
    for result_cpu in result_chr:
        utrs.extend(result_cpu['utrs'])
        list_idx_gene.extend(result_cpu['list_idx_gene'])
        list_value_idx_gene.extend(result_cpu['list_value_idx_gene'])
        list_idx_mRNA.extend(result_cpu['list_idx_mRNA'])
        list_value_idx_mRNA.extend(result_cpu['list_value_idx_mRNA'])
        list_idx_five.extend(result_cpu['list_idx_five'])
        list_value_idx_five.extend(result_cpu['list_value_idx_five'])
        list_idx_three.extend(result_cpu['list_idx_three'])
        list_value_idx_three.extend(result_cpu['list_value_idx_three'])
        records_gene_mRNA.extend(result_cpu['records_gene_mRNA'])
        n_gen_without_utrs += result_cpu['n_gen_without_utrs']
        

    return records_gene_mRNA, utrs, list_idx_gene, list_value_idx_gene, list_idx_mRNA, list_value_idx_mRNA, list_idx_five, list_value_idx_five, list_idx_three, list_value_idx_three, n_gen_without_utrs