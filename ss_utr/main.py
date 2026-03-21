import argparse
import subprocess
import os
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple

def obtener_argumentos():

    parser = argparse.ArgumentParser()

    parser.add_argument('--gff', type=str, required=True, help="Path to GFF file.")
    parser.add_argument('--gtf', type=str, required=False, help="Path to stringtie output.")
    parser.add_argument('--stringtie', action="store_true", help="Execute stringtie.")
    parser.add_argument('--bams', type=str, required=False, help="Path to bams dir.")
    parser.add_argument('--out', type=str, required=True, help='Path for output csv.')
    parser.add_argument('--all_genes', action="store_true", help="")
    return parser.parse_args()



def obtain_gff(route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''Devuelve todos los cromosomas de la especie junto a su fichero GFF3 como un dataframe'''
        data = pd.read_csv(route, comment='#', sep='\t', header=None, encoding= encoding)
        data.columns = ['chr','db','type','start','end','score','strand','phase','attributes']
        data['old_idx'] = data.index
        return data

def obtain_gene_w_mRNA(dataset: pd.DataFrame, all_genes: bool = False):
    '''
    1. Obtiene los genes que poseen mRNA y elimina el resto de genes. Solo elimina los genes que dan lugar a mRNA, el otro tipo de muestras las almacena.
    Una parte muy importante de esta función es que mantiene la estructura interna del gen que da lugar al mRNA, por tanto, se mantienen clases como exon, CDS, UTR, intrón, etc.
    '''

    list_records: List[Dict] = dataset.to_dict(orient="records")

    new_list_records = []
    for record in list_records:
        record['ID'] = dict( part.split("=", 1) for part in record['attributes'].split(";")).get("ID", 'None')  
        record['Parent'] = dict( part.split("=", 1) for part in record['attributes'].split(";")).get("Parent", 'None') 
        new_list_records.append(record)
    list_records = new_list_records

    # 1.2 filters the gene that don't produce mRNA.
    # 1.2.1 Obtain ids and conexions.
    dict_ids_record: Dict = {}
    gene_mRNA_record: Dict = {}
    for record in list_records:
        dict_ids_record[record['ID']] = record
        if record['type'] == 'mRNA':
            gene_mRNA_record[record['Parent']] = record['ID']


    dict_gen_cds = defaultdict(list)
    # 1.2.2 Obtain records belong to mRNA/gene.
    records_genes_produce_mRNA = []
    remove_for_utr = []
    for record in list_records:
        # Tienen que tener padre, así que se eliminan chr, ri y genes; solo quedan exones, CDS, UTRs, intrones y mRNA o cosas varias. Como se pone la condición de que el padre debe de ser mRNA entonces se busca las subpartes del mRNA.
        if dict_ids_record.get(record['Parent'], -1) != -1 and dict_ids_record[record['Parent']]['type'] == 'mRNA':
            if record['type'] == 'CDS' or record['type'] == 'exon':
                gen_record = dict_ids_record[dict_ids_record[record['Parent']]['Parent']]
                dict_gen_cds[gen_record['ID']].append(record)
            elif (record['type'] == 'three_prime_UTR' or record['type'] == 'five_prime_UTR'):
                gen_record = dict_ids_record[dict_ids_record[record['Parent']]['Parent']]
                if gen_record not in remove_for_utr:
                    remove_for_utr.append(gen_record)
        elif record['type'] == 'mRNA':
            dict_gen_cds[dict_ids_record[record['Parent']]['ID']].append(record)
        # Si es gen y produce mRNA también lo metemos.
        elif record['type'] == "gene" and gene_mRNA_record.get(record['ID'], -1) != -1:
            dict_gen_cds[record['ID']].append(record)
            records_genes_produce_mRNA.append(record)

    if not all_genes:
        for record_gene in remove_for_utr:
            records_genes_produce_mRNA.remove(record_gene)
            del dict_gen_cds[record_gene['ID']]

    dict_mRNA_stuff = {}
    dict_idx_gen = {}
    dict_idx_mRNA = {}
    for key in dict_gen_cds.keys():
        dict_mRNA_stuff[key] = {}
        dict_idx_mRNA[key] = {}
        list_stuff = dict_gen_cds[key]
        for record in list_stuff:
            if record['type'] == 'CDS' or record['type'] == 'exon':
                parent = dict_ids_record[record['Parent']]['ID']
                if dict_mRNA_stuff[key].get(parent, -1) == -1:
                    dict_mRNA_stuff[key][parent] = [record]
                else:
                    dict_mRNA_stuff[key][parent].append(record)
            elif record['type'] == 'gene':
                dict_idx_gen[key] = {}
                dict_idx_gen[key]['old_idx'] = record['old_idx']
                dict_idx_gen[key]['start'] = record['start']
                dict_idx_gen[key]['end'] = record['end']
            elif record['type'] == 'mRNA':
                dict_idx_mRNA[key][record['ID']] = {}
                dict_idx_mRNA[key][record['ID']]['old_idx'] = record['old_idx']
                dict_idx_mRNA[key][record['ID']]['start'] = record['start']
                dict_idx_mRNA[key][record['ID']]['end'] = record['end']


    dict_idx_exon_three = {}
    dict_idx_exon_five = {}
    for key in dict_mRNA_stuff.keys():
        dict_idx_exon_three[key] = {}
        dict_idx_exon_five[key] = {}
        for key2 in dict_mRNA_stuff[key].keys():
            min = np.inf
            min_record = None
            max = 0
            max_record = None
            min_exon = np.inf
            min_record_exon = None
            max_exon = 0
            max_record_exon = None
            for record in dict_mRNA_stuff[key][key2]:
                if record['type'] == 'CDS':
                    if record['start'] <= min:
                        min = record['start']
                        min_record = record
                    if record['end'] >= max:
                        max = record['end']
                        max_record = record
            for record in dict_mRNA_stuff[key][key2][:]:
                if record['type'] == 'exon':
                    if record['start'] <= min_exon and min_record['start'] >= record['start'] and min_record['end'] <= record['end']:
                        min_exon = record['start']
                        min_record_exon = record.copy()
                    if record['end'] >= max_exon and max_record['start'] >= record['start'] and max_record['end'] <= record['end']:
                        max_exon = record['end']
                        max_record_exon = record.copy()
                    dict_mRNA_stuff[key][key2].remove(record)
                    
            min_record['five'] = 'yes'
            max_record['three'] = 'yes'
            dict_idx_exon_three[key][key2] = max_record_exon['old_idx']
            dict_idx_exon_five[key][key2] = min_record_exon['old_idx']

    return records_genes_produce_mRNA, dict_mRNA_stuff, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five

def extract_info_gtf(gtf: pd.DataFrame):

    list_gtf: List[Dict] = gtf.to_dict(orient='records')
    dict_gtf: Dict[str, Dict] = defaultdict(list)
    dict_transcript_exon: Dict[str, Dict[str, List]] = {} # gen, isoforma, exones.

    for record in list_gtf:
        record['attributes'] = record['attributes'].strip()
        if record['type'] == 'transcript':
            id_record: str = [ attribute.split(' ')[1] for attribute in record['attributes'].split(';') if attribute.split(' ')[0] == 'gene_id'][0]
            id_transcript: str = [ attribute.strip().split(' ')[1] for attribute in record['attributes'].split(';') if attribute.strip().split(' ')[0] == 'transcript_id'][0]
            record['ID_gene'] = id_record.replace('"', '')
            record['ID_transcript'] = id_transcript.replace('"', '')
            dict_gtf[record['chr']].append(record)
        else:
            id_record: str = [ attribute.split(' ')[1] for attribute in record['attributes'].split(';') if attribute.split(' ')[0] == 'gene_id'][0]
            id_transcript: str = [ attribute.strip().split(' ')[1] for attribute in record['attributes'].split(';') if attribute.strip().split(' ')[0] == 'transcript_id'][0]
            id_record = id_record.replace('"', '')
            id_transcript = id_transcript.replace('"', '')
            record['ID_gene'] = id_record
            record['ID_transcript'] = id_transcript
            if dict_transcript_exon.get(id_record, -1) == -1:
                dict_transcript_exon[id_record] = {}
                dict_transcript_exon[id_record][id_transcript] = [record]
            elif dict_transcript_exon[id_record].get(id_transcript, -1) == -1:
                dict_transcript_exon[id_record][id_transcript] = [record]
            else:
                dict_transcript_exon[id_record][id_transcript].append(record)

    return dict_gtf, dict_transcript_exon

def compare(list_transcript: List[Dict], list_content_isoform: List[Dict]) -> Tuple:

    list_transcript = sorted(list_transcript, key= lambda x: x['start'])
    list_content_isoform = sorted(list_content_isoform, key= lambda x: x['start'])

    total = 0
    total_utrs = 0
    new_records = []
    new_min = np.inf
    new_max = 0
    min_modify_exon = 0
    max_modify_exon = 0

    for cds in list_content_isoform:

        nucleotide_no_overlap = cds['end'] - cds['start']
        
        if cds.get('three', -1) != -1 and cds.get('five', -1) != -1:
            # TODO: calculo del error, del 3'UTR y del 5'UTR
            for exon in list_transcript:
                nucleotide_no_overlap -= calculate_overlap(cds, exon)
                utr_value, n_records, new_min, min_modify_exon = calculate_five_prime_utr(cds, exon, new_min, min_modify_exon)
                total_utrs += utr_value
                new_records.extend(n_records)
                utr_value, n_records, new_max, max_modify_exon = calculate_three_prime_utr(cds, exon, new_max, max_modify_exon)
                total_utrs += utr_value
                new_records.extend(n_records)
        elif cds.get('three', -1) != -1:
            # TODO: calculo del error y del 3'UTR
            for exon in list_transcript:
                nucleotide_no_overlap -= calculate_overlap(cds, exon)
                utr_value, n_records, new_max, max_modify_exon = calculate_three_prime_utr(cds, exon, new_max, max_modify_exon)
                total_utrs += utr_value
                new_records.extend(n_records)
        elif cds.get('five', -1) != -1:
            # TODO: calculo del error y del 5'UTR
            for exon in list_transcript:
                nucleotide_no_overlap -= calculate_overlap(cds, exon)
                utr_value, n_records, new_min, min_modify_exon = calculate_five_prime_utr(cds, exon, new_min, min_modify_exon)
                total_utrs += utr_value
                new_records.extend(n_records)
        else:
            # TODO: calculo del error.
            for exon in list_transcript:
                nucleotide_no_overlap -= calculate_overlap(cds, exon)
        total += nucleotide_no_overlap

    if new_min == np.inf:
        new_min = min_modify_exon
    if new_max == 0:
        new_max = max_modify_exon

    return total, total_utrs, new_records, new_min, new_max, min_modify_exon, max_modify_exon



def calculate_overlap(cds: Dict, transcript: Dict) -> int:
    c_start, c_end = cds['start'], cds['end']
    t_start, t_end = transcript['start'], transcript['end']

    overlap: int = 0

    if c_end <= t_start or c_start >= t_end:
        return overlap
    else:
        start_limit: int = max(c_start, t_start)
        end_limit: int = min(c_end, t_end)
        overlap = end_limit - start_limit
        return overlap
    
def calculate_five_prime_utr(cds: Dict, transcript: Dict, min_value: int, min_modify_exon: int) -> Tuple:
    nucleotide_utr: int = 0
    new_records: List[Dict] = []
    c_start, c_end = cds['start'], cds['end']
    t_start, t_end = transcript['start'], transcript['end']

    if c_end <= t_start:
        return nucleotide_utr, new_records, min_value, min_modify_exon
    elif t_end <= c_start:
        # TODO: generas nuevo exon y utr. El exon se encuentra a la izquierda del cds y no hay solape
        new_exon = cds.copy()
        new_exon['start'], new_exon['end'], new_exon['type'] = t_start, t_end, 'exon'
        new_utr = cds.copy()
        new_utr['start'], new_utr['end'] = t_start, t_end
        if cds['strand'] == '-':
            new_utr['type'] = 'three_prime_UTR'
        else:
            new_utr['type'] = 'five_prime_UTR'
        new_records.append(new_exon)
        new_records.append(new_utr)
        if min_value >= t_start:
            min_value = t_start 
        nucleotide_utr = t_end - t_start
    else:
        # TODO: hay solape. Generas un nuevo utr y alargas el exon.
        # Nos aseguramos de que el exon sea más largo que el cds por el lado izquierdo, implicando la existencia de un UTR.
        if t_start <= c_start:
            new_utr = cds.copy()
            new_utr['start'], new_utr['end'] = t_start, c_start-1
            if cds['strand'] == '-':
                new_utr['type'] = 'three_prime_UTR'
            else:
                new_utr['type'] = 'five_prime_UTR'
            new_records.append(new_utr)
            min_modify_exon = t_start
            nucleotide_utr = c_start - t_start
    return nucleotide_utr, new_records, min_value, min_modify_exon


def calculate_three_prime_utr(cds: Dict, transcript: Dict, max_value: int, max_modify_exon: int) -> Tuple:
    nucleotide_utr: int = 0
    new_records: List[Dict] = []
    c_start, c_end = cds['start'], cds['end']
    t_start, t_end = transcript['start'], transcript['end']

    if t_end <= c_start:
        return nucleotide_utr, new_records, max_value, max_modify_exon
    elif t_start >= c_end:
        # TODO: generas nuevo exon y utr. El exon se encuentra a la derecha del cds y no hay solape.
        new_exon = cds.copy()
        new_exon['start'], new_exon['end'], new_exon['type'] = t_start, t_end, 'exon'
        new_utr = cds.copy()
        new_utr['start'], new_utr['end'] = t_start, t_end
        if cds['strand'] == '-':
            new_utr['type'] = 'five_prime_UTR'
        else:
            new_utr['type'] = 'three_prime_UTR'
        new_records.append(new_exon)
        new_records.append(new_utr)
        if max_value <= t_end:
            max_value = t_end 
        nucleotide_utr = t_end - t_start
    else:
        # TODO: hay solape. Generas un nuevo utr y alargas el exon.
        # Nos aseguramos de que el exon sea más largo que el cds por el lado derecho, implicando la existencia de un UTR.
        if t_end >= c_end:
            new_utr = cds.copy()
            new_utr['start'], new_utr['end'] = c_end+1, t_end
            if cds['strand'] == '-':
                new_utr['type'] = 'five_prime_UTR'
            else:
                new_utr['type'] = 'three_prime_UTR'
            new_records.append(new_utr)
            max_modify_exon = t_end
            nucleotide_utr = t_end - c_end
    return nucleotide_utr, new_records, max_value, max_modify_exon
            

def change_value(data_frame: pd.DataFrame, list_idx: List[int], list_values: List[int], column: str, value_to_ignore: int):
    for idx, val in zip(list_idx, list_values):
        if val != value_to_ignore:
            data_frame.at[idx, column] = val
    return data_frame

def ejecutar():
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

    df_gff: pd.DataFrame = obtain_gff(gff)
    df_gtf: pd.DataFrame = obtain_gff(gtf)

    df_gff_sorted: pd.DataFrame = df_gff.sort_values(by=['chr', 'start'])
    df_gtf: pd.DataFrame = df_gtf.sort_values(by=['chr', 'start'])

    utrs: List[Dict] = []
    list_idx_three: List[int] = []
    list_idx_five: List[int] = []
    list_idx_mRNA: List[int] = []
    list_idx_gene: List[int] = []
    list_value_idx_three: List[int] = []
    list_value_idx_five: List[int] = []
    list_value_idx_mRNA: List[int] = []
    list_value_idx_gene: List[int] = []


    records_transcript, structure_transcript = extract_info_gtf(df_gtf)
    all_genes = False if args.all_genes == None else True
    records_gene_mRNA, structure_gene, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five = obtain_gene_w_mRNA(df_gff_sorted, all_genes)

    n_gen_without_utrs: int = 0

    for gene in records_gene_mRNA:
        list_transcript = records_transcript[gene['chr']]
        j: int = 0
        gene_iso_best = {}
        no_utr: bool = True
        
        while j < len(list_transcript):
            transcript = list_transcript[j]
            # TODO: quizá es mejor comparar las isoformas del gen contra los transcritos, en vez de el tamaño del gen al completo.
            if gene['end'] > transcript['start'] and transcript['end'] > gene['start'] and gene['strand'] == transcript['strand']:

                transcript_exon: List[Dict] = structure_transcript[transcript['ID_gene']][transcript['ID_transcript']]
                isoform_exon_cds: Dict[str, List] = structure_gene[gene['ID']]

                no_utr = False

                for key_iso in isoform_exon_cds.keys():
                    distance, length_utrs, new_records, min_value, max_value, min_modify_exon, max_modify_exon = compare(transcript_exon, isoform_exon_cds[key_iso])
                    # TODO: if args.alt añades todos los elementos de new_records.
                    # TODO: Se selecciona el exon que se debe de ampliar para el nuevo 3'UTR pero ese exon no es el que se debe de ampliar. Este caso sucede cuando sobre el gen ya están anotados los UTRs y los exones entremos no corresponden a CDS.
                    if gene_iso_best.get(key_iso, -1) == -1:
                        gene_iso_best[key_iso] = {}
                        gene_iso_best[key_iso]['distance'] = distance
                        gene_iso_best[key_iso]['length_utrs'] = length_utrs
                        gene_iso_best[key_iso]['new_records'] = new_records
                        gene_iso_best[key_iso]['min'] = min_value
                        gene_iso_best[key_iso]['max'] = max_value
                        gene_iso_best[key_iso]['min_exon'] = min_modify_exon
                        gene_iso_best[key_iso]['max_exon'] = max_modify_exon
                    elif (gene_iso_best[key_iso]['distance'] > distance) or (gene_iso_best[key_iso]['distance'] == distance and gene_iso_best[key_iso]['length_utrs'] < length_utrs):
                        gene_iso_best[key_iso]['distance'] = distance
                        gene_iso_best[key_iso]['length_utrs'] = length_utrs
                        gene_iso_best[key_iso]['new_records'] = new_records
                        gene_iso_best[key_iso]['min'] = min_value
                        gene_iso_best[key_iso]['max'] = max_value
                        gene_iso_best[key_iso]['min_exon'] = min_modify_exon
                        gene_iso_best[key_iso]['max_exon'] = max_modify_exon

            if list_transcript[j]['start'] > gene['end']:
                break
            j += 1
            
        if no_utr:
            n_gen_without_utrs += 1
        best_min = np.inf
        best_max = 0
        for key in gene_iso_best:
            utrs.extend(gene_iso_best[key]['new_records'])
            list_idx_mRNA.append(dict_idx_mRNA[gene['ID']][key]['old_idx'])
            list_value_idx_mRNA.append((min(gene_iso_best[key]['min'], dict_idx_mRNA[gene['ID']][key]['start']), max(gene_iso_best[key]['max'], dict_idx_mRNA[gene['ID']][key]['end'])))
            list_idx_three.append(dict_idx_exon_three[gene['ID']][key])
            list_value_idx_three.append(gene_iso_best[key]['max_exon']) # TODO: Si coincide el cds con el utr en el mismo exón y el nuevo utr es más corto que el que ya había entonces el nuevo valor del exón no cuadrará.
            list_idx_five.append(dict_idx_exon_five[gene['ID']][key])
            list_value_idx_five.append(gene_iso_best[key]['min_exon']) # TODO: Si coincide el cds con el utr en el mismo exón y el nuevo utr es más corto que el que ya había entonces el nuevo valor del exón no cuadrará.
            if best_min > gene_iso_best[key]['min']:
                best_min = gene_iso_best[key]['min']
            if best_max < gene_iso_best[key]['max']:
                best_max = gene_iso_best[key]['max']
        list_idx_gene.append(dict_idx_gen[gene['ID']]['old_idx'])
        list_value_idx_gene.append((min(best_min, dict_idx_gen[gene['ID']]['start']), max(best_max, dict_idx_gen[gene['ID']]['end'])))

    df_gff = change_value(df_gff, list_idx_gene, [v[0] for v in list_value_idx_gene], 'start', 0)
    df_gff = change_value(df_gff, list_idx_gene, [v[1] for v in list_value_idx_gene], 'end', 0)

    df_gff = change_value(df_gff, list_idx_mRNA, [v[0] for v in list_value_idx_mRNA], 'start', 0)
    df_gff = change_value(df_gff, list_idx_mRNA, [v[1] for v in list_value_idx_mRNA], 'end', 0)

    df_gff = change_value(df_gff, list_idx_three, list_value_idx_three, 'end', 0)
    df_gff = change_value(df_gff, list_idx_five, list_value_idx_five, 'start', 0)

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

    with open(args.out, "w") as f:
        f.write("##gff-version 3\n")
        for _, row in df_gff.iterrows():
            line = "\t".join(str(row[col]) for col in df_gff.columns)
            f.write(line + "\n")

    print("Número de genes: ", len(records_gene_mRNA))
    print("Número de genes sin UTRs añadidos: ", n_gen_without_utrs)
    print("Número de 5'UTR añadidos: ", n_five)
    print("Número de 3'UTR añadidos: ", n_three)

