import pandas as pd
from typing import Dict, List, Tuple
import numpy as np
from collections import defaultdict

class HandleGFF:


    def obtain_gff(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''Devuelve todos los cromosomas de la especie junto a su fichero GFF3 como un dataframe'''
        data = pd.read_csv(route, comment='#', sep='\t', header=None, encoding= encoding)
        data.columns = ['chr','db','type','start','end','score','strand','phase','attributes']
        data['old_idx'] = data.index
        return data
    
    def obtain_gene_w_mRNA(self, dataset: pd.DataFrame, all_genes: bool = False) -> Tuple:
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
    
    def change_value(data_frame: pd.DataFrame, list_idx: List[int], list_values: List[int], column: str, value_to_ignore: int):
        for idx, val in zip(list_idx, list_values):
            if val != value_to_ignore:
                data_frame.at[idx, column] = val
        return data_frame
    
    def write_gff(self, gff: pd.DataFrame, route: str):
        with open(route, "w") as f:
            f.write("##gff-version 3\n")
            for _, row in gff.iterrows():
                line = "\t".join(str(row[col]) for col in gff.columns)
                f.write(line + "\n")
    

class HandleGTF:

    def obtain_gtf(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''Devuelve todos los cromosomas de la especie junto a su fichero GFF3 como un dataframe'''
        data = pd.read_csv(route, comment='#', sep='\t', header=None, encoding= encoding)
        data.columns = ['chr','db','type','start','end','score','strand','phase','attributes']
        data['old_idx'] = data.index
        return data
    
    def extract_info_gtf(self, gtf: pd.DataFrame) -> Tuple:

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
    
    def change_value(data_frame: pd.DataFrame, list_idx: List[int], list_values: List[int], column: str, value_to_ignore: int):
        for idx, val in zip(list_idx, list_values):
            if val != value_to_ignore:
                data_frame.at[idx, column] = val
        return data_frame