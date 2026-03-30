import pandas as pd
from typing import Dict, List, Tuple
import numpy as np
from collections import defaultdict

class HandleGFF:
    '''
    - The HandleGFF class contains all the functions associated with handling a GFF3.
    '''


    def obtain_gff(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''
        - Read the GFF3 located in 'route' with the encoding specified in 'encoding'. 
          Additionally, add the old_idx column, indicating the location of each sample 
          in the original dataframe (for future modifications).
        '''
        data = pd.read_csv(route, comment='#', sep='\t', header=None, encoding= encoding)
        data.columns = ['chr','db','type','start','end','score','strand','phase','attributes']
        data['old_idx'] = data.index
        return data
    

    # def _extract_info_cds(self, record_cds: Dict, dict_mRNA_cds: Dict, dict_conexions: Dict) -> Dict:
    #     # necesito la estructura de los cds
    #     id_mRNA: str = record_cds['Parent']
    #     id_gene: str = dict_conexions[record_cds['Parent']]['Parent']
    #     if dict_mRNA_cds.get(id_gene, -1) == -1:
    #         dict_mRNA_cds[id_gene] = {}
    #         dict_mRNA_cds[id_gene][id_mRNA] = [record_cds]
    #     elif dict_mRNA_cds[id_gene].get(id_mRNA, -1) == -1:
    #         dict_mRNA_cds[id_gene][id_mRNA] = [record_cds]
    #     else:
    #         dict_mRNA_cds[id_gene][id_mRNA].append(record_cds)
    #     return dict_mRNA_cds

    # def _know_five_three_cds_prime(self, records_cds: Dict, min_value: int, max_value: int, min_record: Dict, max_record: Dict) -> Dict:
    #     # Necesito anotar el mínimo y el máximo (etiquetarlos)
    #     if min_value 

    
    def obtain_gene_w_mRNA(self, dataset: pd.DataFrame, all_genes: bool = False) -> Tuple:
        '''
        1. Obtiene los genes que poseen mRNA y elimina el resto de genes. Solo elimina los genes que dan lugar a mRNA, el otro tipo de muestras las almacena.
        Una parte muy importante de esta función es que mantiene la estructura interna del gen que da lugar al mRNA, por tanto, se mantienen clases como exon, CDS, UTR, intrón, etc.
        '''

        list_records: List[Dict] = dataset.to_dict(orient="records")

        new_list_records = []
        for record in list_records:
            record['ID'] = dict( part.split("=", 1) for part in record['attributes'].split(";") if part != '').get("ID", 'None')  
            record['Parent'] = dict( part.split("=", 1) for part in record['attributes'].split(";") if part != '').get("Parent", 'None') 
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
    
    def change_value(self, data_frame: pd.DataFrame, list_idx: List[int], list_values: List[int], column: str, value_to_ignore: int):
        '''
        - Change the values in the 'column' column of the samples given by the 'list_idx' indices in the 'data_frame' dataframe. 
          The new values are given in the 'list_values' list.
        '''
        for idx, val in zip(list_idx, list_values):
            if val != value_to_ignore:
                data_frame.at[idx, column] = val
        return data_frame
    
    def add_utrs(self, gff: pd.DataFrame, utrs: List[Dict], clean_columns: bool = True) -> Tuple:
        '''
        - Add the utrs from the 'utrs' list to the 'gff' dataframe. 
          If 'clean_columns' is true, the dataframe is formatted to remove extra columns and be in GFF3 format.
        '''
        list_df_gff = gff.to_dict(orient='records')
        n_five: int = 0
        n_three: int = 0
        for idx, utr in enumerate(utrs):
            if utr['type'] == 'five_prime_UTR':
                n_five += 1
            elif utr['type'] == 'three_prime_UTR':
                n_three += 1
            new_idx = utr['old_idx']
            if clean_columns:
                del utr['ID']
                del utr['Parent']
                if utr.get('five', -1) != -1:
                    del utr['five']
                if utr.get('three', -1) != -1:
                    del utr['three']
            list_df_gff.insert(new_idx+idx, utr)

        df_gff = pd.DataFrame(list_df_gff)
        if clean_columns:
            del df_gff['old_idx']

        return df_gff, n_five, n_three
        
    def write_gff(self, gff: pd.DataFrame, route: str):
        '''
        - Write the 'gff' dataframe to the 'route' in GFF3 format.
        '''
        with open(route, "w") as f:
            f.write("##gff-version 3\n")
            for _, row in gff.iterrows():
                line = "\t".join(str(row[col]) for col in gff.columns)
                f.write(line + "\n")
    

class HandleGTF:
    '''
    - The HandleGTF class contains all the functions associated with handling a GTF.
    '''

    def __init__(self):
        self.transcripts = {}

    def obtain_gtf(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''
        - Read the GTF located in 'route' with the encoding specified in 'encoding'. 
          Additionally, add the old_idx column, indicating the location of each sample 
          in the original dataframe (for future modifications).
        '''
        data = pd.read_csv(route, comment='#', sep='\t', header=None, encoding= encoding)
        data.columns = ['chr','db','type','start','end','score','strand','phase','attributes']
        data['old_idx'] = data.index
        return data
    
    def extract_info_gtf(self, gtf: pd.DataFrame) -> Tuple:
        '''
        - Extract information from the 'gtf' dataframe using the attributes column, specifically the gene and transcript identifiers.

          First, group the transcripts for each chromosome into a single dictionary entry. 
          Second, group the exons associated with each transcript.
        '''

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
                if self.transcripts.get(id_record, -1) != -1:
                    self.transcripts[id_record][id_transcript] = record
                else: 
                    self.transcripts[id_record] = {}
                    self.transcripts[id_record][id_transcript] = record
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
    
    def change_value(self, data_frame: pd.DataFrame, list_idx: List[int], list_values: List[int], column: str, value_to_ignore: int):
        '''
        - Change the values in the 'column' column of the samples given by the 'list_idx' indices in the 'data_frame' dataframe. 
          The new values are given in the 'list_values' list.
        '''
        for idx, val in zip(list_idx, list_values):
            if val != value_to_ignore:
                data_frame.at[idx, column] = val
        return data_frame