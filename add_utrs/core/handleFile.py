import pandas as pd
from typing import Dict, List, Tuple
import numpy as np
from collections import defaultdict

class HandleGFF:
    '''
    - The HandleGFF class contains all the functions associated with handling a GFF3.
    '''

    def __init__(self):
        self.struct_genes_in_CHR_and_strand = defaultdict(lambda: defaultdict(list))

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

    def add_id_parent(self, list_records: List[Dict], inner_structure: bool = False) -> Tuple:
        new_list_records: List[Dict] = []
        ids_records: Dict = {}
        gene_produce_mRNA: Dict = {}
        for record in list_records:
            record['ID'] = dict( part.split("=", 1) for part in record['attributes'].split(";") if part != '').get("ID", 'None')
            record['Parent'] = dict( part.split("=", 1) for part in record['attributes'].split(";") if part != '').get("Parent", 'None')
            new_list_records.append(record)
            if inner_structure:
                ids_records[record['ID']] = record
                if record['type'] == 'mRNA':
                    gene_produce_mRNA[record['Parent']] = record['ID']


        return new_list_records, ids_records, gene_produce_mRNA
    
    def obtain_struct_gene(self, list_records: List[Dict], ids_records: Dict[str, Dict], gene_produce_mRNA: Dict[str, str], obtain_genes_produce_mRNA: bool = False, idx_gen: bool = False, idx_mRNA: bool = False) -> Tuple:

        dict_idx_gen: Dict[str, Dict[str, str]]  = {}
        dict_idx_mRNA: Dict[str, Dict[str, Dict[str, str]]] = {}
        records_genes_produce_mRNA: List[Dict] = []
        remove_for_utr: List[Dict] = []
        dict_cds_isoform: Dict[str, Dict[str, List]] = defaultdict(lambda: defaultdict(list))
        dict_exon_isoform: Dict[str, Dict[str, List]] = defaultdict(lambda: defaultdict(list))

        for record in list_records:
            if ids_records.get(record['Parent'], -1) != -1 and ids_records[record['Parent']]['type'] == 'mRNA':
                if record['type'] == 'CDS':
                    mRNA_record = ids_records[record['Parent']]
                    dict_cds_isoform[mRNA_record['Parent']][mRNA_record['ID']].append(record)
                elif record['type'] == 'exon':
                    mRNA_record = ids_records[record['Parent']]
                    dict_exon_isoform[mRNA_record['Parent']][mRNA_record['ID']].append(record)
                elif (record['type'] == 'three_prime_UTR' or record['type'] == 'five_prime_UTR'):
                    gen_record = ids_records[ids_records[record['Parent']]['Parent']]
                    if gen_record not in remove_for_utr:
                        remove_for_utr.append(gen_record)
            elif record['type'] == 'mRNA' and idx_mRNA:
                if dict_idx_mRNA.get(record['Parent'], -1) == -1:
                    dict_idx_mRNA[record['Parent']] = {}
                dict_idx_mRNA[record['Parent']][record['ID']] = {}
                dict_idx_mRNA[record['Parent']][record['ID']]['old_idx'] = record['old_idx']
                dict_idx_mRNA[record['Parent']][record['ID']]['start'] = record['start']
                dict_idx_mRNA[record['Parent']][record['ID']]['end'] = record['end']
            elif record['type'] == "gene" and gene_produce_mRNA.get(record['ID'], -1) != -1:
                self.struct_genes_in_CHR_and_strand[record['chr']][record['strand']].append(record)
                if obtain_genes_produce_mRNA:
                    records_genes_produce_mRNA.append(record)
                if idx_gen:
                    dict_idx_gen[record['ID']] = {}
                    dict_idx_gen[record['ID']]['old_idx'] = record['old_idx']
                    dict_idx_gen[record['ID']]['start'] = record['start']
                    dict_idx_gen[record['ID']]['end'] = record['end']

        return dict_cds_isoform, dict_exon_isoform, records_genes_produce_mRNA, remove_for_utr, dict_idx_gen, dict_idx_mRNA
            
    def know_utrs(self, dict_cds_isoform: Dict, dict_exon_isoform: Dict):

        dict_idx_exon_three = {}
        dict_idx_exon_five = {}
        for key in dict_cds_isoform.keys():
            dict_idx_exon_three[key] = {}
            dict_idx_exon_five[key] = {}
            for key2 in dict_cds_isoform[key].keys():
                min = np.inf
                min_record = None
                max = 0
                max_record = None
                min_exon = np.inf
                min_record_exon = None
                max_exon = 0
                max_record_exon = None
                for record in dict_cds_isoform[key][key2]:
                    if record['start'] <= min:
                        min = record['start']
                        min_record = record
                    if record['end'] >= max:
                        max = record['end']
                        max_record = record
                for record in dict_exon_isoform[key][key2]:
                    if record['start'] <= min_exon and min_record['start'] >= record['start'] and min_record['end'] <= record['end']:
                        min_exon = record['start']
                        min_record_exon = record.copy()
                    if record['end'] >= max_exon and max_record['start'] >= record['start'] and max_record['end'] <= record['end']:
                        max_exon = record['end']
                        max_record_exon = record.copy()
                        
                min_record['five'] = 'yes'
                max_record['three'] = 'yes'
                dict_idx_exon_three[key][key2] = {} 
                dict_idx_exon_three[key][key2]['old_idx'] = max_record_exon['old_idx']
                dict_idx_exon_three[key][key2]['end'] = max_record_exon['end']
                dict_idx_exon_five[key][key2] = {}
                dict_idx_exon_five[key][key2]['old_idx'] = min_record_exon['old_idx']
                dict_idx_exon_five[key][key2]['start'] = min_record_exon['start']

        return dict_idx_exon_three, dict_idx_exon_five


    def obtain_gene_w_mRNA(self, dataset: pd.DataFrame, all_genes: bool = False) -> Tuple:
        '''
        1. Obtiene los genes que poseen mRNA y elimina el resto de genes. Solo elimina los genes que dan lugar a mRNA, el otro tipo de muestras las almacena.
        Una parte muy importante de esta función es que mantiene la estructura interna del gen que da lugar al mRNA, por tanto, se mantienen clases como exon, CDS, UTR, intrón, etc.
        '''

        list_records: List[Dict] = dataset.to_dict(orient="records")

        list_records, dict_ids_record, gene_mRNA_record = self.add_id_parent(list_records, inner_structure=True)

        dict_cds_isoform, dict_exon_isoform, records_genes_produce_mRNA, remove_for_utr, dict_idx_gen, dict_idx_mRNA = self.obtain_struct_gene(list_records, dict_ids_record, gene_mRNA_record, obtain_genes_produce_mRNA=True, idx_gen=True, idx_mRNA=True)

        if not all_genes:
            for record_gene in remove_for_utr:
                records_genes_produce_mRNA.remove(record_gene)
                del dict_cds_isoform[record_gene['ID']]
                del dict_exon_isoform[record_gene['ID']]
                del dict_idx_gen[record_gene['ID']]
                del dict_idx_mRNA[record_gene['ID']]

        dict_idx_exon_three, dict_idx_exon_five = self.know_utrs(dict_cds_isoform, dict_exon_isoform)

        return records_genes_produce_mRNA, dict_cds_isoform, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five
    
    def extract_all_limits_gene(self, list_records: List[Dict]) -> Dict:
        dict_limits_genes = {}
        struct_genes_in_CHR_and_strand_order = {}
        for chr_key in self.struct_genes_in_CHR_and_strand.keys():
            struct_genes_in_CHR_and_strand_order[chr_key] = {}
            for strand_key in self.struct_genes_in_CHR_and_strand[chr_key].keys():
                struct_genes_in_CHR_and_strand_order[chr_key][strand_key] = {}
                list_chr_strand: List[Dict] = self.struct_genes_in_CHR_and_strand[chr_key][strand_key]
                struct_genes_in_CHR_and_strand_order[chr_key][strand_key]['start'] = sorted(list_chr_strand, key=lambda x: x['start'])
                struct_genes_in_CHR_and_strand_order[chr_key][strand_key]['end'] = sorted(list_chr_strand, key=lambda x: x['end'])
        
        for record in list_records:
            dict_list_start_end: Dict[str, List] = struct_genes_in_CHR_and_strand_order[record['chr']][record['strand']]
            limit_start, limit_end = self.obtain_limits_gene(record, dict_list_start_end['end'], dict_list_start_end['start'])
            dict_limits_genes[record['ID']] = (limit_start, limit_end)

        return dict_limits_genes

    
    def obtain_limits_gene(self, record: Dict, list_end_sort: List[Dict], list_start_sort: List[Dict]):
        # TODO: list_end_sort y list_start_sort deben de ser listas que contienen el mismo chr y mismo strand para el record dado.
        start, end = record['start'], record['end']
        limit_start, limit_end = 0, np.inf
        j: int = 0
        while j < len(list_end_sort):
            record_limit_start: Dict = list_end_sort[j]
            if record_limit_start['end'] > start:
                break
            elif record_limit_start['end'] > limit_start:
                limit_start = record_limit_start['end']
            j += 1
        
        j: int = 0
        while j < len(list_start_sort):
            record_limit_end: Dict = list_start_sort[j]
            if record_limit_end['start'] < end:
                j += 1
                continue
            elif record_limit_end['start'] < limit_end:
                limit_end = record_limit_end['start']
            j += 1

        return limit_start+1, limit_end-1



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
        dict_gtf: Dict[str, Dict] = defaultdict( lambda: defaultdict(list))
        dict_transcript_exon: Dict[str, Dict[str, List]] = {} # gen, isoforma, exones.


        for record in list_gtf:
            record['attributes'] = record['attributes'].strip()
            if record['type'] == 'transcript':
                id_record: str = [ attribute.split(' ')[1] for attribute in record['attributes'].split(';') if attribute.split(' ')[0] == 'gene_id'][0]
                id_transcript: str = [ attribute.strip().split(' ')[1] for attribute in record['attributes'].split(';') if attribute.strip().split(' ')[0] == 'transcript_id'][0]
                record['ID_gene'] = id_record.replace('"', '')
                record['ID_transcript'] = id_transcript.replace('"', '')
                dict_gtf[record['chr']][record['strand']].append(record)
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