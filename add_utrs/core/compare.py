
from typing import List, Dict, Tuple

import numpy as np
from collections import defaultdict
from operator import itemgetter

from .metrics import Metrics



class Compare:
    '''
    - The Compare class contains the functions to find overlaps between the samples in the GFF file and the GTF file.
    '''

    def __init__(self, proportion: float = 4/5, proportion_utrs: float = 1/2, overlap_genes: bool = False):
        self.instance_metrics = Metrics()
        self.transcript_overlap_genes = None
        self.proportion = proportion
        self.proportion_utrs = proportion_utrs
        self.overlap_genes = overlap_genes


    def compare(self, list_transcript: List[Dict], list_content_isoform: List[Dict], start_limit, end_limit) -> Tuple:
        '''
        - Compare the CDS of a gene with the exons of a transcript to extract the UTRs.

           Specifically, this function returns the number of nucleotides in the CDS that do not overlap with the transcript exons (error measure), 
           the size in nucleotides of the generated UTRs, 
           the samples (3'UTR, 5'UTR, and associated exons), 
           the new minimum and maximum size for the isoform, 
           and the new minimum and maximum size for the exon of the first and last CDS.
        '''

        list_transcript = sorted(list_transcript, key= itemgetter('start'))
        list_content_isoform = sorted(list_content_isoform, key= itemgetter('start'))

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
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
                    utr_value, n_records, new_min, min_modify_exon = self.instance_metrics.calculate_five_prime_utr(cds, exon, new_min, min_modify_exon, start_limit)
                    total_utrs += utr_value
                    new_records.extend(n_records)
                    utr_value, n_records, new_max, max_modify_exon = self.instance_metrics.calculate_three_prime_utr(cds, exon, new_max, max_modify_exon, end_limit)
                    total_utrs += utr_value
                    new_records.extend(n_records)
            elif cds.get('three', -1) != -1:
                # TODO: calculo del error y del 3'UTR
                for exon in list_transcript:
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
                    utr_value, n_records, new_max, max_modify_exon = self.instance_metrics.calculate_three_prime_utr(cds, exon, new_max, max_modify_exon, end_limit)
                    total_utrs += utr_value
                    new_records.extend(n_records)
            elif cds.get('five', -1) != -1:
                # TODO: calculo del error y del 5'UTR
                for exon in list_transcript:
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
                    utr_value, n_records, new_min, min_modify_exon = self.instance_metrics.calculate_five_prime_utr(cds, exon, new_min, min_modify_exon, start_limit)
                    total_utrs += utr_value
                    new_records.extend(n_records)
            else:
                # TODO: calculo del error.
                for exon in list_transcript:
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
            total += nucleotide_no_overlap

        if new_min == np.inf:
            new_min = min_modify_exon
        if new_max == 0:
            new_max = max_modify_exon

        return total, total_utrs, new_records, new_min, new_max, min_modify_exon, max_modify_exon
    

    def compare_gff_gtf(self, records_gene_mRNA: List[Dict], records_transcript, structure_transcript, structure_gene, limits_gene, dict_idx_gen, dict_idx_mRNA, dict_idx_exon_three, dict_idx_exon_five) -> Tuple:
        '''
        - It seeks to find the best match between a gene and all possible transcripts. 
          Once found, it stores information on the new UTRs and the updates that need to be made, along with their new values.

          It returns all the newly generated UTRs, 
          along with the indexes of the genes/mRNA/3'UTR/5'UTR that need updating, 
          the new values, 
          and the number of genes for which no UTRs have been added.
        '''
        print("Obtaining UTRs...")
        utrs: List[Dict] = []
        list_idx_three: List[int] = []
        list_idx_five: List[int] = []
        list_idx_mRNA: List[int] = []
        list_idx_gene: List[int] = []
        list_value_idx_three: List[int] = []
        list_value_idx_five: List[int] = []
        list_value_idx_mRNA: List[int] = []
        list_value_idx_gene: List[int] = []

        n_gen_without_utrs: int = 0

        self.transcript_overlap_genes: Dict = defaultdict(list)

        for key_chr in records_transcript.keys():
            for key_strand in records_transcript[key_chr].keys():
                list_transcript_chr_strand = records_transcript[key_chr][key_strand]
                records_transcript[key_chr][key_strand] = sorted(list_transcript_chr_strand, key=itemgetter('start'))

        for gene in records_gene_mRNA:
                
            list_transcript = records_transcript[gene['chr']][gene['strand']]
            if self.overlap_genes:
                start_limit, end_limit = limits_gene[gene['ID']]
            else:
                start_limit, end_limit = -1, np.inf
            length_gene: int = gene['end'] - gene['start']
            condition: float = self.proportion * length_gene
            condition_utrs: float = self.proportion_utrs * length_gene
            j: int = 0
            gene_iso_best = {}
            no_utr: bool = True
            
            while j < len(list_transcript):
                transcript = list_transcript[j]
                length_overlap = max(0, min(gene['end'], transcript['end']) - max(gene['start'], transcript['start']))
                length_five = max(0, gene['start'] - transcript['start'])
                length_three = max(0, transcript['end'] - gene['end'])
                if gene['end'] > transcript['start'] and transcript['end'] > gene['start'] and gene['strand'] == transcript['strand'] and (length_overlap >= condition) and (length_five <= condition_utrs) and (length_three <= condition_utrs):

                    transcript_exon: List[Dict] = structure_transcript[transcript['ID_gene']][transcript['ID_transcript']]
                    isoform_cds: Dict[str, List] = structure_gene[gene['ID']]

                    self.transcript_overlap_genes[transcript['ID_gene']].append(gene)

                    no_utr = False

                    for key_iso in isoform_cds.keys():
                        distance, length_utrs, new_records, min_value, max_value, min_modify_exon, max_modify_exon = self.compare(transcript_exon, isoform_cds[key_iso], start_limit, end_limit)
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
                list_idx_three.append(dict_idx_exon_three[gene['ID']][key]['old_idx'])
                list_value_idx_three.append(max(gene_iso_best[key]['max_exon'], dict_idx_exon_three[gene['ID']][key]['end'])) # TODO: Si coincide el cds con el utr en el mismo exón y el nuevo utr es más corto que el que ya había entonces el nuevo valor del exón no cuadrará.
                list_idx_five.append(dict_idx_exon_five[gene['ID']][key]['old_idx'])
                list_value_idx_five.append(min(gene_iso_best[key]['min_exon'], dict_idx_exon_five[gene['ID']][key]['start'])) # TODO: Si coincide el cds con el utr en el mismo exón y el nuevo utr es más corto que el que ya había entonces el nuevo valor del exón no cuadrará.
                if best_min > gene_iso_best[key]['min']:
                    best_min = gene_iso_best[key]['min']
                if best_max < gene_iso_best[key]['max']:
                    best_max = gene_iso_best[key]['max']
            list_idx_gene.append(dict_idx_gen[gene['ID']]['old_idx'])
            list_value_idx_gene.append((min(best_min, dict_idx_gen[gene['ID']]['start']), max(best_max, dict_idx_gen[gene['ID']]['end'])))


        return utrs, list_idx_gene, list_value_idx_gene, list_idx_mRNA, list_value_idx_mRNA, list_idx_five, list_value_idx_five, list_idx_three, list_value_idx_three, n_gen_without_utrs
    

            
    def get_overlap_transcript_over_all_genes(self):
        '''
        - Groups all overlapping genes against the same transcript. 
        This function should only be called after calling 'compare_gff_gtf'; otherwise, its value will be None.
        '''
        return self.transcript_overlap_genes