
from typing import List, Dict, Tuple
from .metrics import Metrics
import numpy as np


class Compare:

    def __init__(self):
        self.instance_metrics = Metrics()


    def compare(self, list_transcript: List[Dict], list_content_isoform: List[Dict]) -> Tuple:

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
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
                    utr_value, n_records, new_min, min_modify_exon = self.instance_metrics.calculate_five_prime_utr(cds, exon, new_min, min_modify_exon)
                    total_utrs += utr_value
                    new_records.extend(n_records)
                    utr_value, n_records, new_max, max_modify_exon = self.instance_metrics.calculate_three_prime_utr(cds, exon, new_max, max_modify_exon)
                    total_utrs += utr_value
                    new_records.extend(n_records)
            elif cds.get('three', -1) != -1:
                # TODO: calculo del error y del 3'UTR
                for exon in list_transcript:
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
                    utr_value, n_records, new_max, max_modify_exon = self.instance_metrics.calculate_three_prime_utr(cds, exon, new_max, max_modify_exon)
                    total_utrs += utr_value
                    new_records.extend(n_records)
            elif cds.get('five', -1) != -1:
                # TODO: calculo del error y del 5'UTR
                for exon in list_transcript:
                    nucleotide_no_overlap -= self.instance_metrics.calculate_overlap(cds, exon)
                    utr_value, n_records, new_min, min_modify_exon = self.instance_metrics.calculate_five_prime_utr(cds, exon, new_min, min_modify_exon)
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
    

    def compare_gff_gtf(self, records_gene_mRNA: List[Dict], records_transcript: Dict[List], structure_transcript: Dict[Dict[List]], structure_gene: Dict[List], dict_idx_gen: Dict[Dict], dict_idx_mRNA: Dict[Dict], dict_idx_exon_three: Dict[Dict], dict_idx_exon_five: Dict[Dict]) -> Tuple:

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
                        distance, length_utrs, new_records, min_value, max_value, min_modify_exon, max_modify_exon = self.compare(transcript_exon, isoform_exon_cds[key_iso])
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


        return utrs, list_idx_gene, list_value_idx_gene, list_idx_mRNA, list_value_idx_mRNA, list_idx_five, list_value_idx_five, list_idx_three, list_value_idx_three, n_gen_without_utrs