from typing import Dict, List, Tuple


class Metrics:

    def calculate_overlap(self, cds: Dict, transcript: Dict) -> int:
        '''
        - Returns the overlap (in nucleotides) between two given samples.
        '''
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
        
    def calculate_five_prime_utr(self, cds: Dict, transcript: Dict, min_value: int, min_modify_exon: int, start_limit: int) -> Tuple:
        '''
        - It obtains the number of nucleotides that make up the new UTR, 
        the generated samples (exons and 5'UTR), 
        the minimum value for the isoform and the minimum value for the exon of the last CDS.
        '''
        nucleotide_utr: int = 0
        new_records: List[Dict] = []
        c_start, c_end = cds['start'], cds['end']
        t_start, t_end = transcript['start'], transcript['end']

        if c_end <= t_start:
            return nucleotide_utr, new_records, min_value, min_modify_exon
        elif t_end <= c_start:
            # TODO: generas nuevo exon y utr. El exon se encuentra a la izquierda del cds y no hay solape
            t_start = max(start_limit, t_start)
            if t_end - t_start > 0:
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
                t_start = max(start_limit, t_start)
                if (c_start-1) - t_start > 0:
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
    
    def calculate_three_prime_utr(self, cds: Dict, transcript: Dict, max_value: int, max_modify_exon: int, end_limit: int) -> Tuple:
        '''
        - It obtains the number of nucleotides that make up the new UTR, 
        the generated samples (exons and 3'UTR), 
        the minimum value for the isoform and the minimum value for the exon of the last CDS.
        '''
        
        nucleotide_utr: int = 0
        new_records: List[Dict] = []
        c_start, c_end = cds['start'], cds['end']
        t_start, t_end = transcript['start'], transcript['end']

        if t_end <= c_start:
            return nucleotide_utr, new_records, max_value, max_modify_exon
        elif t_start >= c_end:
            # TODO: generas nuevo exon y utr. El exon se encuentra a la derecha del cds y no hay solape.
            t_end = min(end_limit, t_end)
            if t_end - t_start > 0:
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
                t_end = min(end_limit, t_end)
                if t_end - (c_end+1)> 0:
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