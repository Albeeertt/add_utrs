
from typing import Dict, List

import json

class ProcessTranscript:

    def __init__(self, transcript_overlaps_genes: Dict):
        self.transcript_overlaps_genes = transcript_overlaps_genes
    
    def valid_genes(self, write_file: bool = False, route: str = None):
        '''
        Por simplicidad hago uso de de los genes, no de las isoformas de este.
        '''
        new_transcript_overlaps_genes: Dict = {}

        for key_transcript in self.transcript_overlaps_genes.keys():
            list_genes: List[Dict] = self.transcript_overlaps_genes[key_transcript]
            # The condition >1 implies that the number of gene that overlap with the transcript
            # are more than one.
            if len(list_genes) > 1:
                new_transcript_overlaps_genes[key_transcript] = list_genes

        if write_file and route != None:
            self._write_json(new_transcript_overlaps_genes, route)
            
        return new_transcript_overlaps_genes
    
    def _write_json(self, dict_transcript: Dict, route: str) -> None:
        with open(route, 'a') as f:
            json.dump(dict_transcript, f, indent=4)

class ProcessOutputInformation:

    def __init__(self):
        pass