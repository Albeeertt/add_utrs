

from typing import Dict, List


class ProcessTranscript:

    def __init__(self, transcript_overlaps_genes: Dict, transcripts: Dict):
        self.transcript_overlaps_genes = transcript_overlaps_genes
        self.transcripts = transcripts
    
    def valid_genes(self):
        for key_transcript in self.transcript_overlaps_genes.keys():
