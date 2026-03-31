from typing import List, Tuple

import pandas as pd

def split_into_chunks(gtf: pd.DataFrame, gff: pd.DataFrame) -> List[Tuple[pd.DataFrame, pd.DataFrame]]:

    list_chunks = [{'df_gff':g1, 'df_gtf':g2}
                for chr_value in set(gff['chr']).intersection(gtf['chr'])
                for g1 in [gff[gff['chr'] == chr_value]]
                for g2 in [gtf[gtf['chr'] == chr_value]]
                ]
    
    return list_chunks