import pandas as pd
import sys
from pathlib import Path
from io import StringIO
from collections import defaultdict

COLS = [
    'target_name', 'target_accession', 'query_name', 'query_accession',
    'e_value', 'score', 'bias',
    'e_value_best', 'score_best', 'bias_best',
    'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc',
    'description'
]

directory="/data/joscha/output/hmmer_hits_long/" #sys.argv[1]
dfdict=defaultdict(dict)
for filepath in Path(directory).rglob("*.tbl"):
    with open(filepath) as f:
        content = ''.join(l for l in f if not l.startswith('#'))

    df = pd.read_csv(StringIO(content), sep=r'\s+', header=None, names=COLS, usecols=range(19))
    print(filepath)
    # Cast numeric columns
    for col in ['e_value', 'score', 'bias', 'e_value_best', 'score_best', 'bias_best', 'exp']:
        df[col] = df[col].astype(float)
    for col in ['reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc']:
        df[col] = df[col].astype(int)

    print(df.head())
    #print(df)
    splitname=str(filepath).split("_")
    dfdict[splitname[1]][splitname[-1]]=df



# for k1,v1 in dfdict.items():
#     for k2, v2 in v1.items():
#         for i, r in v2.iterrows():
#             if r["length"] >