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

    

    # Extract full multi-word desc= value
    df['desc'] = df['description'].str.extract(r'desc=(.+)$')

    # Parse key=value pairs into a dict per row, then expand into columns
    def parse_kv(s):
        s = re.sub(r'\s*desc=.+$', '', s)  # strip desc= onwards
        return dict(re.findall(r'(\w+)=(\S+)', s))

    desc_df = df['description'].apply(parse_kv).apply(pd.Series)
    df = pd.concat([df.drop(columns='description'), desc_df], axis=1)

    # Cast float columns
    for col in ['e_value', 'score', 'bias', 'e_value_best', 'score_best', 'bias_best', 'exp']:
        df[col] = df[col].astype(float)
    for col in ['reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'length']:
        df[col] = df[col].astype(int)

    print(df.head())
    break
    splitname=str(filepath).split("_")
    dfdict[splitname[1]][splitname[-1]]=df



# for k1,v1 in dfdict.items():
#     for k2, v2 in v1.items():
#         for i, r in v2.iterrows():
            
