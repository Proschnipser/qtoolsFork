import pandas as pd
from io import StringIO

COLS = [
    'target_name', 'target_accession', 'query_name', 'query_accession',
    'e_value', 'score', 'bias',
    'e_value_best', 'score_best', 'bias_best',
    'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc',
    'description'
]

with open('/data/joscha/output/hmmer_hits_long/GCF_020745735.1_sHemOce1.pat.X.cur._genomic_chunkedTANGO1.tbl') as f:
    content = ''.join(l for l in f if not l.startswith('#'))

df = pd.read_csv(StringIO(content), sep=r'\s+', header=None, names=COLS, usecols=range(19))

# Cast numeric columns
for col in ['e_value', 'score', 'bias', 'e_value_best', 'score_best', 'bias_best', 'exp']:
    df[col] = df[col].astype(float)
for col in ['reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc']:
    df[col] = df[col].astype(int)

print(df.head())
print(df)