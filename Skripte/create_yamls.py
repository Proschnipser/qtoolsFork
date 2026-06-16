import pandas as pd
from pathlib import Path
import re

yaml_example="/data/joscha/Downloads/MOTH_no._HMMer_template.yaml"
csv= "/data/joscha/Downloads/hmmer_Cyclostomata_Chondrithyes_SP.csv"
df= pd.read_csv(csv)
with open(yaml_example, "r") as f:
    content = f.read()  
for i, r in df.iterrows():
    if 'X' not in r["Sequence cutted"]:
        new_file = f"/data/joscha/Downloads/yamls/boltz_results_MOTH_{r['no.']}_HMMer__{len(r['Sequence cutted'])}_AA_{r['query_name'].split('strict')[0]}.yaml"
        Path(new_file).parent.mkdir(parents=True, exist_ok=True)
        with open(new_file, "w") as f:
            f.write(content+r["Sequence cutted"])
