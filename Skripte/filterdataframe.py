import pandas as pd
import sys

filepath=sys.argv[0]

df= pd.read_csv("/data/joscha/Data/uniprot_2025_01_08_MIAs_02_MOTHs.csv")

df= df.iloc[""]