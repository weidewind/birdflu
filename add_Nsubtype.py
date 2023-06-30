import pandas as pd
import re


metafile = '/export/home/popova/workspace/birdflu/data/meta_short.csv'
meta = pd.read_csv(metafile, sep = "\t")
if "NSubtype" not in meta.columns:
	nsubtype = meta['Subtype'].str.extract(r'(N[0-9]+)', expand=True)
	meta.insert(2,'NSubtype', nsubtype)
	meta.to_csv(metafile, sep = "\t", index = False)
