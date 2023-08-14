#!/usr/bin/env python3

import pandas as pd
import sys

nanostats_out = sys.argv[1]
isolate_id=sys.argv[2]
output=sys.argv[3]

df = pd.read_csv(nanostats_out, sep="\t")
df['isolate']=isolate_id
#print(df)
df_wide=df.pivot(columns='Metrics',index='isolate',values='dataset',)
df_wide.to_csv(output, sep="\t")