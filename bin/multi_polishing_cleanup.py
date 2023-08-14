#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

polishing_stats_out = sys.argv[1]
output=sys.argv[2]

df = pd.read_csv(polishing_stats_out, sep="\t",header=None)
df.columns=['isolate','contig','base_change','polishing_step']
#print(df)
df_wide=df.pivot(columns=['polishing_step'],index=['isolate','contig'],values='base_change').reset_index()
columns= list(df_wide.columns)
#print(df_wide)
print(columns)
if 'medaka' in columns and 'polca' in columns and 'polypolish' in columns:
	print('yep')
	df_wide = df_wide[['isolate', 'contig', 'medaka', 'polypolish', 'polca']]
elif 'medaka' not in columns and 'polca' in columns and 'polypolish' in columns:
	df_wide = df_wide[['isolate', 'contig', 'polypolish', 'polca']]
df_wide.to_csv(output, sep="\t",index=False)