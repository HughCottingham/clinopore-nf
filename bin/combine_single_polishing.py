#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

polishing_stats_out = sys.argv[1]
polishing_step = sys.argv[2]
output=sys.argv[3]

df = pd.read_csv(polishing_stats_out, sep="\t",header=None)
df.columns=['assembly','contig','length','isolate','polishing_step']
print(df)
df['negative_input'] = np.where(df['assembly']=='assembly_input',-df['length'],df['length'])
changes_df=df.groupby(['isolate','contig'])['negative_input'].sum(numeric_only=True).reset_index()
changes_df['polishing_step']=polishing_step
changes_df.to_csv(output,sep='\t',header=None,index=False)
