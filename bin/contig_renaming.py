#!/usr/bin/env python3

import sys
import subprocess
import re
import csv
import Bio

old_fasta = sys.argv[1]
inter_fasta=sys.argv[2]
new_fasta=sys.argv[3]
final_fasta=sys.argv[4]

with open(old_fasta,'r') as old,open(inter_fasta,'r') as inter,open(new_fasta,'w') as new:
	from Bio import SeqIO
	old_headers=[]
	for row in old:
		if '>' not in row:
			continue
		else:
			#row=row.strip()
			#print(row)
			old_headers.append(row)
	#print(old_headers)
	for line in inter:
		if '>' not in line:
			new.write(line)
		else:
			line=line.strip()
			#print(line)
			for val in old_headers:
				if line in val:
					new.write(val)
with open(new_fasta,'r') as new,open(final_fasta,'w') as final:
	for record in SeqIO.parse(new, "fasta"):
		header=str(record.description)
		seq=str(record.seq)
		length=str(len(record.seq))
		final.write(">"+header+" length="+length+'\n'+seq+'\n')






