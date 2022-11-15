#!/usr/bin/env python3

import sys
import subprocess
import re
import csv
import Bio

old_fasta = sys.argv[1]
inter_fasta1=sys.argv[2]
inter_fasta2=sys.argv[3]
final_fasta=sys.argv[4]

with open(old_fasta,'r') as old,open(inter_fasta1,'r') as inter1,open(inter_fasta2,'w') as inter2:
	old_headers=[]
	for row in old:
		if '>' not in row:
			continue
		else:
			#row=row.strip()
			#print(row)
			old_headers.append(row)
	#print(old_headers)
	for line in inter1:
		if '>' not in line:
			inter2.write(line)
		else:
			line=line.strip()
			line=line+" "
			#print(line)
			for val in old_headers:
				if line in val:
					inter2.write(val)
with open(inter_fasta2,'r') as inter2,open(final_fasta,'w') as final:
	from Bio import SeqIO
	for record in SeqIO.parse(inter2, "fasta"):
		header=str(record.description)
		seq=str(record.seq)
		length=str(len(record.seq))
		final.write(">"+header+" length="+length+'\n'+seq+'\n')






