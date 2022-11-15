#!/usr/bin/env python3
import sys
import subprocess
import re
import csv

old_fasta = sys.argv[1]
new_fasta1 = sys.argv[2]
new_fasta2 = sys.argv[3]
old_gfa = sys.argv[4]
new_gfa = sys.argv[5]
old_assembly_info = sys.argv[6]
new_assembly_info = sys.argv[7]
final_fasta = sys.argv[8]

with open(old_fasta,'r') as old_fasta_file,open(new_fasta1,'w') as new_fasta_file1:
	#print(old_sorted_header_order)
	bashCommand = f"seqkit sort --by-length --reverse {old_fasta} > {new_fasta1}"
	process = subprocess.Popen([bashCommand],shell=True)
	output, error = process.communicate()


with open(new_fasta1,'r') as new_fasta_file1, open(new_fasta2,'w') as new_fasta_file2:
	old_sorted_header_order =[]
	for line in new_fasta_file1:
		line=line.strip()
		if '>' in line:
			line=re.sub('>','',line)
			old_sorted_header_order.append(line+'\t')
		else:
			continue
	bashCommand = f"seqkit replace --pattern '.+' --replacement 'contig_{{nr}}' {new_fasta1} > {new_fasta2}"
	process = subprocess.Popen([bashCommand],shell=True)
	output, error = process.communicate()

with open(new_fasta2,'r') as new_fasta_file2:
	new_header_order =[]
	for line in new_fasta_file2:
		line=line.strip()
		if '>' in line:
			line=re.sub('>','',line)
			new_header_order.append(line+'_new')
		else:
			continue
	#print(new_header_order)

print('old_headers')
print(old_sorted_header_order)
print('new_headers')
print(new_header_order)

with open(old_gfa,'r') as old_gfa_file,open(new_gfa,'w') as new_gfa_file:
	old_gfa_string=old_gfa_file.read()
	new_gfa_string=old_gfa_string
	#print(old_gfa_string)
	for i in range(len(old_sorted_header_order)):
		if old_sorted_header_order[i] in new_gfa_string:
			#print('MATCH!')
			new_gfa_string = re.sub(old_sorted_header_order[i],new_header_order[i]+'\t',new_gfa_string)
		else:
			continue
	new_gfa_string=re.sub('_new','',new_gfa_string)
	#print(old_gfa_string[-200:])
	#print(new_gfa_string[-200:])
	new_gfa_file.write(new_gfa_string)


with open(old_assembly_info,'r') as assembly_info_file,open(new_assembly_info,'w') as new_assembly_info_file:
	assembly_info_string=assembly_info_file.read()
	for i in range(len(old_sorted_header_order)):
		if old_sorted_header_order[i] in assembly_info_string:
			#print('MATCH!')
			assembly_info_string = re.sub(old_sorted_header_order[i],new_header_order[i]+'\t',assembly_info_string)
		else:
			continue
	assembly_info_string=re.sub('_new','',assembly_info_string)
	#print(assembly_info_string)
	new_assembly_info_file.write(assembly_info_string)
	

with open(new_assembly_info,'r') as new_assembly_info_file,open(new_fasta2,'r') as new_fasta_file,open(final_fasta,'w') as final_fasta_file:
	final_header_order = [header.replace('_new', '') for header in new_header_order]
	#print(final_header_order)
	new_assembly_info_table = csv.reader(new_assembly_info_file,delimiter='\t')
	for rows in new_assembly_info_table:
		#print(rows)
		if 'contig' not in rows[0]:
			continue
		else:
			name=rows[0]
			length = int(rows[1])
			coverage = int(rows[2])
			circular = rows[3]
		#print(name)
		for i in range(len(final_header_order)):
			if name in final_header_order[i]:
				final_header_order[i]=">"+name+f" coverage={coverage} circular={circular}\n"
			else:
				continue
	#print(final_header_order)
	for line in new_fasta_file:
		if '>' not in line:
			final_fasta_file.write(line)
		else:
			line=line.strip()
			for i in final_header_order:
				#print(line)
				#print(i)
				if line in i:
					final_fasta_file.write(i)



