#!/usr/bin/env python3
import sys
import shutil
import re
import csv
import os

output_dir = sys.argv[1]
run_medaka=sys.argv[2]=='true'
run_polypolish=sys.argv[3]=='true'
run_polca=sys.argv[4]=='true'
keep_gfa=sys.argv[5]=='true'
keep_intermediates=sys.argv[6]=='true'

print(run_polca)
print(keep_intermediates)

if run_polca and not keep_intermediates:
	print('ALL GOOD')
	source_folder = f'{output_dir}/polca'
	destination_folder = f'{output_dir}'
	for file_name in os.listdir(source_folder):
	    # construct full file path
	    source = source_folder + "/" + file_name
	    destination = destination_folder + "/" +  file_name
	    #print(source)
	    #print(destination)
	    # move only files
	    if os.path.isfile(source):
	        shutil.move(source, destination)
	        #print('Moved:', file_name)
	shutil.rmtree(f'{output_dir}/flye')
	shutil.rmtree(f'{output_dir}/medaka')
	shutil.rmtree(f'{output_dir}/polypolish')
	shutil.rmtree(f'{output_dir}/polca')

if run_polypolish and not run_polca and not keep_intermediates:
	source_folder = f'{output_dir}/polypolish'
	destination_folder = f'{output_dir}'
	for file_name in os.listdir(source_folder):
	    # construct full file path
	    source = source_folder + "/" + file_name
	    destination = destination_folder + "/" +  file_name
	    #print(source)
	    #print(destination)
	    # move only files
	    if os.path.isfile(source):
	        shutil.move(source, destination)
	        #print('Moved:', file_name)
	shutil.rmtree(f'{output_dir}/flye')
	shutil.rmtree(f'{output_dir}/medaka')
	shutil.rmtree(f'{output_dir}/polypolish')

if run_medaka and not run_polca and not run_polypolish and not keep_intermediates:
	source_folder = f'{output_dir}/medaka'
	destination_folder = f'{output_dir}'
	for file_name in os.listdir(source_folder):
	    # construct full file path
	    source = source_folder + "/" + file_name
	    destination = destination_folder + "/" +  file_name
	    #print(source)
	    #print(destination)
	    # move only files
	    if os.path.isfile(source):
	        shutil.move(source, destination)
	        #print('Moved:', file_name)
	shutil.rmtree(f'{output_dir}/flye')
	shutil.rmtree(f'{output_dir}/medaka')

if not run_medaka and not run_polca and not run_polypolish and not keep_intermediates:
	source_folder = f'{output_dir}/flye'
	destination_folder = f'{output_dir}'
	for file_name in os.listdir(source_folder):
	    # construct full file path
	    source = source_folder + "/" + file_name
	    destination = destination_folder + "/" +  file_name
	    #print(source)
	    #print(destination)
	    # move only files
	    if os.path.isfile(source):
	        shutil.move(source, destination)
	        #print('Moved:', file_name)
	shutil.rmtree(f'{output_dir}/flye')

if not keep_gfa:
	shutil.rmtree(f'{output_dir}/gfa')


