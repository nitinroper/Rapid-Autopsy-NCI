import os
import csv
import sys
import re
import glob

from collections import defaultdict

v_types =['nonsyn','gain_lost','splice_site']

all_var = defaultdict(lambda: defaultdict(list))
pat_samples = defaultdict(lambda: [])

pat_variants = defaultdict(lambda: [])

## change to the directory where all samples have their strelka called somatic variants
os.chdir('/data/GuhaData/exome/strelka')
## get the directories of all samples based on a pattern matching
sample_dir = glob.glob("RA00*")

for s in sample_dir:
	if '_T' in s:
		pat = s.split('_')[1]
		pat_samples[pat].append(s)
				
for pat in pat_samples:
	pat_variants = defaultdict(lambda: [])
	#pat_samples[pat].sort()
	# sort will not make any difference, because samples will be the key in all_var
	for s in pat_samples[pat]:
		s_variant = []
		sample = '_'.join(s.split('_')[1:])
		#print (pat + '==' + sample)
		s_variants=read_combined_variants(s +'/results/passed.somatic.snvs_indels_snpEff_on_exome_combined.txt')
		all_var[pat][sample]=s_variants
		

## read the combined somatic variant file
def read_combined_variants(f_name):
	v_list =[]
	with open(f_name) as f:
		reader = csv.DictReader(f, delimiter='\t')
		for row in reader:
			if row['TYPE'] in v_types:
				v_list.append('_'.join([row['CHROM'], row['POS'], row['REF'], row['ALT']]))
	return (v_list)


## export the variants count into a jaccard matrix
##
for pat in pat_samples:
	pat_samples[pat].sort()
	## create matrix for each pat
	num_samples = len(pat_samples[pat])
	## create empty matrix
	pat_var = [[ 0.0 for x in range(num_samples)] for y in range(num_samples)]
	## samples
	samples = []
	r =0
	for s_row in pat_samples[pat]:
		sample_row='_'.join(s_row.split('_')[1:])
		samples.append(sample_row)
		c=0
		for s_col in pat_samples[pat]:
			sample_col='_'.join(s_col.split('_')[1:])	
			common_var = set(all_var[pat][sample_row]).intersection(set(all_var[pat][sample_col]))
			if not (len(common_var) == len(all_var[pat][sample_row]) and len(common_var) == len(all_var[pat][sample_col])):
				jaccard = float(len(common_var))/((len(all_var[pat][sample_row])+len(all_var[pat][sample_col]))-len(common_var))
				#print str(jaccard)
				pat_var[r][c] = jaccard
			c=c+1
		r=r+1
	with open(pat+'_jaccard.txt', 'w') as f_out:
		f_out.write('\t'.join(['']+samples) +'\n')
		for r in range(len(pat_var)):
			f_out.write(samples[r]+'\t'+'\t'.join(map(str,pat_var[r]))+'\n')
	f_out.close()
