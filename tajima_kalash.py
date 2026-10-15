import zarr
import numpy as np 
import numcodecs
import allel
import sys
import csv
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

##This scpits transforms the vcf file to a zarr file and uses the indexing
# vcf_path = 'hgdp_wgs.20190516.full.chr6.vcf.gz'
# zarr_path = 'doing_zarr.zarr'
# allel.vcf_to_zarr(vcf_path, zarr_path, fields='*', log=sys.stdout, overwrite=True)

vcf_path = 'hgdp_wgs.20190516.full.chr6.vcf.gz'
zarr_path = 'doing_zarr.zarr'

#looking at the Zarr file ##################################################
callset = zarr.open_group(zarr_path, mode='r')
callset.tree(expand=True)

#genotype calling in zarr
gt_zarr = callset['calldata/GT']
# print(gt_zarr.info)
# print(gt_zarr)

pos = allel.SortedIndex(callset['variants/POS'])
#this show us the chromosome position of the variant eg. the fist number in the list is the position of the first variant from the vcf file
# print(pos)

#The loc_region object is a slice, which is simply a way of storing a start and stop index. Here the start index is index value for the 
#first variant in the position interval and the stop index is the last variant
loc_region = pos.locate_range(29000000, 33500000)

#to get genotype call from our samples from the variants in our range 
gt_region = allel.GenotypeArray(gt_zarr[loc_region])
# print(gt_region)

#this is the position of all the variants in the loc_region
position_of_variants = pos[loc_region]

#heterzygote call from the region 
heterozygotes = allel.heterozygosity_observed(gt_region)
# print(heterozygotes)

###########Now i would like to filter out the bad variants. #######

#first i make the quality call#
allel_call = callset['/variants/AN'][:]

freq_depth = callset['/variants/AC'][:]
# print(freq_depth)

#now i find the quality for the variants i my locate_range and sort them after threshold
allel_call_region = allel_call[loc_region] > 1400
freq_of_region = freq_depth[loc_region] > 1

def freq_reg(freq_of_region):
	ny_liste = []
	for idx in range(len(freq_of_region)):
		if freq_of_region[idx][0] == True:
			ny_liste.append(True)
		else:
			ny_liste.append(False)
	return ny_liste

frequ = freq_reg(freq_of_region)

# and now i filter them after, so only variants with 1400 called alleles gets through
variant_selection = allel_call_region & frequ
# print(variant_selection)

#this command removes all the variants which are FALSE in the variant_selection from the gt_call
gt_variant_selection = gt_region.compress(variant_selection, axis=0)

# print(gt_variant_selection)

allel_counts = gt_variant_selection.count_alleles()
allel_frequencies = gt_variant_selection.count_alleles().to_frequencies()

modified_count = allel_counts.tolist()
modified_freq = allel_frequencies.tolist()
#now i want to remove the position of all the variants which are FALSE in the variant_selection from the posistion call, so i can align
#pos and gt_call later 
pos_variant_selection = position_of_variants.compress(variant_selection, axis=0)
# print(pos_variant_selection)

observed_heterozygosity= allel.heterozygosity_observed(gt_variant_selection)
expected_heterozygosity = allel.heterozygosity_expected(allel_frequencies, ploidy=2)

####################################### the big file of location of samples
meta_data = 'hgdp_wgs.20190516.metadata.csv'
# and i want to choose on these three criteria 
panel = pd.read_csv(meta_data, usecols=['sample', 'population', 'region'])
show_panel = panel.head()
# print(show_panel)

#this is how i group by population
pop = panel.groupby(by=('population')).count()

# i find the samples from MY DATASET
samples = callset['samples'][:]

#Because the order in the two files are different i have to make a list show the place as index
samples_list = list(samples)
# print(samples_list)
# and the i add a new columm in the metafile with the index number for in the callset
samples_callset_index = [samples_list.index(s) for s in panel['sample']]
panel['callset_index'] = samples_callset_index
index_call = panel.head()
# print(index_call)

##this is how i find all the samples from my data which are from the population under interest. One of the 7
loc_samples_fr = panel[panel.population == 'Kalash'].callset_index.values
# print(loc_samples_fr)

#######Putting together the filters##
gt_fr = gt_variant_selection.take(loc_samples_fr, axis=1)

# print(gt_fr)
#Allele count
fr_allel_counts = gt_fr.count_alleles()

#frequency
fr_allel_frequencies = gt_fr.count_alleles().to_frequencies()

####calculating observed heterzygosity
fr_ob_heterozygosity= allel.heterozygosity_observed(gt_fr)
#And expected
fr_ex_heterozygosity = allel.heterozygosity_expected(fr_allel_frequencies, ploidy=2)

#calcualting Tajimas d
tajimas_d = allel.windowed_tajima_d(pos_variant_selection, fr_allel_counts, size=10000, step=1000)

td_value = tajimas_d[0]
td_window = tajimas_d[1]
td_variants  =tajimas_d[2]
#to csv file

print('t_d','window start','window stop','number of variants', sep = ',')
for idx in range(len(td_value)):
	print(td_value[idx],
		td_window[idx][0],
		td_window[idx][1],
		td_variants[idx],
		sep = ',')