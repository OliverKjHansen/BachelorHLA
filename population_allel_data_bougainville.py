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
vcf_path = 'hgdp_wgs.20190516.full.chr6.vcf.gz'
zarr_path = 'doing_zarr.zarr'
# allel.vcf_to_zarr(vcf_path, zarr_path, fields='*', log=sys.stdout, overwrite=True)


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
#so now we look for which variants, and how many, lay within the 29940000-33100000 range on the chromosome
# print(loc_region)

#to get genotype call from our samples from the variants in our range 
gt_region = allel.GenotypeArray(gt_zarr[loc_region])
# print(gt_region)

#this is the position of all the variants in the loc_region
position_of_variants = pos[loc_region]

#heterzygote call from the region 
heterozygotes = allel.heterozygosity_observed(gt_region)
# print(heterozygotes)

###########Now i would like to filter out the bad variants. the threshold value is QD<2#######

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

# and now i filter them after, so only variants with 1400 called alleles an more then 1 alternative allele per variant gets through
variant_selection = allel_call_region & frequ
# print(variant_selection)

#this command removes all the variants which are FALSE in the variant_selection from the gt_call
gt_variant_selection = gt_region.compress(variant_selection, axis=0)

#now i want to remove the position of all the variants which are FALSE in the variant_selection from the posistion call, so i can align
#pos and gt_call later 
pos_variant_selection = position_of_variants.compress(variant_selection, axis=0)

############now i would like to be able to choose specific samples##

# the big file of location of sample
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

##this is how i find all the samples from my data which are Bedouin in the meta file. I write the name of the population i want to parse from in the brackets 
loc_samples_fr = panel[panel.population == 'Bougainville'].callset_index.values
# print(loc_samples_fr)

#######Putting together the filters into a genotype array over our population ###
gt_fr = gt_variant_selection.take(loc_samples_fr, axis=1)

#counting allels and frequencies 
allel_counts = gt_fr.count_alleles()
allel_frequencies = gt_fr.count_alleles().to_frequencies()

#making it easier to make into a csv file 
modified_count = allel_counts.tolist()
modified_freq = allel_frequencies.tolist()

####calculating observed heterzygosity#####
ob_heterozygosity= allel.heterozygosity_observed(gt_fr)

###expected heterzygosity##But im not going to use this 
ex_heterozygosity = allel.heterozygosity_expected(allel_frequencies, ploidy=2)

print('Position','REF','ALT1','ALT2','ALT3','ALT4','ALT5','ALT6','AN','AC','Fref','F1','F2','F3','F4','F5','F6','1','Hob','Hex', sep = ',')
for idx in range(len(pos_variant_selection)):
  print(pos_variant_selection[idx],
  	modified_count[idx][0],
  	modified_count[idx][1],
  	modified_count[idx][2],
  	modified_count[idx][3],
  	modified_count[idx][4],
  	modified_count[idx][5],
  	modified_count[idx][6],
  	(modified_count[idx][0]+ modified_count[idx][1]+modified_count[idx][2]+ modified_count[idx][3] + modified_count[idx][4]+modified_count[idx][5]+modified_count[idx][6]),
	(modified_count[idx][1]+modified_count[idx][2]+ modified_count[idx][3] + modified_count[idx][4]+modified_count[idx][5]+modified_count[idx][6]),
  	modified_freq[idx][0],
  	modified_freq[idx][1],
  	modified_freq[idx][2],
  	modified_freq[idx][3],
  	modified_freq[idx][4],
  	modified_freq[idx][5],
  	modified_freq[idx][6],
  	(modified_freq[idx][0]+modified_freq[idx][1]+modified_freq[idx][2]+modified_freq[idx][3]+modified_freq[idx][4]+modified_freq[idx][5]+modified_freq[idx][6]),
  	ob_heterozygosity[idx],
  	ex_heterozygosity[idx],
  	sep = ',')