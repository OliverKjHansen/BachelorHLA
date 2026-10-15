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
loc_region = pos.locate_range(91500000, 96000000)
#so now we look for which variants, and how many, lay within the 29000000-33500000 range on the chromosome
# for the control region 91500000-96000000 is used instead
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


#now i find the quality for the variants i my locate_range and sort them after threshold
allel_call_region = allel_call[loc_region] > 1400

# and now i filter them after, so only variants with 1400 called alleles gets through
variant_selection = allel_call_region

#this command removes all the variants which are FALSE in the variant_selection from the gt_call
gt_variant_selection = gt_region.compress(variant_selection, axis=0)

#function for counting the number of allels from the filtered gt call
allel_counts = gt_variant_selection.count_alleles()
# print(allel_counts)

#so the data is easier to make into a csv file
modified_count = allel_counts.tolist()

#now i want to remove the position of all the variants which are FALSE in the variant_selection from the posistion call, so i can align
#pos and gt_call later 
pos_variant_selection = position_of_variants.compress(variant_selection, axis=0)

print('Position','REF','ALT1','ALT2','ALT3','ALT4','ALT5','ALT6','SUM_of_ALT', sep = ',')
for idx in range(len(pos_variant_selection)):
  print(pos_variant_selection[idx],modified_count[idx][0],modified_count[idx][1],modified_count[idx][2],modified_count[idx][3],
  modified_count[idx][4],modified_count[idx][5],modified_count[idx][6],
  (modified_count[idx][1]+modified_count[idx][2]+modified_count[idx][3]+modified_count[idx][4]+modified_count[idx][5]+modified_count[idx][6]), sep = ',')