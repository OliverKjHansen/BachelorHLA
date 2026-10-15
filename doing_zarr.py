
import allel
import sys
import zarr
import numpy as np 
import numcodecs


vcf_path = 'hgdp_wgs.20190516.full.chr6.vcf.gz'
zarr_path = 'doing_zarr.zarr'
allel.vcf_to_zarr(vcf_path, zarr_path, fields='*', log=sys.stdout, overwrite=True)

