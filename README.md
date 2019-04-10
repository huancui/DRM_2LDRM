# DRM and 2L-DRM 
This repository includes source codes for the DRM and 2-layer DRM (2L-DRM) models.

## DRM:
Reference: [Dominguez et al. (2006)](https://journals.ametsoc.org/doi/10.1175/JCLI3691.1) and [Martinez and Dominguez (2014)](https://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00022.1)
#### Input files:
1.	Forcing data: 

|  | File names | Description |
| ------ | ------ | ------ |
| A	| DRMyyyy0601WRFddd_ET.txt | evapotranspiration (ET) at daily resolution |
| B | DRMyyyy0601WRFddd_PP.txt | Precipitation (PP) at daily resolution* |
| C | DRMyyyy0601WRFddd_PW.txt | precipitable water (PW) at daily resolution |
| D | DRMyyyy0601WRFddd_U3.txt | Column vapor weighted U at 3-hourly resolution |
| E | DRMyyyy0601WRFddd_V3.txt | Column vapor weighted V at 3-hourly resolution |

(*Precipitation is used to weight regional contribution, is not necessary if replaced by PW to weight)

2.	Region setup files:


	File names	Description
A	WRF_Tracer_NAM_ij.txt	List of i and j indices for target region grids in the domain
B 	WRF_Tracer_5reg.txt	Regional masks for each sub-region defined in the domain (specify each sub-region with a different value)

Output files:
	File names	description
A 	WRF_OUTPUT_yyyy_GRDR_PW.txt	Currently not used
B	WRF_OUPTUT_yyyy_DAYR_PW.txt	Ratio contributions from each sub-region to target region
C	WRF_OUTPUT_yyyy_BCTJ_PW*.txt	Recording information along each back-trajectory (diagnosing purpose)
      *_PW in output files indicate regional contributions weighted by precipitable water


2L-DRM:
Reference: Dominguez et al. (2019) to be submitted
Input files:
1.	Forcing data:
	File names	Description
A	DRMyyyy0601WRFddd_ET.txt	evapotranspiration (ET) at daily resolution
B	DRMyyyy0601WRFddd_PP.txt	Precipitation (PP) at daily resolution 
(used to weight regional contribution, is not necessary if replaced by PW to weight)
C	DRMyyyy0601WRFddd_PWU.txt	Precipitable water for the upper slab (upper-slab integrated vapor) at daily resolution
D	DRMyyyy0601WRFddd_PWL.txt	Precipitable water for the lower slab (lower-slab integrated vapor) at daily resolution
E	DRMyyyy0601WRFddd_U3U.txt	Upper-slab vapor weighted U at 3-hourly resolution
F	DRMyyyy0601WRFddd_U3L.txt	Lower-slab vapor weighted U at 3-hourly resolution
G	DRMyyyy0601WRFddd_V3U.txt	Upper-slab vapor weighted V at 3-hourly resolution
H 	DRMyyyy0601WRFddd_V3L.txt	Lower-slab vapor weighted V at 3-hourly resolution
I	DRMyyyy0601WRFddd_W3U.txt	Upward vertical wind W at 700mb level at 3-hourly resolution
J 	DRMyyyy0601WRFddd_W3L.txt	Downward vertical wind W at 700mb level at 3-hourly resolution
K 	DRMyyyy0601WRFddd_Q70.txt	Humidity at 700mb level at 3-hourly resolution

2.	Region setup files:
	File names	description
A	WRF_Tracer_NAM_ij.txt	List of i and j indices for target region grids in the domain
B	WRF_Tracer_5reg.txt	Regional masks for each sub-region defined in the domain (specify each sub-region with a different value)
C	700mb_mask.txt	Masks to indicate if the grid has only the upper slab (elevation is above 700mb) or both slabs
(0 – only upper slab, 1 – both slabs)

Output files:
	File names	Description
A	WRF_OUTPUT_2L5R_yyyy_DAYRA_PW.txt	Ratio contributions from each sub-region to target region averaged from upper and lower slabs
B	WRF_OUTPUT_2L5R_yyyy_DAYRU_PW.txt	Ratio contributions from each sub-region to target region only for the upper slab
C	WRF_OUTPUT_2L5R_yyyy_DAYRL_PW*.txt	Ratio contributions from each sub-region to target region only for the lower slab
*_PW in output files indicates regional contributions weighted by precipitable water

