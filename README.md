# DRM and 2L-DRM 
This repository includes source codes for the DRM and 2-layer DRM (2L-DRM) models. The following explains the input files required to run the models and the output files from the models.

## DRM: (drm_WRF25NAM5.f90)
Reference: [Dominguez et al. (2006)](https://journals.ametsoc.org/doi/10.1175/JCLI3691.1) and [Martinez and Dominguez (2014)](https://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00022.1)

![N|Solid](https://github.com/huancui/DRM_2LDRM/blob/master/Picture1.png)

### Input files:
1.	Forcing data: 

|  | File names | Description |
| ------ | ------ | ------ |
| A	| DRMyyyy0601WRFddd_ET.txt | evapotranspiration (ET) with daily resolution |
| B | DRMyyyy0601WRFddd_PP.txt | precipitation (PP) with daily resolution* |
| C | DRMyyyy0601WRFddd_PW.txt | precipitable water (PW) with daily resolution |
| D | DRMyyyy0601WRFddd_U3.txt | column vapor weighted U with 3-hourly resolution |
| E | DRMyyyy0601WRFddd_V3.txt | column vapor weighted V with 3-hourly resolution |

(*Precipitation is used to weight regional contribution, is not necessary if replaced by PW to weight)

2.	Region setup files (both can be found in this repository as examples):


|  | File names | Description |
| ------ | ------ | ------ |
| A | WRF_Tracer_NAM_ij.txt | list of i and j indices for target region grids in the domain |
| B | WRF_Tracer_5reg.txt | regional masks for each sub-region defined in the domain (specify each sub-region with a different value) | 


### Output files:

|  | File names | description |
| ------ | ------| ------ |
| A* | WRF_OUTPUT_yyyy_GRDR_PW.txt | currently not used |
| B* | WRF_OUPTUT_yyyy_DAYR_PW.txt | ratio contributions from each sub-region to target region |
| C* | WRF_OUTPUT_yyyy_BCTJ_PW.txt | recording information along each back-trajectory (diagnosing purpose) |

(*"_PW" in output files indicate regional contributions weighted by precipitable water)

---

## 2L-DRM: (drm2L_WRF25NAM5_PW.f90)
Reference: [Dominguez et al. (2020)](https://journals.ametsoc.org/view/journals/hydr/21/1/jhm-d-19-0101.1.xml)

![N|Solid](https://github.com/huancui/DRM_2LDRM/blob/master/Picture2.png)

### Input files:
1.	Forcing data:

|  | File names | Description |
| ------ | ------ | ------ |
| A | DRMyyyy0601WRFddd_ET.txt | evapotranspiration (ET) with daily resolution |
| B | DRMyyyy0601WRFddd_PP.txt | precipitation (PP) with daily resolution* |
| C | DRMyyyy0601WRFddd_PWU.txt | precipitable water for the upper slab (upper-slab integrated vapor) with daily resolution | 
| D | DRMyyyy0601WRFddd_PWL.txt | precipitable water for the lower slab (lower-slab integrated vapor) with daily resolution | 
| E | DRMyyyy0601WRFddd_U3U.txt | upper-slab vapor weighted U with 3-hourly resolution |
| F | DRMyyyy0601WRFddd_U3L.txt | lower-slab vapor weighted U with 3-hourly resolution |
| G | DRMyyyy0601WRFddd_V3U.txt | upper-slab vapor weighted V with 3-hourly resolution |
| H | DRMyyyy0601WRFddd_V3L.txt | lower-slab vapor weighted V with 3-hourly resolution |
| I | DRMyyyy0601WRFddd_W3U.txt | upward vertical wind W at 700mb level with 3-hourly resolution |
| J | DRMyyyy0601WRFddd_W3L.txt | downward vertical wind W at 700mb level with 3-hourly resolution |
| K | DRMyyyy0601WRFddd_Q70.txt | humidity at 700mb level with 3-hourly resolution |

(*Precipitation is used to weight regional contribution, is not necessary if replaced by PW to weight)

2.	Region setup files:

|  | File names | description |
| ------ | ------ | ------ |
| A | WRF_Tracer_NAM_ij.txt | list of i and j indices for target region grids in the domain |
| B | WRF_Tracer_5reg.txt | regional masks for each sub-region defined in the domain (specify each sub-region with a different value) |
| C | 700mb_mask.txt | masks to indicate if the grid has only the upper slab (elevation is above 700mb) or both slabs (0 – only upper slab, 1 – both slabs) | 


### Output files:

|  | File names | Description |
| ------ | ------ | ------ |
| A* | WRF_OUTPUT_2L5R_yyyy_DAYRA_PW.txt | ratio contributions from each sub-region to target region averaged from upper and lower slabs | 
| B* | WRF_OUTPUT_2L5R_yyyy_DAYRU_PW.txt | ratio contributions from each sub-region to target region only for the upper slab |
| C* | WRF_OUTPUT_2L5R_yyyy_DAYRL_PW.txt | ratio contributions from each sub-region to target region only for the lower slab |

(*"_PW" in output files indicates regional contributions weighted by precipitable water)

