# RADAR_HF_TOOLS

This repository regroup :
- Some tools to rotate, and interpolate modelled NEMO current vectors on the HF radar grids
- Some python scripts to plot:
  - [Maps of KE, with the current vectors overlayed](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/MAP_KE_and_current_vectors_RADAR_HF_per_season.py)
  - [Maps of the RMS with the observed currents](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/MAP_U_RMS_RADAR_HF_per_season.py)
  - [Maps of the correlation with the observations](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/MAP_U_correl_RADAR_HF_per_season.py)
  - [Taylor diagrams](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/RADAR_HF_taylor_diagram_comparaison.py)
  - [Current roses](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/RADAR_HF_currentrose_comparaison.py)



### Requirements to use the scripts

For the scripts to work:
- Make sure you have cdo installed, and the python libraries [skill_metrics](https://github.com/PeterRochford/SkillMetrics) (for the taylor diagrams) and [windrose](https://windrose.readthedocs.io/en/latest/install.html)
- The HF RADAR files must be gridded in a lon/lat grid and their grid format must be compatible with cdo 
=> To make sure the grid is well defined, type 'cdo sinfo gridfile_U.nc'. the 'Grid coordinates' must be 'lonlat' or 'curvilinear'
- The radar data must be stored in 1 directory per radar and per year such as in this example:
- The HF RADAR files also must have a time variable that can be readable by xarray
```
=> for radar IROI:

 ${DATADIR]/IROI/2017/IROI_2017??.nc
 
 ${DATADIR]/IROI/2018/IROI_2018??.nc
 etc.......
```

### How to use the scripts:

- First, change paths in [rot_file.py](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PRE_PROCESSING/rot_file.py) and use this script (`python rot_file.py`) to rotate U & V towards a regular North - South - East - West cartesian coordinate system
- Then change paths in [remap_U.sh](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PRE_PROCESSING/remap_U.sh) and [remap_V.sh](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PRE_PROCESSING/remap_U.sh) and use these scripts (e.g : `./remap_U.sh`) to remap the model data towards the radar HF grids
- At this point, in each radar folder, you will have a new folder created (with the name of the folder given by the variable EQUIVALENT_MODELE_DIR in remap_U.sh and remap_V.sh scripts with the remapped model files inside
- Note that the remapped files are only remapped spatially. The temporal interpolation is made in each python scripts by the xarray interp function

 
