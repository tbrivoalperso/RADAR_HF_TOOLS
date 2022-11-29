# RADAR_HF_TOOLS

This repository regroup :
- Some tools to rotate, and interpolate modelled NEMO current vectors on the HF radar grids
- Some python scripts to plot:
  - [Maps of KE, with the current vectors overlayed](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/MAP_KE_and_current_vectors_RADAR_HF_per_season.py)
  - [Maps of the RMS with the observed currents](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/MAP_U_RMS_RADAR_HF_per_season.py)
  - [Maps of the correlation with the observations](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/MAP_U_correl_RADAR_HF_per_season.py)
  - [Taylor diagrams](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/RADAR_HF_taylor_diagram_comparaison.py)
  - [Current roses](https://github.com/tbrivoalperso/RADAR_HF_TOOLS/blob/master/PLOTS/RADAR_HF_currentrose_comparaison.py)



### How to use the scripts

For the scripts to work:
- Make sure you have cdo installed, and the python libraries [skill_metrics](https://github.com/PeterRochford/SkillMetrics) (for the taylor diagrams) and [windrose](https://windrose.readthedocs.io/en/latest/install.html)
- The HF RADAR files must be gridded in a lon/lat grid (or curvilinear grid) and their grid format must be compatible with cdo
- To make sure the grid is well defined, type 'cdo sinfo gridfile_U.nc'. the 'Grid coordinates' must be 'lonlat' or 'curvilinear'

