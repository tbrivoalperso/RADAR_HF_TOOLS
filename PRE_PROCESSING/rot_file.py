#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 13:04:08 2019

@author: tbrivoal
"""
import pandas
import os
import xarray as xr
import numpy as np
import scipy.stats as st
import multiprocessing
import time
import glob
from pathlib import Path
import netCDF4 as nc

fldr_U = '/scratch/work/brivoalt/RUNS_NEMO/eNEATL36_AGRIF_vvl/EXP02_AGRIF_finaldomain_bathycorrected_qco_boost2_noslip/DATA_IMMERSE/NC4/eNEATL36_1h_gridU/'
fldr_V = '/scratch/work/brivoalt/RUNS_NEMO/eNEATL36_AGRIF_vvl/EXP02_AGRIF_finaldomain_bathycorrected_qco_boost2_noslip/DATA_IMMERSE/NC4/eNEATL36_1h_gridV/'
start_date = '20170104'
end_date = '20180705'
def rot_file(datein):
    lon_lat = xr.open_dataset('/scratch/work/brivoalt/RUNS_NEMO/eNEATL36_AGRIF_vvl/EXP02_AGRIF_finaldomain_bathycorrected_qco_boost2_noslip//domain_cfg.nc', drop_variables={"x", "y",})
    out_fldr = '/scratch/work/brivoalt/DATA/RADAR_HF/AGRIF_rot/' 
    lont = lon_lat['glamt'].squeeze()
    latt = lon_lat['gphit'].squeeze()
    lonv = lon_lat['glamv'].squeeze()
    latv = lon_lat['gphiv'].squeeze()    
    e1u = lon_lat['e1u'].squeeze()
    e2v = lon_lat['e2v'].squeeze()
    e1t = lon_lat['e1t'].squeeze()
    e2t = lon_lat['e2t'].squeeze()

    # open files and read U & V
    ds_U = xr.open_dataset(fldr_U+'eNEATL36_1h_gridU_'+str(datein)+'-'+str(datein)+'.nc')
    ds_V = xr.open_dataset(fldr_V+'eNEATL36_1h_gridV_'+str(datein)+'-'+str(datein)+'.nc')
    U_noproj_u = ds_U.sozocrtx * e1u 
    V_noproj_v = ds_V.somecrty * e2v 
    U_noproj = U_noproj_u.copy() 
    V_noproj = V_noproj_v.copy()
    U_noproj[:,:,:] = np.nan
    V_noproj[:,:,:] = np.nan
  
    # move U and V to gridT 
    U_noproj[:,:,1::] = (U_noproj_u[:,:,0:-1] + U_noproj_u[:,:,1::])/(e1t[:,0:-1] + e1t[:,1::]) 
    V_noproj[:,1::,:] = (V_noproj_v[:,0:-1,:] + V_noproj_v[:,1::,:])/(e2t[0:-1,:] + e2t[1::,:])
    angle = np.arctan2(((lonv-lont)*np.cos(np.radians(latv))),(latv - latt)) 
   
    # Rotate U & V
    # Nb: rotated U & V will be located at gridT
    U = U_noproj * np.cos(angle) - V_noproj * np.sin(angle)
    V = V_noproj * np.cos(angle) + U_noproj * np.sin(angle)
    U.name = ds_U.sozocrtx.name
    V.name = ds_V.somecrty.name
    U = U.assign_coords(ds_U.sozocrtx.coords)
    V = V.assign_coords(ds_V.somecrty.coords) 

    # Save to a new netcdf file
    U.to_netcdf(str(out_fldr) + 'eNEATL36_1h_gridU_rot_'+str(datein)+'-'+str(datein)+'.nc')
    V.to_netcdf(str(out_fldr) + 'eNEATL36_1h_gridV_rot_'+str(datein)+'-'+str(datein)+'.nc')


# Define a list of dates
datelist = [d.strftime('%Y%m%d') for d in pandas.date_range(start_date,end_date)]
print(datelist)

Nparallel = 64 # Number of (daily) files to treat in parallel
pool = multiprocessing.Pool(processes=Nparallel)

#Launch function filter_file in parallel
a=pool.map(rot_file,datelist)
pool.close()

