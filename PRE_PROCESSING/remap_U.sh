#!/bin/sh
#SBATCH -J RMP
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --no-requeue
#SBATCH --time=20:00:00
#SBATCH --account=cmems

# The radar data must be stored such as in this example:
# => for radar IROI:
# ${DATADIR]/IROI/2017/IROI_2017??.nc 
# ${DATADIR]/IROI/2018/IROI_2018??.nc 
# etc.......

# This script does the following operations:
# - Browse through all radars in $DATADIR
# - Select the first .nc radar file in $DATADIR/RADAR_NAME/201?/ => store it in a gridfile_U.nc file for cdo
# 	/!\ The gridfile_U.nc must have a well defined grid for cdo !!! 
# 	To make sure the grid is well defined, type 'cdo sinfo gridfile_U.nc'. the 'Grid coordinates' must be 'lonlat' or 'curvilinear'
# - browse through nemo files in INDIR, and remap each files towards gridfile_U.nc with cdo remapbil 

DATADIR=/scratch/work/brivoalt/DATA/RADAR_HF/DATA_RADAR/ # Folder where the radar data is stored. They must be stored in 
FLDR_in='/scratch/work/brivoalt/RUNS_NEMO/eNEATL36_AGRIF_vvl/EXP02_AGRIF_finaldomain_bathycorrected_qco_boost2_noslip/DATA_IMMERSE/NC4/'
INDIR=/scratch/work/brivoalt/DATA/RADAR_HF/AGRIF_rot/ # Folder where the rotated nemo files are stored
cd $INDIR
pwd
list_grid_U=$( ls *gridU* )


cd ${DATADIR}
echo 'here'
pwd
liste_radar=$( ls --ignore="*rot*" -d */ ) 
echo $liste_radar
for RADAR in $liste_radar ;
do
   if [ "$RADAR" == "AGRIF_rot/" ] || [ "$RADAR" == "TWIN_rot/" ];
   then 
      echo "test"
   else
      pwd
      echo $RADAR
      grid_out=$( ls ${RADAR}/2*/*.nc | head -1 )
      echo $grid_out
      rm gridfile_U.nc
      cdo selvar,EWCT $grid_out gridfile_U.nc
      if [ -f gridfile_U.nc ] 
      then
          echo 'grid file done'
      else
          echo 'EWCT is not in data, tryin u instead'
          cdo selvar,u $grid_out gridfile_U.nc
   
          if [ -f gridfile_U.nc ]
          then
              echo 'grid file done'
          else
              echo 'u is not in data, tryin water_u instead'
              cdo selvar,water_u $grid_out gridfile_U.nc
    
          fi
      fi
   fi 

   
   echo $grid_out
   mkdir ${RADAR}/EQUIVALENT_MODELE

   for file in $list_grid_U   
   do
   if [ "$RADAR" == "AGRIF_rot/" ] || [ "$RADAR" == "TWIN_rot/" ];
   then
      echo "test"
   else

       if [ -f ${RADAR}/EQUIVALENT_MODELE/${file} ]
       then 
       echo "${RADAR}/EQUIVALENT_MODELE/${file} DONE"
       else 
      # echo cdo remapbil,gridfile_U.nc ${INDIR}/${file} ${RADAR}/EQUIVALENT_MODELE/${file} 
       cdo remapbil,gridfile_U.nc ${INDIR}/${file} ${RADAR}/EQUIVALENT_MODELE/${file}
       fi
   fi
   done
   echo $RADAR ' done'
done


