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
# - Select the first .nc radar file in $DATADIR/RADAR_NAME/201?/ => store it in a gridfile_V.nc file for cdo
#       /!\ The gridfile_V.nc must have a well defined grid for cdo !!! 
#       To make sure the grid is well defined, type 'cdo sinfo gridfile_V.nc'. the 'Grid coordinates' must be 'lonlat' or 'curvilinear'
# - browse through nemo files in INDIR, and remap each files towards gridfile_V.nc with cdo remapbil 

DATADIR=/scratch/work/brivoalt/DATA/RADAR_HF/DATA_RADAR/ # Folder where the radar data is stored. They must be stored in 
FLDR_in='/scratch/work/brivoalt/RUNS_NEMO/eNEATL36_AGRIF_vvl/EXP02_AGRIF_finaldomain_bathycorrected_qco_boost2_noslip/DATA_IMMERSE/NC4/'
INDIR=/scratch/work/brivoalt/DATA/RADAR_HF/AGRIF_rot/ # Folder where nemo files are stored
EQUIVALENT_MODELE_DIR=EQUIVALENT_MODELE # Name of the equivalent model folder to be created


cd $INDIR
pwd
list_grid_V=$( ls *gridV* )


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
      rm gridfile_V_AGRIF.nc
      echo $grid_out
      cdo selvar,NSCT $grid_out gridfile_V_AGRIF.nc
      if [ -f gridfile_V_AGRIF.nc ] 
      then
          echo 'grid file done'
      else
          echo 'NSCT is not in data, tryin u instead'
          cdo selvar,v $grid_out gridfile_V_AGRIF.nc
   
          if [ -f gridfile_V_AGRIF.nc ]
          then
              echo 'grid file done'
          else
              echo 'u is not in data, tryin water_v instead'
              cdo selvar,water_v $grid_out gridfile_V_AGRIF.nc
    
          fi
      fi
   fi 

   
   echo $grid_out
   mkdir ${RADAR}/${EQUIVALENT_MODELE_DIR}

   for file in $list_grid_V   
   do
   if [ "$RADAR" == "AGRIF_rot/" ] || [ "$RADAR" == "TWIN_rot/" ];
   then
      echo "test"
   else

       if [ -f ${RADAR}/${EQUIVALENT_MODELE_DIR}/${file} ]
       then 
       echo "${RADAR}/${EQUIVALENT_MODELE_DIR}/${file} DONE"
       else 
#       echo cdo remapbil,gridfile_V_AGRIF.nc ${INDIR}/${file} ${RADAR}/${EQUIVALENT_MODELE_DIR}/${file} 
       cdo remapbil,gridfile_V_AGRIF.nc ${INDIR}/${file} ${RADAR}/${EQUIVALENT_MODELE_DIR}/${file}
       fi
   fi
   done
   echo $RADAR ' done'

done


