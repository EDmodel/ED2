#!/bin/bash

# use an ref lai of 4 for pheno-plasticity
# reduce stem and root respiration

# deal with potential argument
# we need to store whether it is a test
is_submit=false
for i in "$*"
do
    if [ "$i" == "-s" ]; then
        is_submit=true
    fi
done

# test the new PFT setup
# constant CO2
# run for 150 years

# some constant
rean_res='0.5' # half degree resolution

# first site information, comes from MLONGO's dissertation
# from dry to wet
#site_array=(    'PNZ'    'BSB'    'PDG'    'S67'    'RJA'    'CAX'     'GYF'    'M34'    )
#lon_array=(     '-40.37' '-47.71' '-47.65' '-54.96' '-61.93' '-51.46'  '-52.91' '-60.21' )
#lat_array=(     '-9.17'  '-15.60' '-21.62' '-2.86'  '-10.08' '-1.71'   '5.28'   '-2.61'  )
#sand_array=(    '0.82'   '0.13'   '0.85'   '0.39'   '0.8'    '0.78'    '0.56'   '0.2'    )
#clay_array=(    '0.05'   '0.53'   '0.03'   '0.59'   '0.1'    '0.15'    '0.34'   '0.68'   )
#METCYC1_array=( '2004'   '2010'   '2001'   '2001'   '1999'   '1999'    '2004'   '1999'   )
#METCYCF_array=( '2012'   '2012'   '2003'   '2011'   '2002'   '2008'    '2014'   '2006'   )
#soil_array=(    1        1        1        0        1        0         0        0        )

site_array=(      'M34'    )
lon_array=(       '-60.21' )
lat_array=(       '-2.61'  )
sand_array=(      '0.2'    )
clay_array=(      '0.68'   )
METCYC1_array=(   '1999'   )
METCYCF_array=(   '2006'   )
soil_array=(      0        )

rean_METCYC1=1981
rean_METCYCF=2000 # do not use drought years

# soil info
# 0 -> 10 meter
# 1 -> 3 meter
NZG_array=(16 16)
SLZ_array=(
 '-10.00,-9.00,-8.00,-7.00,-6.00,-5.00,-4.00,-3.00,-2.50,-2.00,-1.50,-1.00,-0.60,-0.40,-0.20,-0.10'
 '-3.00,-2.70,-2.40,-2.10,-1.80,-1.50,-1.20,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10'
 )
SLMSTR_array=(
 ' 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,0.90, 0.90, 0.90'
 ' 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90,0.90, 0.90, 0.90'
)
STGOFF_array=(
 ' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
 ' 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
)
#SOILBC_array=(1 1 0)
SOILBC_array=(1 1 1)

# setup array
setup_array=()

# two kind of simulations
# include or not include hydraulic mortality
setup_array=()
for hm in {1..1}; do
    setup_array+=("test")
done

# hr is hydraulic redistribution
# hm is hydraulic failure mortality

site_num=${#site_array[@]}
setup_num=${#setup_array[@]}


# some constants
ED2IN_template="./ED2IN"
xml_template="../template_config.xml"
ed2_exec="./ed2"

queue_name='shared'
core_num=24
sim_time='120:00:00' # 5 days
# loop over each setup
for (( isetup=0; isetup<${setup_num}; isetup++ )); do
    setup=${setup_array[$isetup]}
    # first whether to include hydraulic redistribution

    # second whether to include hydro mortality

    # loop over site

    for (( isite=0; isite<${site_num}; isite++ ));
    do
        site=${site_array[$isite]}
        lon=${lon_array[$isite]}
        lat=${lat_array[$isite]}
        sand=${sand_array[$isite]}
        clay=${clay_array[$isite]}
        NZG=${NZG_array[${soil_array[$isite]}]}
        SLZ=${SLZ_array[${soil_array[$isite]}]}
        SLMSTR=${SLMSTR_array[${soil_array[$isite]}]}
        STGOFF=${STGOFF_array[${soil_array[$isite]}]}
        SOILBC=${SOILBC_array[${soil_array[$isite]}]}

        #prefix of this simulation
        sim_pf="${site}_${setup}"
        job_sf="${sim_pf}"
        #create ED2IN file
        ED2IN_fn="ED2IN_${sim_pf}"
        cp $ED2IN_template $ED2IN_fn
        #create xml file
#        xml_fn="config_${sim_pf}.xml"
#        cp $xml_template $xml_fn

        
        # define ED2IN flags to change
        YEARA=1980
        YEARZ=1980
        RUNTYPE="'INITIAL'"
        IED_INIT_MODE=0

        # first option meteorology
        MET_HEADER="\'\/n\/moorcroftfs5\/xiangtao\/ED_MET\/MERRA2_CHIRPS_sites\/${site}\/REAN_GRND_MET_HEADER\'"
        METCYC1=${rean_METCYC1}
        METCYCF=${rean_METCYCF}


        # other physiological setups
        IPHYSIOL=2
        IALLOM=4 # max hite is 45m
        PLANT_HYDRO_SCHEME=1 # BC paper
        CARBON_MORTALITY_SCHEME=2
        HYDRAULIC_MORTALITY_SCHEME=1
        H2O_PLANT_LIM=4      #
        ISTOMATA_SCHEME=1    #
        ISTRUCT_GROWTH_SCHEME=2 # repro increase with height
        ISTEM_RESPIRATION_SCHEME=1 # new stem respiraition with growth resp
        GROWTHRESP='0.35' # growth respiration fraction
        INCLUDE_THESE_PFT='1,2,3,4'
        IPHEN_SCHEME=4
        TRAIT_PLASTICITY_SCHEME=3

        ################################


        # modify ED2IN
        # runtype
        sed -i "s/   NL%RUNTYPE.*/   NL%RUNTYPE   = ${RUNTYPE}/g" $ED2IN_fn
        sed -i "s/.*NL%IED_INIT_MODE.*/   NL%IED_INIT_MODE   = ${IED_INIT_MODE}/g" $ED2IN_fn

        # Sim time
        sed -i "s/.*NL%IMONTHA.*/   NL%IMONTHA   = 07/g" $ED2IN_fn
        sed -i "s/.*NL%IDATEA.*/   NL%IDATEA   = 01/g" $ED2IN_fn
        sed -i "s/.*NL%IYEARA.*/   NL%IYEARA   = ${YEARA}/g" $ED2IN_fn
        sed -i "s/.*NL%ITIMEA.*/   NL%ITIMEA   = 0000/g" $ED2IN_fn
        sed -i "s/.*NL%IMONTHZ.*/   NL%IMONTHZ   = 07/g" $ED2IN_fn
        sed -i "s/.*NL%IDATEZ.*/   NL%IDATEZ   = 01/g" $ED2IN_fn
        sed -i "s/.*NL%IYEARZ.*/   NL%IYEARZ   = ${YEARZ}/g" $ED2IN_fn
        sed -i "s/.*NL%ITIMEZ.*/   NL%ITIMEZ   = 0000/g" $ED2IN_fn
        sed -i "s/.*NL%IMONTHH.*/   NL%IMONTHH   = 09/g" $ED2IN_fn
        sed -i "s/.*NL%IDATEH.*/   NL%IDATEH   = 01/g" $ED2IN_fn
        sed -i "s/.*NL%IYEARH.*/   NL%IYEARH   = ${YEARA}/g" $ED2IN_fn
        sed -i "s/.*NL%ITIMEH.*/   NL%ITIMEH   = 0000/g" $ED2IN_fn

        # Geographical Coordinates
        sed -i "s/.*NL%POI_LAT.*/   NL%POI_LAT   = ${lat}/g" $ED2IN_fn
        sed -i "s/.*NL%POI_LON.*/   NL%POI_LON   = ${lon}/g" $ED2IN_fn

        # I/O
        sed -i "s/.*NL%IYOUTPUT.*/   NL%IYOUTPUT   = 3/g" $ED2IN_fn
        sed -i "s/.*NL%IMOUTPUT.*/   NL%IMOUTPUT   = 3/g" $ED2IN_fn
        sed -i "s/.*NL%IQOUTPUT.*/   NL%IQOUTPUT   = 3/g" $ED2IN_fn
        sed -i "s/.*NL%ISOUTPUT.*/   NL%ISOUTPUT   = 3/g" $ED2IN_fn
        sed -i "s/.*NL%FRQSTATE.*/   NL%FRQSTATE   = 5/g" $ED2IN_fn
        sed -i "s/.*NL%IADD_PATCH_MEANS.*/   NL%IADD_PATCH_MEANS  = 1/g" $ED2IN_fn
        sed -i "s/.*NL%IADD_COHORT_MEANS.*/   NL%IADD_COHORT_MEANS = 1/g" $ED2IN_fn
        sed -i "s/.*NL%IADD_SITE_MEANS.*/   NL%IADD_SITE_MEANS  = 1/g" $ED2IN_fn

        output_sed=".\/test"
        sed -i "s/^   NL%FFILOUT.*/   NL%FFILOUT = \'${output_sed}\'/g" $ED2IN_fn
        sed -i "s/^   NL%SFILOUT.*/   NL%SFILOUT = \'${output_sed}\'/g" $ED2IN_fn
        sed -i "s/^   NL%SFILIN.*/   NL%SFILIN = \'${output_sed}\'/g" $ED2IN_fn

        # physiology and water stress
        sed -i "s/.*NL%IPHYSIOL.*/   NL%IPHYSIOL = ${IPHYSIOL}/g" $ED2IN_fn
        sed -i "s/.*NL%IPHEN_SCHEME.*/   NL%IPHEN_SCHEME = ${IPHEN_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%H2O_PLANT_LIM.*/   NL%H2O_PLANT_LIM = ${H2O_PLANT_LIM}/g" $ED2IN_fn
        sed -i "s/.*NL%PLANT_HYDRO_SCHEME.*/   NL%PLANT_HYDRO_SCHEME = ${PLANT_HYDRO_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%ISTRUCT_GROWTH_SCHEME.*/   NL%ISTRUCT_GROWTH_SCHEME = ${ISTRUCT_GROWTH_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%ISTEM_RESPIRATION_SCHEME.*/   NL%ISTEM_RESPIRATION_SCHEME = ${ISTEM_RESPIRATION_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%GROWTHRESP.*/   NL%GROWTHRESP = ${GROWTHRESP}/g" $ED2IN_fn
        sed -i "s/.*NL%ISTOMATA_SCHEME.*/   NL%ISTOMATA_SCHEME = ${ISTOMATA_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%TRAIT_PLASTICITY_SCHEME.*/   NL%TRAIT_PLASTICITY_SCHEME = ${TRAIT_PLASTICITY_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%INCLUDE_THESE_PFT.*/   NL%INCLUDE_THESE_PFT = ${INCLUDE_THESE_PFT}/g" $ED2IN_fn
        sed -i "s/.*NL%CARBON_MORTALITY_SCHEME.*/   NL%CARBON_MORTALITY_SCHEME = ${CARBON_MORTALITY_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%HYDRAULIC_MORTALITY_SCHEME.*/   NL%HYDRAULIC_MORTALITY_SCHEME = ${HYDRAULIC_MORTALITY_SCHEME}/g" $ED2IN_fn
        sed -i "s/.*NL%IALLOM.*/   NL%IALLOM = ${IALLOM}/g" $ED2IN_fn
        #sed -i "s/.*NL%IEDCNFGF.*/   NL%IEDCNFGF = '${xml_fn}'/g" $ED2IN_fn

        # soil
        sed -i "s/.*NL%NZG.*/   NL%NZG = ${NZG}/g" $ED2IN_fn
        sed -i "s/.*NL%SLXCLAY.*/   NL%SLXCLAY = ${clay}/g" $ED2IN_fn
        sed -i "s/.*NL%SLXSAND.*/   NL%SLXSAND = ${sand}/g" $ED2IN_fn
        sed -i "s/.*NL%SLZ.*/   NL%SLZ     = ${SLZ}/g" $ED2IN_fn
        sed -i "s/.*NL%SLMSTR.*/   NL%SLMSTR  = ${SLMSTR}/g" $ED2IN_fn
        sed -i "s/.*NL%STGOFF.*/   NL%STGOFF  = ${STGOFF}/g" $ED2IN_fn
        sed -i "s/.*NL%ISOILBC.*/   NL%ISOILBC = ${SOILBC}/g" $ED2IN_fn

        # Meteorology
        sed -i "s/.*NL%ED_MET_DRIVER_DB.*/   NL%ED_MET_DRIVER_DB  = $MET_HEADER/g" $ED2IN_fn
        sed -i "s/.*NL%IMETAVG.*/   NL%IMETAVG      = 3/g" $ED2IN_fn
        sed -i "s/.*NL%METCYC1.*/   NL%METCYC1      = ${METCYC1}/g" $ED2IN_fn
        sed -i "s/.*NL%INITIAL_CO2.*/   NL%INITIAL_CO2 = 380./g" $ED2IN_fn
        sed -i "s/.*NL%METCYCF.*/   NL%METCYCF      = ${METCYCF}/g" $ED2IN_fn

        # OTHER
        #sed -i "s/.*NL%TREEFALL_DISTURBANCE_RATE.*/   NL%TREEFALL_DISTURBANCE_RATE  = 0.0/g" $ED2IN_fn
        #sed -i "s/.*NL%MAXPATCH.*/   NL%MAXPATCH      = 12/g" $ED2IN_fn



        # submit the job
        cmd_str="ulimit -s unlimited ; export OMP_NUM_THREADS=${core_num} ; cd \$SLURM_SUBMIT_DIR ; srun -n 1 --cpus-per-task=${core_num} ${ed2_exec} -f ${ED2IN_fn}"

        slurm_opts="-o ${job_sf}.out -e ${job_sf}.err -J ${job_sf} -t ${sim_time} --mem-per-cpu=1000 -n 1 -c ${core_num} -p ${queue_name} --mail-type=END --mail-user=xu.withoutwax@gmail.com"
        echo ${slurm_opts}
        echo ${cmd_str}
        
        # submit job
        if [ "$is_submit" == true ]; then
            jid_last=$(sbatch ${slurm_opts} --wrap="${cmd_str}")
            sleep 0.5
            echo ${jid_last}
        else
            echo "Only a test; Job not submitted. Use -s to submit the jobs"
        fi


    done
done

