#!/bin/bash
#Please the header of this file carefully!
#Running this MWR_PRO script will produce daily netcdf files and plots of
#1.) level1 (measurement products (TBs))
#2.) level2a (integrated retrieval products)
#(needs LWP, IWV, ZWD or ATT retrievals)
#3.) level2b&c (T&q profile retrieval products)
#(needs HUM, TZE or TEL retrievals)
#from raw radiometer data on a monthly basis. You can process (test-wise) only one day by setting the variable "day_files" 
#to one specific day (see below).
#--> Look for ####SPECIFY begin and ####SPECIFY end in this file to specify the dates you would like to process.
#--> You must specify the retrieval products and the rest of processing options
#    in the file par_mwr_pro.pro in this directory. Getting aquainted with par_mwr_pro.pro is imperative for succesful MWR data processing!

####SPECIFY begin
#3 letter code for identification of measurement setup
meas_site=mosaic
# instrument name
mwr_name=hatpro
#times to process:
#set the following 1 line if you want to reprocess certain months (mmdd)
#months=(1301 1302 1303 1304 1305 1306 1307 1308 1309 1310 1311 1312)
months=(1909 1910 1911 1912 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010)
#months=(1101 1102 1103 1104 1105 1106 1107 1108 1109 1110 1111 1112)
#months=(1305)
#months=1101
#months=1307
#set the following 4 lines if you want to process only the current day (online-processing mode)
#tm=$(date +%y%m)
#ym=$(date --date 'yesterday' +%y%m)

#thisday=$(date +%y%m%d)
#yesterday=$(date --date 'yesterday' +%y%m%d)
#twodaysago=$(date --date '2 days ago' +%y%m%d)
#months=$ym
####SPECIFY end

cd /home/hatpro/mwr_pro_mosaic/mwr_pro/scripts

#copy parameter file to source directory where it will be executed
#cp par_mwr_pro_jue.pro par_mwr_pro.pro
cp par_mwr_pro_pol.pro par_mwr_pro.pro
cp par_mwr_pro.pro ../source/

#specify add mwr_pro directory to IDL path
IDL_PATH="<IDL_DEFAULT>"
path_idl="!path +': ../source'"

#get path parameters from path_info.txt
declare -a ARRAY

for LINE in `cat path_info_nds.txt`; do
 ARRAY[$count]=$LINE
 ((count++))
done

mwr_data_path=${ARRAY[0]}
idl_exec=${ARRAY[1]}
home_path=${ARRAY[2]}
www_path=${ARRAY[3]}/${meas_site}
mwr_dir=${mwr_data_path}/${meas_site}/${mwr_name}

# check whether month directories exist, if not create
for month in ${months[*]}; do
  year=20${month:0:2}
  monthcc=${month:2}
 if [ ! -d ${mwr_dir}/l0/${year} ]; then
  mkdir -p ${mwr_dir}/l0/${year}
 fi
 if [ ! -d ${mwr_dir}/l1/${year} ]; then
  mkdir -p ${mwr_dir}/l1/${year}
 fi
 if [ ! -d ${mwr_dir}/l2/${year} ]; then
  mkdir -p ${mwr_dir}/l2/${year}
 fi

 if [ ! -d ${mwr_dir}/l0/${year}/${monthcc}  ]; then
  mkdir -p ${mwr_dir}/l0/${year}/${monthcc}
 fi
 if [ ! -d ${mwr_dir}/l1/${year}/${monthcc}  ]; then
  mkdir -p ${mwr_dir}/l1/${year}/${monthcc}
 fi
 if [ ! -d ${mwr_dir}/l2/${year}/${monthcc}  ]; then
  mkdir -p ${mwr_dir}/l2/${year}/${monthcc}
 fi

 if [ ! -d ${www_path}/TB/${month}   ]; then
  mkdir -p ${www_path}/TB/${month}
 fi
 if [ ! -d ${www_path}/IWV_LWP/${month}   ]; then
  mkdir -p ${www_path}/IWV_LWP/${month}
 fi
 if [ ! -d ${www_path}/VAPs/radar_mwr/${month}   ]; then
  mkdir -p ${www_path}/VAPs/radar_mwr/${month}
 fi
if [ ! -d ${www_path}/VAPs/scans/${month}   ]; then
  mkdir -p ${www_path}/VAPs/scans/${month}
fi
 if [ ! -d ${www_path}/VAPs/Tq_profiles/${month}   ]; then
  mkdir -p ${www_path}/VAPs/Tq_profiles/${month}
fi
if [ ! -d ${www_path}/VAPs/cili/${month}   ]; then
  mkdir -p ${www_path}/VAPs/cili/${month}
 fi


# data_path=${mwr_dir}/data/raw/${month}
# plot_path0=${mwr_dir}/plots/level0/${month}
# plot_path1=${mwr_dir}/plots/level1/${month}
# plot_path2=${mwr_dir}/plots/level2/${month}

####SPECIFY begin
#set the following line if you would like to process only current day
#day_files=`echo $twodaysago $yesterday `
#set day_files to data path if you are reprocessing a whole month
#day_files=$(ls ${data_path})
#set day_files "by hand" to process only this day (e.g. testing purpose)
declare -a days=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
#declare -a days=("10", "11", "12")
#declare -a days=("11", "19")
#declare -a days=("12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
#declare -a days=("06" "07" "08" "09" "11" "12" "13")
#declare -a days=("17" "19" "22" "23" "24" "25" "28" "30" "31")
#declare -a days=("01" "05" "06" "07" "09" "15" "16" "18" "19" "20" "30")
#declare -a days=("06" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")

####SPECIFY end
# for day_file in $day_files; do
 for dayn in "${days[@]}"; do
  day_file=${month}${dayn} 
  day_file=`echo $day_file | cut -c -6`
  echo $day_file

# check whether day directories exist, if not create
  daycc=${day_file:4}
  if [ ! -d ${mwr_dir}/l0/${year}/${monthcc}/${daycc}  ]; then
  mkdir -p ${mwr_dir}/l0/${year}/${monthcc}/${daycc}
  fi
  if [ ! -d ${mwr_dir}/l1/${year}/${monthcc}/${daycc}  ]; then
  mkdir -p ${mwr_dir}/l1/${year}/${monthcc}/${daycc}
  fi
  if [ ! -d ${mwr_dir}/l2/${year}/${monthcc}/${daycc}  ]; then
  mkdir -p ${mwr_dir}/l2/${year}/${monthcc}/${daycc}
  fi

#copy parameter file to source directory where it will be executed
#cp par_mwr_pro_jue.pro par_mwr_pro.pro
  cp par_mwr_pro_pol.pro par_mwr_pro.pro

#call pl_mk.pro
  echo 'calling pl_mk_nds for: '$day_file
  pl_mk_call='pl_mk_nds, date='"'"${day_file}"'"',mwr_name='"'"${mwr_name}"'"',meas_site='"'"${meas_site}"'"',mwr_dir='"'"${mwr_dir}"'"',home_path='"'"${home_path}"'"
  echo "!path="${path_idl} > @l0_${meas_site}
  echo ${pl_mk_call} >> @l0_${meas_site}
  echo 'exit' >> @l0_${meas_site}
  ${idl_exec} @l0_${meas_site}
  rm @l0_${meas_site}

  daycc=${day_file:4}
  plot_path0=${mwr_dir}/l0/${year}/${monthcc}/${daycc}
  plot_path1=${mwr_dir}/l1/${year}/${monthcc}/${daycc}
  plot_path2=${mwr_dir}/l2/${year}/${monthcc}/${daycc}

#convert ps-files file to png-files:
#Please comment out lines, where no plots are produced
cd ${plot_path1} #NEW FIX PSTOIMG ISSUES
#1.) daily time series tbs
  ps_file=${plot_path1}/${day_file}_${meas_site}_l1.ps
  echo $ps_file
  if [ -f ${ps_file} ]; then
   png_file=${plot_path1}/${day_file}_${meas_site}_l1.png
   echo 'creating: ' ${png_file} 
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l1.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/TB/${month}
  fi
#2.) mean daily tb spectrum
  ps_file=${plot_path1}/${day_file}_${meas_site}_l1_sp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path1}/${day_file}_${meas_site}_l1_sp.png
   echo 'creating: ' ${png_file} 
   #pstoimg -density 200 ${ps_file}
   pstoimg -density 200 -type png ./${day_file}_${meas_site}_l1_sp.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/TB/${month}
 fi
#3.) daily flag time series
  ps_file=${plot_path1}/${day_file}_${meas_site}_flag.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path1}/${day_file}_${meas_site}_flag.png
   echo 'creating: ' ${png_file} 
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_flag.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/TB/${month}
 fi
cd ${plot_path2} #NEW FIX PSTOIMG ISSUES
#4.) daily time series lwp/iwv/wdl
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2a.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2a.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2a.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/IWV_LWP/${month}
  fi
#5.) daily time series of iwv/lwp azimuth-time contours
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2a_iwv_aztp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2a_iwv_aztp.png
   echo 'creating: ' ${png_file}
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2a_iwv_aztp.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/scans/${month}
  fi
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2a_lwp_aztp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2a_lwp_aztp.png
   echo 'creating: ' ${png_file}
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2a_lwp_aztp.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/scans/${month}
  fi
#6.) daily time series of IRT, ceilometer, LWP/IWV
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2_cili.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2_cili.png
   echo 'creating: ' ${png_file}
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2_cili.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/cili/${month}
  fi
#7.) daily time series 2b profiles (T and q along line of sight)
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2b_tze.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2b_tze.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2b_tze.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
  fi
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2b_hze.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2b_hze.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2b_hze.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
  fi
#8.) daily time series 2c profiles (T, Tpot, qiwv, rhlwp)
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2c_t.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2c_t.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2c_t.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2c_tpot.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2c_tpot.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2c_tpot.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2c_qiwv.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2c_qiwv.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2c_qiwv.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2c_rhlwp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2c_rhlwp.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2c_rhlwp.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi
#9.) daily time series radar measurements and 2a LWP              
  ps_file=${plot_path2}/${day_file}_${meas_site}_l2_radar.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2_radar.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2_radar.ps
   #NEW FIX PSTOIMG ISSUES
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/radar_mwr/${month}
  fi

  ps_file=${plot_path2}/${day_file}_${meas_site}_l2_radar_lwp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${meas_site}_l2_radar_lwp.png
   echo 'creating: ' ${png_file}       
   #pstoimg -flip r90 -density 200 ${ps_file}
   pstoimg -flip r90 -density 200 -type png ./${day_file}_${meas_site}_l2_radar_lwp.ps
   #NEW FIX PSTOIMG ISSUE
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/radar_mwr/${month}
  fi
cd /home/hatpro/mwr_pro_mosaic/mwr_pro/scripts
 done #end day files
done #end loop over months
