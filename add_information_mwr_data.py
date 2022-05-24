import numpy as np
import xarray as xr
import glob
import sys
import os
import datetime as dt
import pdb


"""
	Script to make changes to the HATPRO and MiRAC-P data products as suggested by reviewers.
	The changes merely address variable names and attributes.
	- import daily files
	- modify data set
	- export daily files
"""


# Variables that may not have fill values:
exclude_vars_fill_value = ['time', 'time_bnds', 'lat', 'lon', 'zsl', 'freq_sb', 'wl_irp']

valid_min = {'tb': np.array([2.7]).astype(np.float32)[0],
				'prw': np.array([0.0]).astype(np.float32)[0],
				'clwvi': np.array([-0.2]).astype(np.float32)[0],
				'ta': np.array([180.0]).astype(np.float32)[0],
				'hua': np.array([-0.5]).astype(np.float32)[0]}
valid_max = {'tb': np.array([330.0]).astype(np.float32)[0],
				'prw': np.array([100.0]).astype(np.float32)[0],
				'clwvi': np.array([3.0]).astype(np.float32)[0],
				'ta': np.array([330.0]).astype(np.float32)[0],
				'hua': np.array([30.0]).astype(np.float32)[0]}


def modify_nc(
	files,
	var,
	path_output_dir,
	files_HKD=None):

	"""
	Add some attributes (and variables) to existing netCDF files that have been
	processed by MWR_PRO (or NN_retrieval_miracp.py). Variable names will be 
	translated to CF1.6 conventions.

	Parameters:
	-----------
	files : list of str
		Potential netCDF files to modify. Length of the list should be 1.
	var : str
		Identifier of the main variable the file is about.
	path_output_dir : str
		Path where modified netCDF file is to be saved to.
	files_HKD : list of str
		Hourly housekeeping files of microwave radiometers from RPG.
	"""

	if var not in ['tb', 'prw', 'clwvi', 'hua', 'ta']:
		raise ValueError("'var' must be in ['tb', 'prw', 'clwvi', 'hua', 'ta'].")

	assert len(files) in [0,1]	# should be only one daily file each
	if len(files) == 1:
		files = files[0]

		# import: 
		DS = xr.open_dataset(files, decode_times=False)

		# modify and export datasets:
		# - add valid_min, valid_max as variable attributes to tb, ta, prw, clwvi, hua
		DS[var].attrs['valid_min'] = valid_min[var]
		DS[var].attrs['valid_max'] = valid_max[var]

		# - global attributes: time_start, time_end
		# DS.attrs['time_start'] = str(DS.time.values[0])
		# DS.attrs['time_end'] = str(DS.time.values[-1])
		DS.attrs['time_start'] = dt.datetime.utcfromtimestamp(DS.time.values[0]).strftime("%Y-%m-%d %H:%M:%SZ")
		DS.attrs['time_end'] = dt.datetime.utcfromtimestamp(DS.time.values[-1]).strftime("%Y-%m-%d %H:%M:%SZ")

		"""
		# add housekeeping data if available:
		if files_HKD:
			# concat hourly files
			DS_HKD = xr.open_mfdataset(files_HKD, concat_dim='time', combine='nested')

			# check if time axis is similar:
			pdb.set_trace()
			DS_HKD.close()
		"""

		# export:
		# adapt fill values: Make sure that _FillValue is not added to certain variables!
		if var in ['tb', 'prw', 'clwvi']:
			for kk in DS.variables:
				if kk in exclude_vars_fill_value:
					DS[kk].encoding["_FillValue"] = None
				elif kk != 'flag':
					DS[kk].encoding["_FillValue"] = float(-999.)
				else:
					DS[kk].encoding["_FillValue"] = np.array([0]).astype(np.short)[0]

		else:
			for kk in DS.variables:
				if kk in exclude_vars_fill_value:
					DS[kk].encoding["_FillValue"] = None
				elif kk == 'flag':
					DS[kk].encoding["_FillValue"] = np.array([0]).astype(np.short)[0]
				elif kk == f'{var}_offset':
					DS[kk].encoding["_FillValue"] = float(0.)
				else:
					DS[kk].encoding["_FillValue"] = float(-999.)

		# encode time:
		DS['time'].attrs['units'] = "seconds since 1970-01-01 00:00:00 UTC"
		DS['time'].encoding['units'] = 'seconds since 1970-01-01 00:00:00'
		DS['time'].encoding['dtype'] = 'double'

		# if "hatpro" in path_output_dir:
			# DS.to_netcdf(path_output_dir + os.path.basename(files), mode='w', format="NETCDF4_CLASSIC")
		# elif "mirac" in path_output_dir:
		DS.to_netcdf(path_output_dir + os.path.basename(files), mode='w', format="NETCDF4")
		DS.close()



###################################################################################################
###################################################################################################


# Paths:
path_output_hatpro = {'l1': "/data/obs/campaigns/mosaic/hatpro/l1/",
					'l2': "/data/obs/campaigns/mosaic/hatpro/l2/"}
path_output_mirac = {'l1': "/data/obs/campaigns/mosaic/mirac-p/l1/",
					'l2': "/data/obs/campaigns/mosaic/mirac-p/l2/"}
path_hatpro = {'l1': "/net/blanc/awalbroe/Data/MOSAiC_radiometers/backup_20220518/HATPRO_l1_v01/",
				'l2': "/net/blanc/awalbroe/Data/MOSAiC_radiometers/backup_20220518/HATPRO_l2_v01/"}
path_mirac = {'l1': "/net/blanc/awalbroe/Data/MOSAiC_radiometers/backup_20220518/MiRAC-P_l1_v01/",
				'l2': "/net/blanc/awalbroe/Data/MOSAiC_radiometers/backup_20220518/MiRAC-P_l2_v01/"}


# check existence of output paths:
for k in path_output_hatpro.keys():
	if not os.path.exists(path_output_hatpro[k]): os.makedirs(path_output_hatpro[k])
for k in path_output_mirac.keys():
	if not os.path.exists(path_output_mirac[k]): os.makedirs(path_output_mirac[k])


# date range:
date_start = dt.datetime(2019,9,20)
date_end = dt.datetime(2020,10,12)

# cycle through dates:
now_date = date_start
while now_date <= date_end:

	path_addition = f"{now_date.year:04}/{now_date.month:02}/{now_date.day:02}/"
	print(now_date)


	# HATPRO L1 zenith TBs:
	# check if paths exist:
	path_dir = path_hatpro['l1']
	path_output_dir = path_output_hatpro['l1'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwr00_l1_tb_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	# files_HKD = sorted(glob.glob(path_dir + "*.HKD.NC"))
	modify_nc(files, 'tb', path_output_dir)


	# HATPRO L1 BL TBs:
	# check if paths exist:
	path_dir = path_hatpro['l1']
	path_output_dir = path_output_hatpro['l1'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwrBL00_l1_tb_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	# files_HKD = sorted(glob.glob(path_dir + "*.HKD.NC"))
	modify_nc(files, 'tb', path_output_dir)


	# HATPRO L2 PRW:
	# check if paths exist:
	path_dir = path_hatpro['l2']
	path_output_dir = path_output_hatpro['l2'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwr00_l2_prw_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	modify_nc(files, 'prw', path_output_dir)


	# HATPRO L2 CLWVI:
	# check if paths exist:
	path_dir = path_hatpro['l2']
	path_output_dir = path_output_hatpro['l2'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwr00_l2_clwvi_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	modify_nc(files, 'clwvi', path_output_dir)


	# HATPRO L2 TA:
	# check if paths exist:
	path_dir = path_hatpro['l2']
	path_output_dir = path_output_hatpro['l2'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwr00_l2_ta_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	modify_nc(files, 'ta', path_output_dir)


	# HATPRO L2 TA_BL:
	# check if paths exist:
	path_dir = path_hatpro['l2']
	path_output_dir = path_output_hatpro['l2'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwrBL00_l2_ta_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	modify_nc(files, 'ta', path_output_dir)


	# HATPRO L2 HUA:
	# check if paths exist:
	path_dir = path_hatpro['l2']
	path_output_dir = path_output_hatpro['l2'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"ioppol_tro_mwr00_l2_hua_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	modify_nc(files, 'hua', path_output_dir)


	# MIRAC-P L1 TBs:
	# check if paths exist:
	path_dir = path_mirac['l1']
	path_output_dir = path_output_mirac['l1'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"MOSAiC_uoc_lhumpro-243-340_l1_tb_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	# files_HKD = sorted(glob.glob(path_dir + "*.HKD.NC"))
	modify_nc(files, 'tb', path_output_dir)


	# MIRAC-P L2 PRW:
	# check if paths exist:
	path_dir = path_mirac['l2']
	path_output_dir = path_output_mirac['l2'] + path_addition
	if not os.path.exists(path_output_dir):	os.makedirs(path_output_dir)

	# identify correct files:
	files = sorted(glob.glob(path_dir + f"MOSAiC_uoc_lhumpro-243-340_l2_prw_v01_{now_date.year:04}{now_date.month:02}{now_date.day:02}*.nc"))
	modify_nc(files, 'prw', path_output_dir)


	now_date += dt.timedelta(days=1)
