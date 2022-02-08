import numpy as np
import netCDF4 as nc
import datetime as dt
# import pandas as pd
# import xarray as xr
import copy
import pdb
import os
import glob
import sys
import warnings
import csv
from met_tools import *
from data_tools import *



def import_mirac_BRT(
	filename,
	keys='basic'):

	"""
	Importing automatically created MiRAC-P BRT hourly files
	with the ending .BRT.NC in the level 1 folder. Time will be
	converted to seconds since 1970-01-01 00:00:00 UTC.

	Parameters:
	-----------
	filename : str
		Path and filename of MiRAC-P .BRT.NC data.
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (keys='all') or import basic keys
		that the author considers most important (keys='basic')
		or leave this argument out.
	"""

	# BRT in K!

	file_nc = nc.Dataset(filename)

	if keys == 'basic':
		keys = ['time', 'RF', 'TBs', 'Freq']

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	reftime = dt.datetime(1970,1,1,0,0,0)
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in MiRAC-P .BRT.NC file." % key)

		mwr_dict[key] = np.asarray(file_nc.variables[key])

		if key == 'time':	# convert to sec since 1970-01-01 00:00:00 UTC (USE FLOAT64)
			mwr_dict['time'] = (np.float64(datetime_to_epochtime(dt.datetime(2001,1,1,0,0,0))) +
								mwr_dict[key].astype(np.float64))

	return mwr_dict


def import_mirac_BRT_RPG_daterange(
	path_data,
	date_start,
	date_end,
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the BRT time
	series of each day so that you'll have one dictionary that will contain the TBs
	for the entire date range period.

	Parameters:
	-----------
	path_data : str
		Base path of MiRAC-P level 1 data. This directory contains subfolders representing the 
		year, which, in turn, contain months, which contain day subfolders. Example:
		path_data = "/data/obs/campaigns/mosaic/mirac-p/l1/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")

	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'RF']	# keys with time as coordinate
	mwr_time_freq_keys = ['TBs']	# keys with time and frequency as coordinates
	mwr_freq_keys = ['Freq']		# keys with frequency as coordinate

	# mwr_master_dict (output) will contain all desired variables on time axis for entire date range:
	mwr_master_dict = dict()
	n_seconds = n_days*86400		# max number of seconds: n_days*86400
	n_freq = 8						# number of frequencies (inquired from .BRT.NC file)
	for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_seconds,), np.nan)
	for mtfk in mwr_time_freq_keys: mwr_master_dict[mtfk] = np.full((n_seconds, n_freq), np.nan)
	for mfk in mwr_freq_keys: mwr_master_dict[mfk] = np.full((n_freq,), np.nan)


	# Load the TBs into mwr_master_dict:
	# cycle through all years, all months and days:
	time_index = 0	# this index will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	day_index = 0	# will increase for each day
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on RPG TBs, MiRAC-P BRT, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# list of .BRT.NC files: Sorting is important as this will
		# ensure automatically that the time series of each hour will
		# be concatenated appropriately!
		mirac_nc = sorted(glob.glob(day_path + "*.BRT.NC"))
		if len(mirac_nc) == 0:
			if verbose >= 2:
				warnings.warn("No .BRT.NC files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for l1_file in mirac_nc: 
			mwr_dict = import_mirac_BRT(l1_file)

			n_time = len(mwr_dict['time'])
			time_shape = mwr_dict['time'].shape

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				if mwr_key in mwr_time_keys:
					mwr_master_dict[mwr_key][time_index:time_index + n_time] = mwr_dict[mwr_key]

				elif mwr_key in mwr_time_freq_keys:
					mwr_master_dict[mwr_key][time_index:time_index + n_time,:] = mwr_dict[mwr_key]

				elif mwr_dict[mwr_key].shape == ():
					mwr_master_dict[mwr_key][day_index] = mwr_dict[mwr_key]

				elif mwr_key in mwr_freq_keys:	# frequency will be handled after the for loop
					continue

				else:
						raise ValueError("The length of one used variable ('%s') of MiRAC-P .BRT.NC data "%(mwr_key) +
							"neither equals the length of the time axis nor equals 1.")

			time_index = time_index + n_time
		day_index = day_index + 1


	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))

	else:
		# Save frequency array into mwr_master_dict:
		for fkey in mwr_freq_keys: mwr_master_dict[fkey] = mwr_dict[fkey]

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		time_freq_shape_old = mwr_master_dict[mwr_time_freq_keys[0]].shape
		for mwr_key in mwr_master_dict.keys():
			shape_new = mwr_master_dict[mwr_key].shape
			if shape_new == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]
			elif shape_new == time_freq_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1,:]

	return mwr_master_dict


def import_mirac_level2a(
	filename,
	keys='basic',
	minute_avg=False):

	"""
	Importing MiRAC-P level 2a (integrated quantities, e.g. IWV, LWP).

	Parameters:
	-----------
	filename : str
		Path and filename of mwr data (level2a).
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (keys='all') or import basic keys
		that the author considers most important (keys='basic')
		or leave this argument out.
	minute_avg : bool
		If True: averages over one minute are computed and returned instead of False when all
		data points are returned (more outliers, higher memory usage).
	"""

	file_nc = nc.Dataset(filename)

	if keys == 'basic': 
		keys = ['time', 'lat', 'lon', 'zsl', 'azi', 'ele', 'flag']
		if 'clwvi_' in filename:
			for add_key in ['clwvi']: keys.append(add_key)
		if 'prw_' in filename:
			for add_key in ['prw']: keys.append(add_key)

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in level 2a file." % key)
		mwr_dict[key] = np.asarray(file_nc.variables[key])


	if 'time' in keys:	# avoid nasty digita after decimal point
		mwr_dict['time'] = np.rint(mwr_dict['time']).astype(float)
		time_shape_old = mwr_dict['time'].shape

		if minute_avg:
			# start the timer at the first time, when seconds is 00 (e.g. 09:31:00):
			time0 = mwr_dict['time'][0]		# start time in sec since 1970-01-01...
			dt_time0 = dt.datetime.utcfromtimestamp(mwr_dict['time'][0])
			dt_time0_Y = dt_time0.year
			dt_time0_M = dt_time0.month
			dt_time0_D = dt_time0.day
			dt_time0_s = dt_time0.second
			dt_time0_m = dt_time0.minute
			dt_time0_h = dt_time0.hour
			if dt_time0_s != 0:		# then the array mwr_dict['time'] does not start at second 0
				start_time = datetime_to_epochtime(dt.datetime(dt_time0_Y, dt_time0_M,
													dt_time0_D, dt_time0_h, dt_time0_m+1, 0))
			else:
				start_time = time0

			if np.abs(start_time - time0) >= 60:
				print("Start time is far off the first time point in this file.")
				pdb.set_trace()
			# compute minute average
			n_minutes = int(np.ceil((mwr_dict['time'][-1] - start_time)/60))	# number of minutes
			min_time_idx_save = 0		# saves the last min_time_index value to speed up computation
			for min_count in range(n_minutes):
				# find time_idx when time is in the correct minute:
				# slower version:
				# # # # min_time_idx = np.argwhere((mwr_dict['time'] >= (start_time + min_count*60)) & 
								# # # # (mwr_dict['time'] < (start_time + (min_count+1)*60))).flatten()
				# faster version:
				min_time_idx = np.argwhere((mwr_dict['time'][min_time_idx_save:min_time_idx_save+180] >= (start_time + min_count*60)) & 
								(mwr_dict['time'][min_time_idx_save:min_time_idx_save+180] < (start_time + (min_count+1)*60))).flatten()

				# it may occur that no measurement exists in a certain minute-range. Then
				# we cannot compute the average but simply set that minute to nan.
				if len(min_time_idx) == 0:
					for key in keys:
						if key == 'time':
							mwr_dict['time'][min_count] = start_time + min_count*60
						elif mwr_dict[key].shape == time_shape_old and key != 'flag':
							mwr_dict[key][min_count] = np.nan
						elif mwr_dict[key].shape == time_shape_old and key == 'flag':
							mwr_dict[key][min_count] = 99		# np.nan not possible because int is required
				else:
					min_time_idx = min_time_idx + min_time_idx_save		# also belonging to the 'faster version'
					min_time_idx_save = min_time_idx[-1]				# also belonging to the 'faster version'
					for key in keys:
						if key == 'time':
							mwr_dict['time'][min_count] = start_time + min_count*60
						elif mwr_dict[key].shape == time_shape_old and key != 'flag':
							if min_time_idx[-1] < len(mwr_dict['time']):
								mwr_dict[key][min_count] = np.nanmean(mwr_dict[key][min_time_idx])
							else:
								mwr_dict[key][min_count] = np.nanmean(mwr_dict[key][min_time_idx[0]:])
						elif mwr_dict[key].shape == time_shape_old and key == 'flag':
							# find out how many entries show flag > 0. Then if it exceeds a threshold
							# the whole minute is flagged. If threshold not exceeded, minute is not
							# flagged!
							if min_time_idx[-1] < len(mwr_dict['time']):
								if np.count_nonzero(mwr_dict[key][min_time_idx]) > len(min_time_idx)/10:	
									# then there are too many flags set... so flag the whole minute:
									mwr_dict[key][min_count] = 99
								else:
									mwr_dict[key][min_count] = 0
							else:
								if np.count_nonzero(mwr_dict[key][min_time_idx[0]:]) > len(min_time_idx)/10:
									# then there are too many flags set... so flag the whole minute:
									mwr_dict[key][min_count] = 99
								else:
									mwr_dict[key][min_count] = 0

			# truncate time arrays to reduce memory usage!
			for key in keys:
				if mwr_dict[key].shape == time_shape_old:
					mwr_dict[key] = mwr_dict[key][:n_minutes]

	else:
		if minute_avg:
			raise KeyError("'time' must be included in the list of keys that will be imported for minute averages.")

	return mwr_dict


def import_mirac_level1b_daterange(
	path_data,
	date_start,
	date_end,
	vers='v01',
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the level 1b TB time
	series of each day so that you'll have one dictionary, whose 'TB' will contain the IWV
	for the entire date range period.

	Parameters:
	-----------
	path_data : str
		Base path of level 1 data. This directory contains subfolders representing the year, which,
		in turn, contain months, which contain day subfolders. Example path_data:
		"/data/obs/campaigns/mosaic/mirac-p/l1/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	vers : str
		Indicates the mwr_pro output version number. Valid options: 'v00', and 'v01'. Default: 'v01'
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	if vers not in ['v00', 'v01']:
		raise ValueError("In import_hatpro_level1b_daterange, the argument 'vers' must be one of the" +
							" following options: 'v00', 'v01'")

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")


	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1
	n_freq = 8			# inquired from level 1 data, number of frequencies

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'flag', 'ta', 'pa', 'hur']				# keys with time as coordinate
	mwr_freq_keys = ['freq_sb', 'freq_shift', 'tb_absolute_accuracy']	# keys with frequency as coordinate
	mwr_time_freq_keys = ['tb', 'tb_bias_estimate']						# keys with frequency and time as coordinates

	# mwr_master_dict (output) will contain all desired variables on specific axes:
	mwr_master_dict = dict()

	# max number of seconds: n_days*86400
	n_seconds = n_days*86400
	for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_seconds,), np.nan)
	for mtkab in mwr_time_freq_keys: mwr_master_dict[mtkab] = np.full((n_seconds, n_freq), np.nan)
	for mtkab in mwr_freq_keys: mwr_master_dict[mtkab] = np.full((n_freq,), np.nan)
	mwr_master_dict['tb_cov'] = np.full((n_freq,n_freq), np.nan)		# has got a special shape -> handled manually


	# cycle through all years, all months and days:
	time_index = 0	# this index (for lvl 1) will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	day_index = 0	# same as above, but only increases by 1 for each day
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on MiRAC-P Level 1, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# list of v01 tb files with zenith scan:
		file_nc = sorted(glob.glob(day_path + "*_mwr01_l1_tb_%s_*.nc"%vers))

		if len(file_nc) == 0:
			if verbose >= 2:
				warnings.warn("No netcdf files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for lvl_nc in file_nc: 
			mwr_dict = import_hatpro_level1b(lvl_nc, keys='all')

			n_time = len(mwr_dict['time'])
			cur_time_shape = mwr_dict['time'].shape

			# save to mwr_master_dict
			for mwr_key in mwr_time_keys:
				mwr_master_dict[mwr_key][time_index:time_index + n_time] = mwr_dict[mwr_key]
			for mwr_key in mwr_time_freq_keys:
				mwr_master_dict[mwr_key][time_index:time_index + n_time,:] = mwr_dict[mwr_key]


		time_index = time_index + n_time
		day_index = day_index + 1

	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))
	else:
		# assign frequency to master dict:
		for mwr_key in mwr_freq_keys: mwr_master_dict[mwr_key] = mwr_dict[mwr_key]
		mwr_master_dict['tb_cov'] = mwr_dict['tb_cov']

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		time_freq_shape_old = mwr_master_dict['tb'].shape
		for mwr_key in mwr_master_dict.keys():
			shape_new = mwr_master_dict[mwr_key].shape
			if shape_new == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]
			elif shape_new == time_freq_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1, :]

	return mwr_master_dict


def import_mirac_level2a_daterange(
	path_data,
	date_start,
	date_end,
	which_retrieval='both',
	vers='v01',
	minute_avg=False,
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the level 2a data time
	series of each day so that you'll have one dictionary, whose e.g. 'IWV' will contain the IWV
	for the entire date range period.

	Parameters:
	-----------
	path_data : str
		Base path of level 2a data. This directory contains subfolders representing the year, which,
		in turn, contain months, which contain day subfolders. Example path_data:
		"/data/obs/campaigns/mosaic/mirac-p/l2/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	which_retrieval : str, optional
		This describes which variable(s) will be loaded. Options: 'iwv' or 'prw' will load the
		integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. 'both' will 
		load both. Default: 'both'
	vers : str
		Indicates the mwr_pro output version number. Valid options: 'v01'
	minute_avg : bool
		If True: averages over one minute are computed and returned instead of False when all
		data points are returned (more outliers, higher memory usage).
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	# check if the input of the retrieval variable is okay:
	if not isinstance(which_retrieval, str):
			raise TypeError("Argument 'which_retrieval' must be a string. Options: 'iwv' or 'prw' will load the " +
				"integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. " +
				"'both' will load both. Default: 'both'")

	else:
		if which_retrieval not in ['prw', 'iwv', 'clwvi', 'lwp', 'both']:
			raise ValueError("Argument 'which_retrieval' must be one of the following options: 'iwv' or 'prw' will load the " +
				"integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. " +
				"'both' will load both. Default: 'both'")

		else:
				if which_retrieval == 'iwv':
					which_retrieval = ['prw']
				elif which_retrieval == 'lwp':
					which_retrieval = ['clwvi']
				elif which_retrieval == 'both':
					which_retrieval = ['prw', 'clwvi']
				else:
					raise ValueError("Argument '" + which_retrieval + "' not recognised. Please use one of the following options: " +
						"'iwv' or 'prw' will load the " +
						"integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. " +
						"'both' will load both. Default: 'both'")
					

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")


	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1
	n_ret = 86			# inquired from level 2a data, number of available elevation angles in retrieval

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'azi', 'ele', 'flag', 'lat', 'lon', 'zsl']				# keys with time as coordinate

	# mwr_master_dict (output) will contain all desired variables on specific axes:
	# e.g. level 2a and 2b have got the same time axis (according to pl_mk_nds.pro)
	# and azimuth and elevation angles.
	mwr_master_dict = dict()
	if minute_avg:	# max number of minutes: n_days*1440
		n_minutes = n_days*1440
		for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_minutes,), np.nan)

		if 'prw' in which_retrieval:
			mwr_master_dict['prw'] = np.full((n_minutes,), np.nan)

		if 'clwvi' in which_retrieval:
			mwr_master_dict['clwvi'] = np.full((n_minutes,), np.nan)

	else:			# max number of seconds: n_days*86400
		n_seconds = n_days*86400
		for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_seconds,), np.nan)

		if 'prw' in which_retrieval:
			mwr_master_dict['prw'] = np.full((n_seconds,), np.nan)

		if 'clwvi' in which_retrieval:
			mwr_master_dict['clwvi'] = np.full((n_seconds,), np.nan)


	# cycle through all years, all months and days:
	time_index = 0	# this index (for lvl 2a) will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	day_index = 0	# same as above, but only increases by 1 for each day
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on MiRAC-P Level 2a, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# list of files:
		mirac_level2_nc = sorted(glob.glob(day_path + "*.nc"))
		# filter for i01 files:
		mirac_level2_nc = [lvl2_nc for lvl2_nc in mirac_level2_nc if vers in lvl2_nc]

		if len(mirac_level2_nc) == 0:
			if verbose >= 2:
				warnings.warn("No netcdf files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# identify level 2a files:
		mirac_level2a_nc = []
		for lvl2_nc in mirac_level2_nc:
			for wr in which_retrieval:
				if wr + '_' in lvl2_nc:
					mirac_level2a_nc.append(lvl2_nc)

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for lvl2_nc in mirac_level2a_nc: 
			mwr_dict = import_mirac_level2a(lvl2_nc, minute_avg=minute_avg)

			n_time = len(mwr_dict['time'])
			cur_time_shape = mwr_dict['time'].shape

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				mwr_key_shape = mwr_dict[mwr_key].shape
				if mwr_key_shape == cur_time_shape:	# then the variable is on time axis:
					mwr_master_dict[mwr_key][time_index:time_index + n_time] = mwr_dict[mwr_key]

				elif len(mwr_dict[mwr_key]) == 1:
					mwr_master_dict[mwr_key][day_index:day_index + 1] = mwr_dict[mwr_key]

				else:
					raise RuntimeError("Something went wrong in the " +
						"import_mirac_level2a_daterange routine. Unexpected MWR variable dimension. " + 
						"The length of one used variable ('%s') of level 2a data "%(mwr_key) +
							"neither equals the length of the time axis nor equals 1.")


		time_index = time_index + n_time
		day_index = day_index + 1

	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))
	else:

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		for mwr_key in mwr_master_dict.keys():
			if mwr_master_dict[mwr_key].shape == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]

	return mwr_master_dict


def import_hatpro_level1b(
	filename,
	keys='basic'):

	"""
	Importing HATPRO level 1b (zenith TBs in K). Can also be used to import
	MiRAC-P level 1b (zenith TBs in K) data if it was processed with mwr_pro.

	Parameters:
	-----------
	filename : str
		Path and filename of mwr data (level2b).
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (specify keys='all') or import basic keys
		that the author considers most important (specify keys='basic')
		or leave this argument out.).
	"""

	file_nc = nc.Dataset(filename)

	if keys == 'basic': 
		keys = ['time', 'freq_sb', 'flag', 'tb']

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key '%s'. Key not found in level 1 file." % key)
		mwr_dict[key] = np.asarray(file_nc.variables[key])


	if 'time' in keys:	# avoid nasty digita after decimal point
		mwr_dict['time'] = np.rint(mwr_dict['time']).astype(float)

	return mwr_dict


def import_hatpro_level2a(
	filename,
	keys='basic',
	minute_avg=False):

	"""
	Importing HATPRO level 2a (integrated quantities, e.g. IWV, LWP).

	Parameters:
	-----------
	filename : str
		Path and filename of mwr data (level2a).
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (keys='all') or import basic keys
		that the author considers most important (keys='basic')
		or leave this argument out.
	minute_avg : bool
		If True: averages over one minute are computed and returned instead of False when all
		data points are returned (more outliers, higher memory usage).
	"""

	file_nc = nc.Dataset(filename)

	if keys == 'basic': 
		keys = ['time', 'lat', 'lon', 'zsl', 'flag']
		if 'clwvi_' in filename:
			for add_key in ['clwvi', 'clwvi_err', 'clwvi_offset']: keys.append(add_key)
		if 'prw_' in filename:
			for add_key in ['prw', 'prw_err', 'prw_offset']: keys.append(add_key)

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in level 2a file." % key)
		mwr_dict[key] = np.asarray(file_nc.variables[key])


	if 'time' in keys:	# avoid nasty digita after decimal point
		mwr_dict['time'] = np.rint(mwr_dict['time']).astype(float)
		time_shape_old = mwr_dict['time'].shape

		if minute_avg:
			# start the timer at the first time, when seconds is 00 (e.g. 09:31:00):
			time0 = mwr_dict['time'][0]		# start time in sec since 1970-01-01...
			dt_time0 = dt.datetime.utcfromtimestamp(mwr_dict['time'][0])
			dt_time0_Y = dt_time0.year
			dt_time0_M = dt_time0.month
			dt_time0_D = dt_time0.day
			dt_time0_s = dt_time0.second
			dt_time0_m = dt_time0.minute
			dt_time0_h = dt_time0.hour
			if dt_time0_s != 0:		# then the array mwr_dict['time'] does not start at second 0
				start_time = datetime_to_epochtime(dt.datetime(dt_time0_Y, dt_time0_M, dt_time0_D,
													dt_time0_h, dt_time0_m+1, 0))
			else:
				start_time = time0

			if np.abs(start_time - time0) >= 60:
				print("Start time is far off the first time point in this file.")
				pdb.set_trace()
			# compute minute average
			n_minutes = int(np.ceil((mwr_dict['time'][-1] - start_time)/60))	# number of minutes
			min_time_idx_save = 0		# saves the last min_time_index value to speed up computation
			for min_count in range(n_minutes):
				# find time_idx when time is in the correct minute:
				# slower version:
				# # # # min_time_idx = np.argwhere((mwr_dict['time'] >= (start_time + min_count*60)) & 
								# # # # (mwr_dict['time'] < (start_time + (min_count+1)*60))).flatten()
				# faster version:
				min_time_idx = np.argwhere((mwr_dict['time'][min_time_idx_save:min_time_idx_save+180] >= (start_time + min_count*60)) & 
								(mwr_dict['time'][min_time_idx_save:min_time_idx_save+180] < (start_time + (min_count+1)*60))).flatten()

				# it may occur that no measurement exists in a certain minute-range. Then
				# we cannot compute the average but simply set that minute to nan.
				if len(min_time_idx) == 0:
					for key in keys:
						if key == 'time':
							mwr_dict['time'][min_count] = start_time + min_count*60
						elif mwr_dict[key].shape == time_shape_old and key != 'flag':
							mwr_dict[key][min_count] = np.nan
						elif mwr_dict[key].shape == time_shape_old and key == 'flag':
							mwr_dict[key][min_count] = 99		# np.nan not possible because int is required
				else:
					min_time_idx = min_time_idx + min_time_idx_save		# also belonging to the 'faster version'
					min_time_idx_save = min_time_idx[-1]				# also belonging to the 'faster version'
					for key in keys:
						if key == 'time':
							mwr_dict['time'][min_count] = start_time + min_count*60
						elif mwr_dict[key].shape == time_shape_old and key != 'flag':
							if min_time_idx[-1] < len(mwr_dict['time']):
								mwr_dict[key][min_count] = np.nanmean(mwr_dict[key][min_time_idx])
							else:
								mwr_dict[key][min_count] = np.nanmean(mwr_dict[key][min_time_idx[0]:])
						elif mwr_dict[key].shape == time_shape_old and key == 'flag':
							# find out how many entries show flag > 0. Then if it exceeds a threshold
							# the whole minute is flagged. If threshold not exceeded, minute is not
							# flagged!
							if min_time_idx[-1] < len(mwr_dict['time']):
								if np.count_nonzero(mwr_dict[key][min_time_idx]) > len(min_time_idx)/10:	
									# then there are too many flags set... so flag the whole minute:
									mwr_dict[key][min_count] = 99
								else:
									mwr_dict[key][min_count] = 0
							else:
								if np.count_nonzero(mwr_dict[key][min_time_idx[0]:]) > len(min_time_idx)/10:
									# then there are too many flags set... so flag the whole minute:
									mwr_dict[key][min_count] = 99
								else:
									mwr_dict[key][min_count] = 0

			# truncate time arrays to reduce memory usage!
			for key in keys:
				if mwr_dict[key].shape == time_shape_old:
					mwr_dict[key] = mwr_dict[key][:n_minutes]

	else:
		if minute_avg:
			raise KeyError("'time' must be included in the list of keys that will be imported for minute averages.")

	return mwr_dict


def import_hatpro_level2b(
	filename,
	keys='basic',
	minute_avg=False):

	"""
	Importing HATPRO level 2b (zenith profiles, temperature or humidity 
	(in K or kg m^-3, respectively).

	Parameters:
	-----------
	filename : str
		Path and filename of mwr data (level2b).
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (specify keys='all') or import basic keys
		that the author considers most important (specify keys='basic')
		or leave this argument out.
	minute_avg : bool
		If True: averages over one minute are computed and returned. False: all
		data points are returned (more outliers, higher memory usage but may result in
		long computation time).
	"""

	file_nc = nc.Dataset(filename)

	if keys == 'basic': 
		keys = ['time', 'lat', 'lon', 'zsl', 'height', 'flag']
		if 'hua_' in filename:
			for add_key in ['hua']: keys.append(add_key)
		if 'ta_' in filename:
			for add_key in ['ta']: keys.append(add_key)

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key '%s'. Key not found in level 2b file." % key)
		mwr_dict[key] = np.asarray(file_nc.variables[key])


	if 'time' in keys:	# avoid nasty digita after decimal point
		mwr_dict['time'] = np.rint(mwr_dict['time']).astype(float)
		time_shape_old = mwr_dict['time'].shape

		if minute_avg:
			# start the timer at the first time, when seconds is 00 (e.g. 09:31:00):
			time0 = mwr_dict['time'][0]		# start time in sec since 1970-01-01...
			dt_time0 = dt.datetime.utcfromtimestamp(mwr_dict['time'][0])
			dt_time0_Y = dt_time0.year
			dt_time0_M = dt_time0.month
			dt_time0_D = dt_time0.day
			dt_time0_s = dt_time0.second
			dt_time0_m = dt_time0.minute
			dt_time0_h = dt_time0.hour
			if dt_time0_s != 0:		# then the array mwr_dict['time'] does not start at second 0
				start_time = datetime_to_epochtime(dt.datetime(dt_time0_Y, dt_time0_M,
								dt_time0_D, dt_time0_h, dt_time0_m+1, 0))
			else:
				start_time = time0

			if np.abs(start_time - time0) >= 60:
				print("Start time is far off the first time point in this file.")
				pdb.set_trace()
			# compute minute average
			n_minutes = int(np.ceil((mwr_dict['time'][-1] - start_time)/60))	# number of minutes
			min_time_idx_save = 0		# saves the last min_time_index value to speed up computation
			for min_count in range(n_minutes):
				# find time_idx when time is in the correct minute:
				# slower version:
				# # # # min_time_idx = np.argwhere((mwr_dict['time'] >= (start_time + min_count*60)) & 
								# # # # (mwr_dict['time'] < (start_time + (min_count+1)*60))).flatten()
				# faster version:
				min_time_idx = np.argwhere((mwr_dict['time'][min_time_idx_save:min_time_idx_save+180] >= (start_time + min_count*60)) & 
								(mwr_dict['time'][min_time_idx_save:min_time_idx_save+180] < (start_time + (min_count+1)*60))).flatten()

				# it may occur that no measurement exists in a certain minute-range. Then
				# we cannot compute the average but simply set that minute to nan.
				if len(min_time_idx) == 0:
					for key in keys:
						if key == 'time':
							mwr_dict['time'][min_count] = start_time + min_count*60
						elif mwr_dict[key].shape == time_shape_old and key != 'flag':
							mwr_dict[key][min_count] = np.nan
						elif mwr_dict[key].shape == time_shape_old and key == 'flag':
							mwr_dict[key][min_count] = 99		# np.nan not possible because int is required
				else:
					min_time_idx = min_time_idx + min_time_idx_save		# also belonging to the 'faster version'
					min_time_idx_save = min_time_idx[-1]				# also belonging to the 'faster version'
					for key in keys:
						if key == 'time':
							mwr_dict['time'][min_count] = start_time + min_count*60
						elif mwr_dict[key].shape == time_shape_old and key != 'flag':
							if min_time_idx[-1] < len(mwr_dict['time']):
								mwr_dict[key][min_count] = np.nanmean(mwr_dict[key][min_time_idx])
							else:
								mwr_dict[key][min_count] = np.nanmean(mwr_dict[key][min_time_idx[0]:])
						elif mwr_dict[key].shape == time_shape_old and key == 'flag':
							# find out how many entries show flag > 0. Then if it exceeds a threshold
							# the whole minute is flagged. If threshold not exceeded, minute is not
							# flagged!
							if min_time_idx[-1] < len(mwr_dict['time']):
								if np.count_nonzero(mwr_dict[key][min_time_idx]) > len(min_time_idx)/10:	
									# then there are too many flags set... so flag the whole minute:
									mwr_dict[key][min_count] = 99
								else:
									mwr_dict[key][min_count] = 0
							else:
								if np.count_nonzero(mwr_dict[key][min_time_idx[0]:]) > len(min_time_idx)/10:
									# then there are too many flags set... so flag the whole minute:
									mwr_dict[key][min_count] = 99
								else:
									mwr_dict[key][min_count] = 0

			# truncate time arrays to reduce memory usage!
			for key in keys:
				if mwr_dict[key].shape == time_shape_old:
					mwr_dict[key] = mwr_dict[key][:n_minutes]

	else:
		if minute_avg:
			raise KeyError("'time' must be included in the list of keys that will be imported for minute averages.")


	return mwr_dict


def import_hatpro_level2c(
	filename,
	keys='basic'):

	"""
	Importing HATPRO level 2c (boundary layer profiles, temperature (or humidity)
	(in K or kg m^-3, respectively).

	Parameters:
	-----------
	filename : str
		Path and filename of mwr data (level2c).
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (specify keys='all') or import basic keys
		that the author considers most important (specify keys='basic')
		or leave this argument out.
	"""

	file_nc = nc.Dataset(filename)

	if keys == 'basic': 
		keys = ['time', 'lat', 'lon', 'zsl', 'height', 'flag']
		if 'hua_' in filename:
			for add_key in ['hua']: keys.append(add_key)
		if 'ta_' in filename:
			for add_key in ['ta']: keys.append(add_key)

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key '%s'. Key not found in level 2c file." % key)
		mwr_dict[key] = np.asarray(file_nc.variables[key])


	if 'time' in keys:	# avoid nasty digita after decimal point
		mwr_dict['time'] = np.rint(mwr_dict['time']).astype(float)
		time_shape_old = mwr_dict['time'].shape

	return mwr_dict


def import_hatpro_level1b_daterange(
	path_data,
	date_start,
	date_end,
	vers='v01',
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the level 1b TB time
	series of each day so that you'll have one dictionary, whose 'TB' will contain the IWV
	for the entire date range period.

	Parameters:
	-----------
	path_data : str
		Base path of level 1 data. This directory contains subfolders representing the year, which,
		in turn, contain months, which contain day subfolders. Example path_data:
		"/data/obs/campaigns/mosaic/hatpro/l1/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	vers : str
		Indicates the mwr_pro output version number. Valid options: 'i01', 'v00', and 'v01'. Default: 'v01'
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	if vers not in ['i01', 'v00', 'v01']:
		raise ValueError("In import_hatpro_level1b_daterange, the argument 'vers' must be one of the" +
							" following options: 'i01', 'v00', 'v01'")

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")


	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1
	n_freq = 14			# inquired from level 1 data, number of frequencies

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'flag']				# keys with time as coordinate
	mwr_freq_keys = ['freq_sb']						# keys with frequency as coordinate
	mwr_time_freq_keys = ['tb']						# keys with frequency and time as coordinates

	# mwr_master_dict (output) will contain all desired variables on specific axes:
	mwr_master_dict = dict()

	# max number of seconds: n_days*86400
	n_seconds = n_days*86400
	for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_seconds,), np.nan)
	for mtkab in mwr_time_freq_keys: mwr_master_dict[mtkab] = np.full((n_seconds, n_freq), np.nan)
	for mtkab in mwr_freq_keys: mwr_master_dict[mtkab] = np.full((n_freq,), np.nan)


	# cycle through all years, all months and days:
	time_index = 0	# this index (for lvl 1) will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	day_index = 0	# same as above, but only increases by 1 for each day
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on HATPRO Level 1, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# list of v01 tb files with zenith scan:
		hatpro_nc = sorted(glob.glob(day_path + "*_mwr00_*_%s_*.nc"%vers))

		if len(hatpro_nc) == 0:
			if verbose >= 2:
				warnings.warn("No netcdf files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for lvl_nc in hatpro_nc: 
			mwr_dict = import_hatpro_level1b(lvl_nc)

			n_time = len(mwr_dict['time'])
			cur_time_shape = mwr_dict['time'].shape

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				mwr_key_shape = mwr_dict[mwr_key].shape
				if mwr_key_shape == cur_time_shape:	# then the variable is on time axis:
					mwr_master_dict[mwr_key][time_index:time_index + n_time] = mwr_dict[mwr_key]

				elif mwr_key_shape == (n_time, n_freq):
					mwr_master_dict[mwr_key][time_index:time_index + n_time,:] = mwr_dict[mwr_key]

				elif mwr_key in mwr_freq_keys:	# will be handled after the loop
					continue

				else:
					raise RuntimeError("Something went wrong in the " +
						"import_hatpro_level1b_daterange routine. Unexpected MWR variable dimension of '%s'. "%(mwr_key))


		time_index = time_index + n_time
		day_index = day_index + 1

	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))
	else:
		# assign frequency to master dict:
		for mwr_key in mwr_freq_keys: mwr_master_dict[mwr_key] = mwr_dict[mwr_key]

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		time_freq_shape_old = mwr_master_dict['tb'].shape
		for mwr_key in mwr_master_dict.keys():
			shape_new = mwr_master_dict[mwr_key].shape
			if shape_new == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]
			elif shape_new == time_freq_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1, :]

	return mwr_master_dict


def import_hatpro_level2a_daterange(
	path_data,
	date_start,
	date_end,
	which_retrieval='both',
	minute_avg=False,
	vers='v01',
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the level 2a data time
	series of each day so that you'll have one dictionary, whose e.g. 'IWV' will contain the IWV
	for the entire date range period.

	Parameters:
	-----------
	path_data : str
		Base path of level 2a data. This directory contains subfolders representing the year, which,
		in turn, contain months, which contain day subfolders. Example path_data:
		"/data/obs/campaigns/mosaic/hatpro/l2/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	which_retrieval : str, optional
		This describes which variable(s) will be loaded. Options: 'iwv' or 'prw' will load the
		integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. 'both' will 
		load both. Default: 'both'
	minute_avg : bool
		If True: averages over one minute are computed and returned instead of False when all
		data points are returned (more outliers, higher memory usage).
	vers : str
		Indicates the mwr_pro output version number. Valid options: 'i01', 'v00', and 'v01'. Default: 'v01'
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	if vers not in ['i01', 'v00', 'v01']:
		raise ValueError("In import_hatpro_level2a_daterange, the argument 'vers' must be one of the" +
							" following options: 'i01', 'v00', 'v01'")

	# check if the input of the retrieval variable is okay:
	if not isinstance(which_retrieval, str):
			raise TypeError("Argument 'which_retrieval' must be a string. Options: 'iwv' or 'prw' will load the " +
				"integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. " +
				"'both' will load both. Default: 'both'")

	else:
		if which_retrieval not in ['prw', 'iwv', 'clwvi', 'lwp', 'both']:
			raise ValueError("Argument 'which_retrieval' must be one of the following options: 'iwv' or 'prw' will load the " +
				"integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. " +
				"'both' will load both. Default: 'both'")

		else:
				if which_retrieval == 'iwv':
					which_retrieval = ['prw']
				elif which_retrieval == 'lwp':
					which_retrieval = ['clwvi']
				elif which_retrieval == 'both':
					which_retrieval = ['prw', 'clwvi']
				else:
					raise ValueError("Argument '" + which_retrieval + "' not recognised. Please use one of the following options: " +
						"'iwv' or 'prw' will load the " +
						"integrated water vapour. 'lwp' or 'clwvi' will load the liquid water path. " +
						"'both' will load both. Default: 'both'")
					

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")


	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1
	n_ret = 86			# inquired from level 2a data, number of available elevation angles in retrieval

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'flag', 'lat', 'lon', 'zsl']				# keys with time as coordinate

	# mwr_master_dict (output) will contain all desired variables on specific axes:
	# e.g. level 2a and 2b have got the same time axis (according to pl_mk_nds.pro)
	# and azimuth and elevation angles.
	mwr_master_dict = dict()
	if minute_avg:	# max number of minutes: n_days*1440
		n_minutes = n_days*1440
		for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_minutes,), np.nan)

		if 'prw' in which_retrieval:
			mwr_master_dict['prw'] = np.full((n_minutes,), np.nan)
			mwr_master_dict['prw_offset'] = np.full((n_minutes,), np.nan)							####### could be reduced to a daily value or even one value for the entire period
			# mwr_master_dict['prw_err'] = np.full((n_days, n_ret), np.nan)								####### could be reduced to one value array for the entire period
			mwr_master_dict['prw_err'] = np.full((n_ret,), np.nan)

		if 'clwvi' in which_retrieval:
			mwr_master_dict['clwvi'] = np.full((n_minutes,), np.nan)
			mwr_master_dict['clwvi_offset'] = np.full((n_minutes,), np.nan)							####### could be reduced to a daily value or even one value for the entire period
			# mwr_master_dict['clwvi_err'] = np.full((n_days, n_ret), np.nan)								###### could be reduced to one value for the entire period
			mwr_master_dict['clwvi_err'] = np.full((n_ret,), np.nan)

	else:			# max number of seconds: n_days*86400
		n_seconds = n_days*86400
		for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_seconds,), np.nan)

		if 'prw' in which_retrieval:
			mwr_master_dict['prw'] = np.full((n_seconds,), np.nan)
			mwr_master_dict['prw_offset'] = np.full((n_seconds,), np.nan)							####### could be reduced to a daily value or even one value for the entire period
			# mwr_master_dict['prw_err'] = np.full((n_days, n_ret), np.nan)								####### could be reduced to one value array for the entire period
			mwr_master_dict['prw_err'] = np.full((n_ret,), np.nan)

		if 'clwvi' in which_retrieval:
			mwr_master_dict['clwvi'] = np.full((n_seconds,), np.nan)
			mwr_master_dict['clwvi_offset'] = np.full((n_seconds,), np.nan)							####### could be reduced to a daily value or even one value for the entire period
			# mwr_master_dict['clwvi_err'] = np.full((n_days, n_ret), np.nan)								###### could be reduced to one value for the entire period
			mwr_master_dict['clwvi_err'] = np.full((n_ret,), np.nan)


	# cycle through all years, all months and days:
	time_index = 0	# this index (for lvl 2a) will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	day_index = 0	# same as above, but only increases by 1 for each day
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on HATPRO Level 2a, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# list of files:
		hatpro_level2_nc = sorted(glob.glob(day_path + "*.nc"))
		# filter for v01 files:
		hatpro_level2_nc = [lvl2_nc for lvl2_nc in hatpro_level2_nc if vers in lvl2_nc]

		if len(hatpro_level2_nc) == 0:
			if verbose >= 2:
				warnings.warn("No netcdf files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# identify level 2a files:
		hatpro_level2a_nc = []
		for lvl2_nc in hatpro_level2_nc:
			for wr in which_retrieval:
				if wr + '_' in lvl2_nc:
					hatpro_level2a_nc.append(lvl2_nc)

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for lvl2_nc in hatpro_level2a_nc: 
			mwr_dict = import_hatpro_level2a(lvl2_nc, minute_avg=minute_avg)

			n_time = len(mwr_dict['time'])
			cur_time_shape = mwr_dict['time'].shape

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				mwr_key_shape = mwr_dict[mwr_key].shape
				if mwr_key_shape == cur_time_shape:	# then the variable is on time axis:
					mwr_master_dict[mwr_key][time_index:time_index + n_time] = mwr_dict[mwr_key]

				elif mwr_key == 'prw_err' or mwr_key == 'clwvi_err': 	# these variables are nret x 1 arrays
					# mwr_master_dict[mwr_key][day_index:day_index + 1, :] = mwr_dict[mwr_key]			## for the case that we leave _err a daily value
					mwr_master_dict[mwr_key][:] = mwr_dict[mwr_key]

				elif len(mwr_dict[mwr_key]) == 1:
					mwr_master_dict[mwr_key][day_index:day_index + 1] = mwr_dict[mwr_key]

				else:
					raise RuntimeError("Something went wrong in the " +
						"import_hatpro_level2a_daterange routine. Unexpected MWR variable dimension. " + 
						"The length of one used variable ('%s') of level 2a data "%(mwr_key) +
							"neither equals the length of the time axis nor equals 1.")


		time_index = time_index + n_time
		day_index = day_index + 1

	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))
	else:

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		for mwr_key in mwr_master_dict.keys():
			if mwr_master_dict[mwr_key].shape == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]

	return mwr_master_dict


def import_hatpro_level2b_daterange(
	path_data,
	date_start,
	date_end,
	which_retrieval='both',
	vers='v01',
	around_radiosondes=True,
	path_radiosondes="",
	s_version='level_2',
	mwr_avg=0,
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the level 2b data time
	series of each day so that you'll have one dictionary, whose e.g. 'ta' will contain the
	temperature profile for the entire date range period with samples around the radiosonde
	launch times or alternatively 4 samples per day at fixed times: 05, 11, 17 and 23 UTC.

	Parameters:
	-----------
	path_data : str
		Base path of level 2b data. This directory contains subfolders representing the year, which,
		in turn, contain months, which contain day subfolders. Example path_data:
		"/data/obs/campaigns/mosaic/hatpro/l2/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	which_retrieval : str, optional
		This describes which variable(s) will be loaded. Options: 'ta' or 'hus' will load either the
		temperature or the specific humidity profile. 'both' will load both. Default: 'both'
	vers : str, optional
		Indicates the mwr_pro output version number. Valid options: 'i01', 'v00', and 'v01'. Default: 'v01'
	around_radiosondes : bool, optional
		If True, data will be limited to the time around radiosonde launches. If False, something else
		(e.g. around 4 times a day) might be done. Default: True
	path_radiosondes : str, optional
		Path to radiosonde data (Level 2). Default: ""
	s_version : str, optional
		Specifies the radiosonde version that is to be imported. Must be 'level_2' to work properly.
		Other versions have not been implemeted because they are considered to be inferior to level_2
		radiosondes.
	mwr_avg : int, optional
		If > 0, an average over mwr_avg seconds will be performed from sample_time to sample_time + 
		mwr_avg seconds. If == 0, no averaging will be performed.
	verbose : int, optional
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	if vers not in ['i01', 'v00', 'v01']:
		raise ValueError("In import_hatpro_level2b_daterange, the argument 'vers' must be one of the" +
							" following options: 'i01', 'v00', 'v01'")

	if mwr_avg < 0:
		raise ValueError("mwr_avg must be an int >= 0.")
	elif type(mwr_avg) != type(1):
		raise TypeError("mwr_avg must be int.")

	# check if the input of the retrieval variable is okay:
	if not isinstance(which_retrieval, str):
		raise TypeError("Argument 'which_retrieval' must be a string. Options: 'ta' or 'hus' will load either the " +
			"temperature or the specific humidity profile. 'both' will load both. Default: 'both'")

	elif which_retrieval not in ['ta', 'hus', 'both']:
		raise ValueError("Argument 'which_retrieval' must be one of the following options: 'ta' or 'hus' will load either the " +
			"temperature or the specific humidity profile. 'both' will load both. Default: 'both'")

	else:
		which_retrieval_dict = {'ta': ['ta'],
								'hus': ['hus'],
								'both': ['ta', 'hus']}
		level2b_dataID_dict = {'ta': ['ta'],
								'hus': ['hua'],
								'both': ['ta', 'hua']}
		level2b_dataID = level2b_dataID_dict[which_retrieval]			# to find correct file names
		which_retrieval = which_retrieval_dict[which_retrieval]


	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")

	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1
	n_ret = 1			# inquired from level 2b data, number of available elevation angles in retrieval
	n_hgt = 43			# inquired from level 2b data, number of vertical retrieval levels (height levels)

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'flag', 'lat', 'lon', 'zsl']				# keys with time as coordinate
	mwr_height_keys = ['height']							# keys with height as coordinate

	# Create an array that includes the radiosonde launch times:
	if around_radiosondes:
		if not path_radiosondes:
			raise ValueError("If 'around_radiosondes' is True, the path to the radiosonde level 2 data ('pathradiosondes') " +
								"must be given.")

		if s_version != 'level_2':
			raise ValueError("Radiosonde version 's_version' must be 'level_2' if around_radiosondes is True because " +
								"for this version, the launch time is directly read from the filename. This has not " +
								"been implemeted for other radiosonde versions ('mossonde', 'psYYMMDDwHH') because these " +
								"are considered to be inferior.")
		else:
			add_files = sorted(glob.glob(path_radiosondes + "*.nc"))		# filenames only; filter path
			add_files = [os.path.basename(a_f) for a_f in add_files]
			
			# identify launch time:
			n_samp = len(add_files)		# number of radiosondes
			launch_times = np.full((n_samp,), dt.datetime(1970,1,1))
			kk = 0
			for a_f in add_files:
				ltt = dt.datetime.strptime(a_f[-19:-4], "%Y%m%d_%H%M%S")
				# only save those that are in the considered period
				if ltt >= date_start and ltt < (date_end + dt.timedelta(days=1)):
					launch_times[kk] = ltt
					kk += 1
			
			# truncate launch_times and convert to sec since 1970-01-01 00:00:00 UTC:
			launch_times = launch_times[:kk]
			sample_times = datetime_to_epochtime(launch_times)
			n_samp_tot = len(sample_times)

	else:
		# max number of samples: n_days*4
		sample_times = [5, 11, 17, 23]		# UTC on each day
		n_samp = len(sample_times)
		n_samp_tot = n_days*n_samp

	# mwr_master_dict (output) will contain all desired variables on specific axes:
	# e.g. level 2b has got a time axis (according to pl_mk_nds.pro) for flag,
	# azimuth, elevation angles and the data.
	mwr_master_dict = dict()

	# save import keys for each retrieval option in a dict:
	import_keys = dict()
	mwr_time_height_keys = []
	for l2b_ID in level2b_dataID: mwr_time_height_keys.append(l2b_ID)

	if 'ta' in which_retrieval:
		mwr_master_dict['ta_err'] = np.full((n_hgt, n_ret), np.nan)

		# define the keys that will be imported via import_hatpro_level2b:
		import_keys['ta'] = (mwr_time_keys + mwr_height_keys +
						['ta', 'ta_err'])

	if 'hus' in which_retrieval:
		# here, we can only import and concat absolute humidity (hua) because
		# the conversion requires temperature and pressure
		mwr_master_dict['hua_err'] = np.full((n_hgt, n_ret), np.nan)

		# define the keys that will be imported via import_hatpro_level2b:
		import_keys['hua'] = (mwr_time_keys + mwr_height_keys +
						['hua', 'hua_err'])

	for mthk in mwr_time_height_keys: mwr_master_dict[mthk] = np.full((n_samp_tot, n_hgt), np.nan)
	for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_samp_tot,), np.nan)
	for mhk in mwr_height_keys: mwr_master_dict[mhk] = np.full((n_hgt,), np.nan)

	# cycle through all years, all months and days:
	time_index = 0	# this index (for lvl 2b) will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	sample_time_tolerance = 900		# sample time tolerance in seconds: mwr time must be within this
									# +/- tolerance of a sample_time to be accepted
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on HATPRO Level 2b, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# specify sample times as time: sec since 1970-01-01 00:00:00 UTC:
		if around_radiosondes:
			now_date_date = now_date.date()
			sample_mask = np.full((n_samp_tot,), False)
			for kk, l_t in enumerate(launch_times):
				sample_mask[kk] = l_t.date() == now_date_date

			sample_times_t = sample_times[sample_mask]

		else:
			sample_times_t = np.asarray([datetime_to_epochtime(dt.datetime(yyyy, mm, dd, st, 0, 0)) for st in sample_times])

		# list of v01 files:
		hatpro_level2_nc = sorted(glob.glob(day_path + "*_%s_*.nc"%vers))

		if len(hatpro_level2_nc) == 0:
			if verbose >= 2:
				warnings.warn("No netcdf files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# identify level 2b files:
		# also save the dataID into the list to access the correct keys to be imported (import_keys)
		# later on.
		hatpro_level2b_nc = []
		for lvl2_nc in hatpro_level2_nc:
			for dataID in level2b_dataID:
				# must avoid including the boundary layer scan
				if (dataID + '_' in lvl2_nc) and ('BL00_' not in lvl2_nc):
					hatpro_level2b_nc.append([lvl2_nc, dataID])

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for lvl2_nc in hatpro_level2b_nc:
			mwr_dict = import_hatpro_level2b(lvl2_nc[0], import_keys[lvl2_nc[1]])

			# it may occur that the whole day is flagged. If so, skip this file:
			if not np.any(mwr_dict['flag'] == 0):
				n_samp_real = 0
				continue

			# remove values where flag > 0:
			for mthk in mwr_time_height_keys:
				if mthk in lvl2_nc[1]:
					mwr_dict[mthk] = mwr_dict[mthk][mwr_dict['flag'] == 0,:]
			for mtkab in mwr_time_keys:
				if mtkab != 'flag':
					mwr_dict[mtkab] = mwr_dict[mtkab][mwr_dict['flag'] == 0]
			mwr_dict['flag'] = mwr_dict['flag'][mwr_dict['flag'] == 0]

			# # # update the flag by taking the manually detected outliers into account:
			# # # (not needed if v01 or later is used)
			# mwr_dict['flag'] = outliers_per_eye(mwr_dict['flag'], mwr_dict['time'], instrument='hatpro')

			# find the time slice where the mwr time is closest to the sample_times.
			# The identified index must be within 15 minutes, otherwise it will be discarded
			# Furthermore, it needs to be respected, that the flag value must be 0 for that case.
			if mwr_avg == 0:
				sample_idx = []
				for st in sample_times_t:
					idx = np.argmin(np.abs(mwr_dict['time'] - st))
					if np.abs(mwr_dict['time'][idx] - st) < sample_time_tolerance:
						sample_idx.append(idx)
				sample_idx = np.asarray(sample_idx)
				n_samp_real = len(sample_idx)	# number of samples that are valid to use; will be equal to n_samp in most cases

			else:
				sample_idx = []
				for st in sample_times_t:
					idx = np.where((mwr_dict['time'] >= st) & (mwr_dict['time'] <= st + mwr_avg))[0]
					if len(idx) > 0:	# then an overlap has been found
						sample_idx.append(idx)
				sample_idx = np.asarray(sample_idx)
				n_samp_real = len(sample_idx)	# number of samples that are valid to use; will be equal to n_samp in most cases

			if n_samp_real == 0: continue

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				mwr_key_shape = mwr_dict[mwr_key].shape

				if (mwr_key_shape == mwr_dict['time'].shape) and (mwr_key in mwr_time_keys):	# then the variable is on time axis:
					if mwr_avg > 0:				# these values won't be averaged because they don't contain "data"
						sample_idx_idx = [sii[0] for sii in sample_idx]
						mwr_master_dict[mwr_key][time_index:time_index + n_samp_real] = mwr_dict[mwr_key][sample_idx_idx]
					
					else:
						mwr_master_dict[mwr_key][time_index:time_index + n_samp_real] = mwr_dict[mwr_key][sample_idx]

				elif mwr_key == 'hua_err' or mwr_key == 'ta_err': 	# these variables are n_hgt x n_ret arrays
					mwr_master_dict[mwr_key] = mwr_dict[mwr_key]

				elif mwr_key in mwr_height_keys:	# handled after the for loop
					continue

				elif mwr_key in mwr_time_height_keys:
					if mwr_avg > 0:
						for k, sii in enumerate(sample_idx):
							mwr_master_dict[mwr_key][time_index+k:time_index+k + 1,:] = np.nanmean(mwr_dict[mwr_key][sii,:], axis=0)
					else:
						mwr_master_dict[mwr_key][time_index:time_index + n_samp_real,:] = mwr_dict[mwr_key][sample_idx,:]

				else:
					raise RuntimeError("Something went wrong in the " +
						"import_hatpro_level2b_daterange routine. Unexpected MWR variable dimension for " + mwr_key + ".")


		time_index = time_index + n_samp_real

	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))
	else:
		# save non height dependent variables to master dict:
		for mwr_key in mwr_height_keys: mwr_master_dict[mwr_key] = mwr_dict[mwr_key]

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		time_height_shape_old = mwr_master_dict[mwr_time_height_keys[0]].shape
		for mwr_key in mwr_master_dict.keys():
			shape_new = mwr_master_dict[mwr_key].shape
			if shape_new == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]
			elif shape_new == time_height_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1, :]

	return mwr_master_dict


def import_hatpro_level2c_daterange(
	path_data,
	date_start,
	date_end,
	which_retrieval='both',
	vers='v01',
	around_radiosondes=True,
	path_radiosondes="",
	s_version='level_2',
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the level 2c data time
	series of each day so that you'll have one dictionary, whose e.g. 'ta' will contain the
	temperature profile for the entire date range period with samples around the radiosonde
	launch times or alternatively 4 samples per day at fixed times: 05, 11, 17 and 23 UTC.

	Parameters:
	-----------
	path_data : str
		Base path of level 2c data. This directory contains subfolders representing the year, which,
		in turn, contain months, which contain day subfolders. Example path_data:
		"/data/obs/campaigns/mosaic/hatpro/l2/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	which_retrieval : str, optional
		This describes which variable(s) will be loaded. Options: 'ta' will load the temperature 
		profile (boundary layer scan). 'both' will also load temperature only because humidity profile
		boundary layer scan does not exist. Default: 'both'
	vers : str
		Indicates the mwr_pro output version number. Valid options: 'i01', 'v00', and 'v01'. Default: 'v01'
	around_radiosondes : bool, optional
		If True, data will be limited to the time around radiosonde launches. If False, something else
		(e.g. around 4 times a day) might be done. Default: True
	path_radiosondes : str, optional
		Path to radiosonde data (Level 2). Default: ""
	s_version : str
		Specifies the radiosonde version that is to be imported. Must be 'level_2' to work properly.
		Other versions have not been implemeted because they are considered to be inferior to level_2
		radiosondes.
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	if vers not in ['i01', 'v00', 'v01']:
		raise ValueError("In import_hatpro_level2c_daterange, the argument 'vers' must be one of the" +
							" following options: 'i01', 'v00', 'v01'")

	# check if the input of the retrieval variable is okay:
	if not isinstance(which_retrieval, str):
		raise TypeError("Argument 'which_retrieval' must be a string. Options: 'ta' will load the temperature " +
		"profile (boundary layer scan). 'both' will also load temperature only because humidity profile" +
		"boundary layer scan does not exist. Default: 'both'")

	elif which_retrieval not in ['ta', 'hus', 'both']:
		raise ValueError("Argument 'which_retrieval' must be one of the following options: 'ta' will load the temperature " +
		"profile (boundary layer scan). 'both' will also load temperature only because humidity profile" +
		"boundary layer scan does not exist. Default: 'both'")

	else:
		which_retrieval_dict = {'ta': ['ta'],
								'both': ['ta']}
		level2c_dataID_dict = {'ta': ['ta'],
								'both': ['ta']}
		level2c_dataID = level2c_dataID_dict[which_retrieval]
		which_retrieval = which_retrieval_dict[which_retrieval]

	# check if around_radiosondes is the right type:
	if not isinstance(around_radiosondes, bool):
		raise TypeError("Argument 'around_radiosondes' must be either True or False (boolean type).")

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")

	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1
	n_ret = 1			# inquired from level 2c data, number of available elevation angles in retrieval
	n_hgt = 43			# inquired from level 2c data, number of vertical retrieval levels (height levels)

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'flag', 'lat', 'lon', 'zsl']				# keys with time as coordinate
	mwr_height_keys = ['height']						# keys with height as coordinate

	# Create an array that includes the radiosonde launch times:
	if around_radiosondes:
		if not path_radiosondes:
			raise ValueError("If 'around_radiosondes' is True, the path to the radiosonde level 2 data ('pathradiosondes') " +
								"must be given.")

		if s_version != 'level_2':
			raise ValueError("Radiosonde version 's_version' must be 'level_2' if around_radiosondes is True because " +
								"for this version, the launch time is directly read from the filename. This has not " +
								"been implemeted for other radiosonde versions ('mossonde', 'psYYMMDDwHH') because these " +
								"are considered to be inferior.")
		else:
			add_files = sorted(glob.glob(path_radiosondes + "*.nc"))		# filenames only; filter path
			add_files = [os.path.basename(a_f) for a_f in add_files]
			
			# identify launch time:
			n_samp = len(add_files)		# number of radiosondes
			launch_times = np.full((n_samp,), dt.datetime(1970,1,1))
			kk = 0
			for a_f in add_files:
				ltt = dt.datetime.strptime(a_f[-19:-4], "%Y%m%d_%H%M%S")
				# only save those that are in the considered period
				if ltt >= date_start and ltt < (date_end + dt.timedelta(days=1)):
					launch_times[kk] = ltt
					kk += 1
			
			# truncate launch_times and convert to sec since 1970-01-01 00:00:00 UTC:
			launch_times = launch_times[:kk]
			sample_times = datetime_to_epochtime(launch_times)
			n_samp_tot = len(sample_times)

	else:
		# max number of samples: n_days*4
		sample_times = [5, 11, 17, 23]		# UTC on each day
		n_samp = len(sample_times)
		n_samp_tot = n_days*n_samp

	# mwr_master_dict (output) will contain all desired variables on specific axes:
	# e.g. level 2c has got a time axis (according to pl_mk_nds.pro) for flag,
	# and the data.
	mwr_master_dict = dict()

	# save import keys for each retrieval option in a dict:
	import_keys = dict()
	mwr_time_height_keys = []
	for l2b_ID in level2c_dataID: mwr_time_height_keys.append(l2b_ID)

	if 'ta' in which_retrieval:
		mwr_master_dict['ta_err'] = np.full((n_hgt,), np.nan)

		# define the keys that will be imported via import_hatpro_level2b:
		import_keys['ta'] = (mwr_time_keys + mwr_height_keys +
						['ta', 'ta_err'])

	for mthk in mwr_time_height_keys: mwr_master_dict[mthk] = np.full((n_samp_tot, n_hgt), np.nan)
	for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_samp_tot,), np.nan)
	for mhk in mwr_height_keys: mwr_master_dict[mhk] = np.full((n_hgt,), np.nan)

	# cycle through all years, all months and days:
	time_index = 0	# this index (for lvl 2c) will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	sample_time_tolerance = 1800		# sample time tolerance in seconds: mwr time must be within this
										# +/- tolerance of a sample_time to be accepted
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on HATPRO Level 2c, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# specify sample times as time: sec since 1970-01-01 00:00:00 UTC:
		if around_radiosondes:
			now_date_date = now_date.date()
			sample_mask = np.full((n_samp_tot,), False)
			for kk, l_t in enumerate(launch_times):
				sample_mask[kk] = l_t.date() == now_date_date

			sample_times_t = sample_times[sample_mask]

		else:
			sample_times_t = np.asarray([datetime_to_epochtime(dt.datetime(yyyy, mm, dd, st, 0, 0)) for st in sample_times])

		# list of v01 files:
		hatpro_level2_nc = sorted(glob.glob(day_path + "*_%s_*.nc"%vers))

		if len(hatpro_level2_nc) == 0:
			if verbose >= 2:
				warnings.warn("No netcdf files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# identify level 2c files:
		# also save the dataID into the list to access the correct keys to be imported (import_keys)
		# later on.
		hatpro_level2c_nc = []
		for lvl2_nc in hatpro_level2_nc:
			for dataID in level2c_dataID:
				# must include the boundary layer scan
				if (dataID + '_' in lvl2_nc) and ('BL00_' in lvl2_nc):
					hatpro_level2c_nc.append([lvl2_nc, dataID])

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for lvl2_nc in hatpro_level2c_nc:
			mwr_dict = import_hatpro_level2c(lvl2_nc[0], import_keys[lvl2_nc[1]])

			# it may occur that the whole day is flagged. If so, skip this file:
			if not np.any(mwr_dict['flag'] == 0):
				n_samp_real = 0
				continue

			# remove values where flag > 0:
			for mthk in mwr_time_height_keys: mwr_dict[mthk] = mwr_dict[mthk][mwr_dict['flag'] == 0,:]
			for mtkab in mwr_time_keys:
				if mtkab != 'flag':
					mwr_dict[mtkab] = mwr_dict[mtkab][mwr_dict['flag'] == 0]
			mwr_dict['flag'] = mwr_dict['flag'][mwr_dict['flag'] == 0]


			# # # update the flag by taking the manually detected outliers into account:
			# # # (not needed if v01 or later is used)
			# mwr_dict['flag'] = outliers_per_eye(mwr_dict['flag'], mwr_dict['time'], instrument='hatpro')

			# find the time slice where the mwr time is closest to the sample_times.
			# The identified index must be within 30 minutes, otherwise it will be discarded.
			# Furthermore, it needs to be respected, that the flag value must be 0 for that case.
			sample_idx = []
			for st in sample_times_t:
				idx = np.argmin(np.abs(mwr_dict['time'] - st))
				if np.abs(mwr_dict['time'][idx] - st) < sample_time_tolerance:
					sample_idx.append(idx)
			sample_idx = np.asarray(sample_idx)
			n_samp_real = len(sample_idx)	# number of samples that are valid to use; will be equal to n_samp in most cases

			if n_samp_real == 0: continue

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				mwr_key_shape = mwr_dict[mwr_key].shape

				if (mwr_key_shape == mwr_dict['time'].shape) and (mwr_key in mwr_time_keys):	# then the variable is on time axis:
					mwr_master_dict[mwr_key][time_index:time_index + n_samp_real] = mwr_dict[mwr_key][sample_idx]

				elif mwr_key == 'ta_err': 	# these variables are n_hgt x n_ret arrays
					mwr_master_dict[mwr_key] = mwr_dict[mwr_key]

				elif mwr_key in mwr_height_keys: # handled after the for loop
					continue

				elif mwr_key in mwr_time_height_keys:
					# first: filter for non-flagged values
					mwr_master_dict[mwr_key][time_index:time_index + n_samp_real,:] = mwr_dict[mwr_key][sample_idx,:]

				else:
					raise RuntimeError("Something went wrong in the " +
						"import_hatpro_level2c_daterange routine. Unexpected MWR variable dimension for %s."%mwr_key)


		time_index = time_index + n_samp_real

	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))
	else:
		# save non time dependent variables in master dict
		for mwr_key in mwr_height_keys: mwr_master_dict[mwr_key] = mwr_dict[mwr_key]

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		time_height_shape_old = mwr_master_dict[mwr_time_height_keys[0]].shape
		for mwr_key in mwr_master_dict.keys():
			shape_new = mwr_master_dict[mwr_key].shape
			if shape_new == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]
			elif shape_new == time_height_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1, :]

	return mwr_master_dict


def import_mirac_MET(
	filename,
	keys='basic'):

	"""
	Importing automatically created MiRAC-P MET hourly files
	with the ending .MET.NC in the level 1 folder. Time will be
	converted to seconds since 1970-01-01 00:00:00 UTC.
	For now, only surface pressure (in Pa) will be loaded.

	Parameters:
	-----------
	filename : str
		Path and filename of MiRAC-P .MET.NC data.
	keys : list of str or str, optional
		Specify which variables are to be imported. Another option is
		to import all keys (keys='all') or import basic keys
		that the author considers most important (keys='basic')
		or leave this argument out.
	"""


	file_nc = nc.Dataset(filename)

	if keys == 'basic':
		keys = ['time', 'RF', 'Surf_P']

	elif keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and ((keys != 'all') and (keys != 'basic')):
		raise ValueError("Argument 'keys' must either be a string ('all' or 'basic') or a list of variable names.")

	mwr_dict = dict()
	reftime = dt.datetime(1970,1,1,0,0,0)
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in MiRAC-P .MET.NC file." % key)

		mwr_dict[key] = np.asarray(file_nc.variables[key])

		if key == 'time':	# convert to sec since 1970-01-01 00:00:00 UTC (USE FLOAT64)
			mwr_dict['time'] = (np.float64(datetime_to_epochtime(dt.datetime(2001,1,1,0,0,0))) +
								mwr_dict[key].astype(np.float64))

	if "Surf_P" in mwr_dict.keys():
		mwr_dict['pres'] = mwr_dict['Surf_P']*100		# convert to Pa

	return mwr_dict


def import_mirac_MET_RPG_daterange(
	path_data,
	date_start,
	date_end,
	verbose=0):

	"""
	Runs through all days between a start and an end date. It concats the MET time
	series of each day so that you'll have one dictionary that will contain the pressure data
	for the entire date range period.

	Parameters:
	-----------
	path_data : str
		Base path of MiRAC-P level 1 data. This directory contains subfolders representing the 
		year, which, in turn, contain months, which contain day subfolders. Example:
		path_data = "/data/obs/campaigns/mosaic/mirac-p/l1/"
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")

	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days + 1

	# basic variables that should always be imported:
	mwr_time_keys = ['time', 'RF', 'pres']	# keys with time as coordinate

	# mwr_master_dict (output) will contain all desired variables on time axis for entire date range:
	mwr_master_dict = dict()
	n_seconds = n_days*86400		# max number of seconds: n_days*86400
	for mtkab in mwr_time_keys: mwr_master_dict[mtkab] = np.full((n_seconds,), np.nan)


	# Load the data into mwr_master_dict:
	# cycle through all years, all months and days:
	time_index = 0	# this index will be increased by the length of the time
						# series of the current day (now_date) to fill the mwr_master_dict time axis
						# accordingly.
	day_index = 0	# will increase for each day
	for now_date in (date_start + dt.timedelta(days=n) for n in range(n_days)):

		if verbose >= 1: print("Working on RPG retrieval, MiRAC-P MET, ", now_date)

		yyyy = now_date.year
		mm = now_date.month
		dd = now_date.day

		day_path = path_data + "%04i/%02i/%02i/"%(yyyy,mm,dd)

		if not os.path.exists(os.path.dirname(day_path)):
			continue

		# list of .MET.NC files: Sorting is important as this will
		# ensure automatically that the time series of each hour will
		# be concatenated appropriately!
		mirac_nc = sorted(glob.glob(day_path + "*.MET.NC"))
		if len(mirac_nc) == 0:
			if verbose >= 2:
				warnings.warn("No .MET.NC files found for date %04i-%02i-%02i."%(yyyy,mm,dd))
			continue

		# load one retrieved variable after another from current day and save it into the mwr_master_dict
		for l1_file in mirac_nc: 
			mwr_dict = import_mirac_MET(l1_file)

			n_time = len(mwr_dict['time'])
			time_shape = mwr_dict['time'].shape

			# save to mwr_master_dict
			for mwr_key in mwr_dict.keys():
				if mwr_key in mwr_time_keys:
					mwr_master_dict[mwr_key][time_index:time_index + n_time] = mwr_dict[mwr_key]


			time_index = time_index + n_time
		day_index = day_index + 1


	if time_index == 0 and verbose >= 1: 	# otherwise no data has been found
		raise ValueError("No data found in date range: " + dt.datetime.strftime(date_start, "%Y-%m-%d") + " - " + 
				dt.datetime.strftime(date_end, "%Y-%m-%d"))

	else:

		# truncate the mwr_master_dict to the last nonnan time index:
		last_time_step = np.argwhere(~np.isnan(mwr_master_dict['time']))[-1][0]
		time_shape_old = mwr_master_dict['time'].shape
		for mwr_key in mwr_master_dict.keys():
			shape_new = mwr_master_dict[mwr_key].shape
			if shape_new == time_shape_old:
				mwr_master_dict[mwr_key] = mwr_master_dict[mwr_key][:last_time_step+1]

	return mwr_master_dict


def import_single_PS122_mosaic_radiosonde_level2(
	filename,
	keys='all',
	verbose=0):

	"""
	Imports single level 2 radiosonde data created with PANGAEA_tab_to_nc.py 
	('PS122_mosaic_radiosonde_level2_yyyymmdd_hhmmssZ.nc'). Converts to SI units
	and interpolates to a height grid with 5 m resolution from 0 to 15000 m. 

	Parameters:
	-----------
	filename : str
		Name (including path) of radiosonde data file.
	keys : list of str or str, optional
		This describes which variable(s) will be loaded. Specifying 'all' will import all variables.
		Specifying 'basic' will load the variables the author consideres most useful for his current
		analysis.
		Default: 'all'
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	"""
		Loaded values are imported in the following units:
		T: in deg C, will be converted to K
		P: in hPa, will be converted to Pa
		RH: in %, will be converted to [0-1]
		Altitude: in m
		q: in kg kg^-1 (water vapor specific humidity)
		time: in sec since 1970-01-01 00:00:00 UTC
	"""

	file_nc = nc.Dataset(filename)

	if (not isinstance(keys, str)) and (not isinstance(keys, list)):
		raise TypeError("Argument 'key' must be a list of strings or 'all'.")

	if keys == 'all':
		keys = file_nc.variables.keys()
	elif keys == 'basic':
		keys = ['time', 'T', 'P', 'RH', 'q', 'Altitude']

	sonde_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in radiosonde file." % key)

		sonde_dict[key] = np.asarray(file_nc.variables[key])
		if key != "IWV" and len(sonde_dict[key]) == 0: # 'and': second condition only evaluated if first condition True
			return None

		if key in ['Latitude', 'Longitude']:	# only interested in the first lat, lon position
			sonde_dict[key] = sonde_dict[key][0]
		if key == 'IWV':
			sonde_dict[key] = np.float64(sonde_dict[key])

	# convert units:
	if 'RH' in keys:	# from percent to [0, 1]
		sonde_dict['RH'] = sonde_dict['RH']*0.01
	if 'T' in keys:		# from deg C to K
		sonde_dict['T'] = sonde_dict['T'] + 273.15
	if 'P' in keys:		# from hPa to Pa
		sonde_dict['P'] = sonde_dict['P']*100
	if 'time' in keys:	# from int64 to float64
		sonde_dict['time'] = np.float64(sonde_dict['time'])
		sonde_dict['launch_time'] = sonde_dict['time'][0]

	keys = [*keys]		# converts dict_keys to a list
	for key in keys:
		if sonde_dict[key].shape == sonde_dict['time'].shape:
			if key not in ['time', 'Latitude', 'Longitude', 'ETIM', 'Altitude']:
				sonde_dict[key + "_ip"] = np.interp(np.arange(0,15001,5), sonde_dict['Altitude'], sonde_dict[key], right=np.nan)
			elif key == 'Altitude':
				sonde_dict[key + "_ip"] = np.arange(0, 15001,5)


	# Renaming variables: ['Lat', 'Lon', 'p', 'T', 'RH', 'GeopHgt', 'qv', 'time', ...]
	renaming = {'T': 'temp', 	'P': 'pres', 	'RH': 'rh',
				'Altitude': 'height', 'h_geom': 'height_geom',
				'Latitude': 'lat', 	'Longitude': 'lon',
				'T_ip': 'temp_ip', 'P_ip': 'pres_ip', 'RH_ip': 'rh_ip',
				'Altitude_ip': 'height_ip', 'h_geom_ip': 'height_geom_ip',
				'IWV': 'iwv'}
	for ren_key in renaming.keys():
		if ren_key in sonde_dict.keys():
			sonde_dict[renaming[ren_key]] = sonde_dict[ren_key]

	# height check: how high does the data reach:
	sonde_dict['height_check'] = sonde_dict['height'][-1]

	return sonde_dict


def import_single_NYA_RS_radiosonde(
	filename,
	keys='all',
	verbose=0):

	"""
	Imports single NYA-RS radiosonde data for Ny Alesund. Converts to SI units
	and interpolates to a height grid with 5 m resolution from 0 to 15000 m. 

	Parameters:
	-----------
	filename : str
		Name (including path) of radiosonde data file.
	keys : list of str or str, optional
		This describes which variable(s) will be loaded. Specifying 'all' will import all variables.
		Specifying 'basic' will load the variables the author consideres most useful for his current
		analysis.
		Default: 'all'
	verbose : int
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	"""
		Loaded values are imported in the following units:
		T: in K
		P: in hPa, will be converted to Pa
		RH: in [0-1]
		Altitude: in m
		time: will be converted to sec since 1970-01-01 00:00:00 UTC
	"""

	file_nc = nc.Dataset(filename)

	if (not isinstance(keys, str)) and (not isinstance(keys, list)):
		raise TypeError("Argument 'key' must be a list of strings or 'all'.")

	if keys == 'all':
		keys = file_nc.variables.keys()
	elif keys == 'basic':
		keys = ['time', 'temp', 'press', 'rh', 'alt']

	sonde_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in radiosonde file." % key)

		sonde_dict[key] = np.asarray(file_nc.variables[key])
		if key != "IWV" and len(sonde_dict[key]) == 0: # 'and': second condition only evaluated if first condition True
			return None

		if key in ['lat', 'lon']:	# only interested in the first lat, lon position
			sonde_dict[key] = sonde_dict[key][0]

	# convert units:
	if 'P' in keys:		# from hPa to Pa
		sonde_dict['P'] = sonde_dict['P']*100
	if 'time' in keys:	# from int64 to float64
		time_unit = file_nc.variables['time'].units
		time_offset = (dt.datetime.strptime(time_unit[-19:], "%Y-%m-%dT%H:%M:%S") - dt.datetime(1970,1,1)).total_seconds()
		sonde_dict['time'] = np.float64(sonde_dict['time']) + time_offset
		sonde_dict['launch_time'] = sonde_dict['time'][0]

	keys = [*keys]		# converts dict_keys to a list
	for key in keys:
		if sonde_dict[key].shape == sonde_dict['time'].shape:
			if key not in ['time', 'lat', 'lon', 'alt']:
				sonde_dict[key + "_ip"] = np.interp(np.arange(0,15001,5), sonde_dict['alt'], sonde_dict[key])
			elif key == 'alt':
				sonde_dict[key + "_ip"] = np.arange(0, 15001,5)


	# Renaming variables to a standard convention
	renaming = {'press': 'pres', 'alt': 'height', 'press_ip': 'pres_ip', 'alt_ip': 'height_ip'}
	for ren_key in renaming.keys():
		if ren_key in sonde_dict.keys():
			sonde_dict[renaming[ren_key]] = sonde_dict[ren_key]

	return sonde_dict


def import_radiosonde_daterange(
	path_data,
	date_start,
	date_end,
	s_version='level_2',
	with_wind=False,
	remove_failed=False,
	verbose=0):

	"""
	Imports radiosonde data 'mossonde-curM1' and concatenates the files into time series x height.
	E.g. temperature profile will have the dimension: n_sondes x n_height

	Parameters:
	-----------
	path_data : str
		Path of radiosonde data.
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	s_version : str, optional
		Specifies the radiosonde version that is to be imported. Possible options: 'mossonde',
		'psYYMMDDwHH', 'level_2', 'nya-rs'. Default: 'level_2' (published by Marion Maturilli)
	with_wind : bool, optional
		This describes if wind measurements are included (True) or not (False). Does not work with
		s_version='psYYMMDDwHH'. Default: False
	remove_failed : bool, optional
		If True, failed sondes with unrealistic IWV values will be removed (currently only implmented
		for s_version == 'level_2'). It also includes "height_check" to avoid sondes that burst before
		reaching > 10000 m.
	verbose : int, optional
		If 0, output is suppressed. If 1, basic output is printed. If 2, more output (more warnings,...)
		is printed.
	"""

	if not isinstance(s_version, str): raise TypeError("s_version in import_radiosonde_daterange must be a string.")

	# extract day, month and year from start date:
	date_start = dt.datetime.strptime(date_start, "%Y-%m-%d")
	date_end = dt.datetime.strptime(date_end, "%Y-%m-%d")

	if s_version == 'mossonde':
		all_radiosondes_nc = sorted(glob.glob(path_data + "mossonde-curM1" + "*.nc"))

		# inquire the number of radiosonde files (date and time of launch is in filename):
		# And fill a list which will include the relevant radiosonde files.
		radiosondes_nc = []
		for rs_nc in all_radiosondes_nc:
			rs_date = rs_nc[-16:-8]		# date of radiosonde from filename
			yyyy = int(rs_date[:4])
			mm = int(rs_date[4:6])
			dd = int(rs_date[6:])
			rs_date_dt = dt.datetime(yyyy,mm,dd)
			if rs_date_dt >= date_start and rs_date_dt <= date_end:
				radiosondes_nc.append(rs_nc)

	elif s_version == 'psYYMMDDwHH':
		all_radiosondes_nc = sorted(glob.glob(path_data + "ps*.w*.nc"))[:-1]	# exclude last file because it's something about Ozone

		# inquire the number of radiosonde files (date and time of launch is in filename):
		# And fill a list which will include the relevant radiosonde files.
		radiosondes_nc = []
		for rs_nc in all_radiosondes_nc:
			rs_date = rs_nc[-13:-3]		# date of radiosonde from filename
			yyyy = 2000 + int(rs_date[:2])
			mm = int(rs_date[2:4])
			dd = int(rs_date[4:6])
			rs_date_dt = dt.datetime(yyyy,mm,dd)
			if rs_date_dt >= date_start and rs_date_dt <= date_end:
				radiosondes_nc.append(rs_nc)

	elif s_version == 'nya-rs':
		all_radiosondes_nc = sorted(glob.glob(path_data + "NYA-RS_*.nc"))

		# inquire the number of radiosonde files (date and time of launch is in filename):
		# And fill a list which will include the relevant radiosonde files.
		radiosondes_nc = []
		for rs_nc in all_radiosondes_nc:
			rs_date = rs_nc[-15:-3]		# date of radiosonde from filename
			yyyy = int(rs_date[:4])
			mm = int(rs_date[4:6])
			dd = int(rs_date[6:8])
			rs_date_dt = dt.datetime(yyyy,mm,dd)
			if rs_date_dt >= date_start and rs_date_dt <= date_end:
				radiosondes_nc.append(rs_nc)

	elif s_version == 'level_2':
		all_radiosondes_nc = sorted(glob.glob(path_data + "PS122_mosaic_radiosonde_level2*.nc"))

		# inquire the number of radiosonde files (date and time of launch is in filename):
		# And fill a list which will include the relevant radiosonde files.
		radiosondes_nc = []
		for rs_nc in all_radiosondes_nc:
			rs_date = rs_nc[-19:-3]		# date of radiosonde from filename
			yyyy = int(rs_date[:4])
			mm = int(rs_date[4:6])
			dd = int(rs_date[6:8])
			rs_date_dt = dt.datetime(yyyy,mm,dd)
			if rs_date_dt >= date_start and rs_date_dt <= date_end:
				radiosondes_nc.append(rs_nc)


	# number of sondes:
	n_sondes = len(radiosondes_nc)

	# count the number of days between start and end date as max. array size:
	n_days = (date_end - date_start).days

	# basic variables that should always be imported:
	if s_version == 'mossonde':
		geoinfo_keys = ['lat', 'lon', 'alt', 'launch_time']
		time_height_keys = ['pres', 'temp', 'rh', 'height', 'rho_v', 'q']		# keys with time and height as coordinate
		if with_wind: time_height_keys = time_height_keys + ['wspeed', 'wdir']
	elif s_version == 'psYYMMDDwHH':
		geoinfo_keys = ['lat', 'lon', 'launch_time']
		time_height_keys = ['pres', 'temp', 'rh', 'height', 'rho_v', 'q']
		if with_wind:
			print("No direct wind calculation for s_version='%s'."%s_version)
	elif s_version == 'nya-rs':
		geoinfo_keys = ['lat', 'lon', 'launch_time']
		time_height_keys = ['pres', 'temp', 'rh', 'height']
		if with_wind: time_height_keys = time_height_keys + ['wspeed', 'wdir']
	elif s_version == 'level_2':
		geoinfo_keys = ['lat', 'lon', 'launch_time', 'iwv']
		time_height_keys = ['pres', 'temp', 'rh', 'height', 'rho_v', 'q']
		if with_wind: time_height_keys = time_height_keys + ['wspeed', 'wdir']
	else:
		raise ValueError("s_version in import_radiosonde_daterange must be 'mossonde', 'psYYMMDDwHH', 'nya-rs', or 'level_2'.")
	all_keys = geoinfo_keys + time_height_keys

	# sonde_master_dict (output) will contain all desired variables on specific axes:
	# Time axis (one sonde = 1 timestamp) = axis 0; height axis = axis 1
	n_height = len(np.arange(0,15001,5))	# length of the interpolated height grid
	sonde_master_dict = dict()
	for gk in geoinfo_keys: sonde_master_dict[gk] = np.full((n_sondes,), np.nan)
	for thk in time_height_keys: sonde_master_dict[thk] = np.full((n_sondes, n_height), np.nan)

	if s_version == 'mossonde':
		all_keys_import = geoinfo_keys + time_height_keys + ['time', 'geopheight']	# 'time' required to create 'launch_time'
		all_keys_import.remove('launch_time')		# because this key is not saved in the radiosonde files
		all_keys_import.remove('rho_v')				# because this key is not saved in the radiosonde files
		all_keys_import.remove('q')					# because this key is not saved in the radiosonde files
		all_keys_import.remove('height')					# because this key is not saved in the radiosonde files
		if with_wind: all_keys_import = all_keys_import + ['wspeed', 'wdir']

		# cycle through all relevant sonde files:
		for rs_idx, rs_nc in enumerate(radiosondes_nc):

			if verbose >= 1:
				# rs_date = rs_nc[-16:-8]
				# print("Working on Radiosonde, ", 
					# dt.datetime(int(rs_date[:4]), int(rs_date[4:6]), int(rs_date[6:])))
				print("Working on Radiosonde, " + rs_nc)

			sonde_dict = import_single_mossonde_curM1(rs_nc, keys=all_keys_import)

			# save to sonde_master_dict:
			for key in all_keys:
				if key in geoinfo_keys:
					sonde_master_dict[key][rs_idx] = sonde_dict[key]

				elif key in time_height_keys:
					sonde_master_dict[key][rs_idx, :] = sonde_dict[key + "_ip"]		# must use the interpolated versions!

				else:
					raise KeyError("Key '" + key + "' not found in radiosonde dictionary after importing it with import_single_mossonde_curM1")

	elif s_version == 'psYYMMDDwHH':
		all_keys_import = ['Lat', 'Lon', 'p', 'T', 'RH', 'GeopHgt', 'qv', 'time']	# 'time' required to create 'launch_time'


		# cycle through all relevant sonde files:
		for rs_idx, rs_nc in enumerate(radiosondes_nc):

			if verbose >= 1: 
				# rs_date = rs_nc[-16:-8]
				print("Working on Radiosonde, " + rs_nc)

			sonde_dict = import_single_psYYMMDD_wHH_sonde(rs_nc, keys=all_keys_import)
			if not sonde_dict:	# then the imported sonde file appears to be empty
				continue

			else:
				# save to sonde_master_dict:
				for key in all_keys:
					if key in geoinfo_keys:
						sonde_master_dict[key][rs_idx] = sonde_dict[key]

					elif key in time_height_keys:
						sonde_master_dict[key][rs_idx, :] = sonde_dict[key + "_ip"]		# must use the interpolated versions!

					else:
						raise KeyError("Key '" + key + "' not found in radiosonde dictionary after importing it with import_single_mossonde_curM1")

		# As there are empty files among the current psYYMMDD.wHH sondes, they have to be filtered out:
		not_corrupted_sondes = ~np.isnan(sonde_master_dict['launch_time'])
		# not_corrupted_sondes_idx = np.where(~np.isnan(sonde_master_dict['launch_time']))[0]
		for key in sonde_master_dict.keys():
			if key in geoinfo_keys:
				sonde_master_dict[key] = sonde_master_dict[key][not_corrupted_sondes]
			else:
				sonde_master_dict[key] = sonde_master_dict[key][not_corrupted_sondes,:]

	elif s_version == 'nya-rs':
		all_keys_import = ['lat', 'lon', 'press', 'temp', 'rh', 'alt', 'time']
		if with_wind: all_keys_import = all_keys_import + ['wdir', 'wspeed']


		# cycle through all relevant sonde files:
		for rs_idx, rs_nc in enumerate(radiosondes_nc):
			
			if verbose >= 1:
				# rs_date = rs_nc[-19:-3]
				print("Working on Radiosonde, " + rs_nc)

			sonde_dict = import_single_NYA_RS_radiosonde(rs_nc, keys=all_keys_import)
			
			# save to sonde_master_dict:
			for key in all_keys:
				if key in geoinfo_keys:
					sonde_master_dict[key][rs_idx] = sonde_dict[key]

				elif key in time_height_keys:
					sonde_master_dict[key][rs_idx, :] = sonde_dict[key + "_ip"]		# must use the interpolated versions!

				else:
					raise KeyError("Key '" + key + "' not found in radiosonde dictionary after importing it with " +
									"import_single_NYA_RS_radiosonde")

	elif s_version == 'level_2':
		all_keys_import = ['Latitude', 'Longitude', 'P', 'T', 'RH', 'Altitude', 'rho_v', 'q', 'time', 'IWV']
		if with_wind: all_keys_import = all_keys_import + ['wdir', 'wspeed']

		if remove_failed:
			failed_sondes_t, failed_sondes_dt = time_prematurely_bursted_sondes()		# load times of failed sondes


		# cycle through all relevant sonde files:
		rs_idx = 0
		for rs_nc in radiosondes_nc:
			
			if verbose >= 1:
				# rs_date = rs_nc[-19:-3]
				print("Working on Radiosonde, " + rs_nc)

			sonde_dict = import_single_PS122_mosaic_radiosonde_level2(rs_nc, keys=all_keys_import)
			if (remove_failed and ((sonde_dict['iwv'] == 0.0) or (np.isnan(sonde_dict['iwv'])) or
				(sonde_dict['height_check'] < 10000) or (np.any(np.abs(sonde_dict['launch_time'] - failed_sondes_t) < 7200)))):
				continue
			
			# save to sonde_master_dict:
			for key in all_keys:
				if key in geoinfo_keys:
					sonde_master_dict[key][rs_idx] = sonde_dict[key]

				elif key in time_height_keys:
					sonde_master_dict[key][rs_idx, :] = sonde_dict[key + "_ip"]		# must use the interpolated versions!

				else:
					raise KeyError("Key '" + key + "' not found in radiosonde dictionary after importing it with " +
									"import_single_PS122_mosaic_radiosonde_level2")

			rs_idx += 1

		# Truncate number of sondes:
		if remove_failed and (rs_idx < n_sondes):
			for key in geoinfo_keys: sonde_master_dict[key] = sonde_master_dict[key][:rs_idx]
			for key in time_height_keys: sonde_master_dict[key] = sonde_master_dict[key][:rs_idx,:]

	return sonde_master_dict


def import_concat_IWV_LWP_mwr_master_time(
	filename,
	date_start,
	date_end):

	"""
	Simple importer to get the IWV or LWP data, that was stored on a 'master time axis',
	of a radiometer (HATPRO, MiRAC-P, ARM).

	Parameters:
	-----------
	filename : str
		Filename and path of the IWV or LWP data of radiometers on the master time axis.
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	"""

	date_start = datetime_to_epochtime(dt.datetime.strptime(date_start, "%Y-%m-%d"))
	date_end = (dt.datetime.strptime(date_end, "%Y-%m-%d") - dt.datetime(1969,12,31)).total_seconds()
	# 1969,12,31 because then the date_end day is INCLUDED

	file_nc = nc.Dataset(filename)

	mwr_dict = dict()
	for key in file_nc.variables.keys():
		if key == 'time':
			master_time = np.asarray(file_nc.variables[key])
		elif key in ['IWV', 'LWP']:
			data_mt = np.asarray(file_nc.variables[key])
		else:
			raise KeyError("Key %s was not found in file containing the concatenated IWV/LWP data on master time axis"%key)

	# Trim the data and time arrays:
	# case: time[0] < date_start: trim lower end
	# case: time[0] >= date_start: dont trim lower end
	# case: time[-1] > date_end: trim upper end
	# case: time[-1] <= date_end: dont trim upper end
	if master_time[0] < date_start:
		trim_low_idx = master_time >= date_start
		master_time = master_time[trim_low_idx]
		data_mt = data_mt[trim_low_idx]
	if master_time[-1] > date_end:
		trim_high_idx = master_time < date_end
		master_time = master_time[trim_high_idx]
		data_mt = data_mt[trim_high_idx]

	return master_time, data_mt


def import_concat_IWV_LWP_mwr_running_mean(
	filename,
	date_start,
	date_end,
	instrument):

	"""
	Simple importer to get the IWV or LWP data, to which a moving average over a certain time
	span was performed, of a radiometer (HATPRO, MiRAC-P, ARM). The loaded data will be trimmed
	according to the specified date_start and date_end.

	Parameters:
	-----------
	filename : str
		Filename and path of the IWV or LWP data of radiometers with moving average (running mean).
	date_start : str
		Marks the first day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	date_end : str
		Marks the last day of the desired period. To be specified in yyyy-mm-dd (e.g. 2021-01-14)!
	instrument : str
		Specifies which instrument is mwr_dict. Can only be 'hatpro', 'mirac' or 'arm'.
	"""

	if instrument not in ['hatpro', 'mirac', 'arm']:
		raise ValueError("'instrument' must be either 'hatpro', 'mirac' or 'arm'.")

	date_start = datetime_to_epochtime(dt.datetime.strptime(date_start, "%Y-%m-%d"))
	date_end = (dt.datetime.strptime(date_end, "%Y-%m-%d") - dt.datetime(1969,12,31)).total_seconds()
	# 1969,12,31 because then the date_end day is INCLUDED

	file_nc = nc.Dataset(filename)

	mwr_dict = dict()
	for key in file_nc.variables.keys():
		if key == 'rm_window':
			rm_window = int(np.asarray(file_nc.variables[key]))
		elif key == 'time':
			time_rm = np.asarray(file_nc.variables[key])
		elif key in ['IWV', 'LWP']:
			data_rm = np.asarray(file_nc.variables[key])

	# Trim the data and time arrays:
	# case: time_rm[0] < date_start: trim lower end
	# case: time_rm[0] >= date_start: dont trim lower end
	# case: time_rm[-1] > date_end: trim upper end
	# case: time_rm[-1] <= date_end: dont trim upper end
	if instrument in ['hatpro', 'arm', 'mirac']:
		if time_rm[0] < date_start:
			trim_low_idx = time_rm >= date_start
			time_rm = time_rm[trim_low_idx]
			data_rm = data_rm[trim_low_idx]
		if time_rm[-1] > date_end:
			trim_high_idx = time_rm < date_end
			time_rm = time_rm[trim_high_idx]
			data_rm = data_rm[trim_high_idx]
		

	return rm_window, data_rm, time_rm


def import_PS_mastertrack_tab(filename):

	"""
	Imports Polarstern master track data during MOSAiC published on PANGAEA. Time
	will be given in seconds since 1970-01-01 00:00:00 UTC and datetime. It also
	returns global attributes in the .tab file so that the information can be
	forwarded to the netcdf version of the master tracks.

	Leg 1, Version 2:
	Rex, Markus (2020): Links to master tracks in different resolutions of POLARSTERN
	cruise PS122/1, Tromsø - Arctic Ocean, 2019-09-20 - 2019-12-13 (Version 2). Alfred
	Wegener Institute, Helmholtz Centre for Polar and Marine Research, Bremerhaven, 
	PANGAEA, https://doi.org/10.1594/PANGAEA.924668

	Leg 2, Version 2:
	Haas, Christian (2020): Links to master tracks in different resolutions of POLARSTERN
	cruise PS122/2, Arctic Ocean - Arctic Ocean, 2019-12-13 - 2020-02-24 (Version 2).
	Alfred Wegener Institute, Helmholtz Centre for Polar and Marine Research,
	Bremerhaven, PANGAEA, https://doi.org/10.1594/PANGAEA.924674

	Leg 3, Version 2:
	Kanzow, Torsten (2020): Links to master tracks in different resolutions of POLARSTERN
	cruise PS122/3, Arctic Ocean - Longyearbyen, 2020-02-24 - 2020-06-04 (Version 2).
	Alfred Wegener Institute, Helmholtz Centre for Polar and Marine Research, 
	Bremerhaven, PANGAEA, https://doi.org/10.1594/PANGAEA.924681

	Leg 4:
	Rex, Markus (2021): Master tracks in different resolutions of POLARSTERN cruise
	PS122/4, Longyearbyen - Arctic Ocean, 2020-06-04 - 2020-08-12. Alfred Wegener 
	Institute, Helmholtz Centre for Polar and Marine Research, Bremerhaven, PANGAEA,
	https://doi.org/10.1594/PANGAEA.926829

	Leg 5:
	Rex, Markus (2021): Master tracks in different resolutions of POLARSTERN cruise
	PS122/5, Arctic Ocean - Bremerhaven, 2020-08-12 - 2020-10-12. Alfred Wegener
	Institute, Helmholtz Centre for Polar and Marine Research, Bremerhaven, PANGAEA,
	https://doi.org/10.1594/PANGAEA.926910

	Parameters:
	-----------
	filename : str
		Filename + path of the Polarstern Track data (.tab) downloaded from the DOI
		given above.
	"""

	n_prel = 20000		# just a preliminary assumption of the amount of data entries
	reftime = dt.datetime(1970,1,1)
	pstrack_dict = {'time_sec': np.full((n_prel,), np.nan),		# in seconds since 1970-01-01 00:00:00 UTC
					'time': np.full((n_prel,), reftime),		# datetime object
					'Latitude': np.full((n_prel,), np.nan),		# in deg N
					'Longitude': np.full((n_prel,), np.nan),	# in deg E
					'Speed': np.full((n_prel,), np.nan),		# in knots
					'Course': np.full((n_prel,), np.nan)}		# in deg

	f_handler = open(filename, 'r')
	list_of_lines = list()

	# identify header size and save global attributes:
	attribute_info = list()
	for k, line in enumerate(f_handler):
		attribute_info.append(line.strip().split("\t"))	# split by tabs
		if line.strip() == "*/":
			break
	attribute_info = attribute_info[1:-1]	# first and last entry are "*/"

	m = 0		# used as index to save the entries into pstrack_dict
	for k, line in enumerate(f_handler):
		if k > 0:		# skip header
			current_line = line.strip().split()		# split by tabs

			# convert time stamp to seconds since 1970-01-01 00:00:00 UTC:
			pstrack_dict['time_sec'][m] = datetime_to_epochtime(dt.datetime.strptime(current_line[0], "%Y-%m-%dT%H:%M"))

			# extract other info:
			pstrack_dict['Latitude'][m] = float(current_line[1])
			pstrack_dict['Longitude'][m] = float(current_line[2])
			pstrack_dict['Speed'][m] = float(current_line[3])
			pstrack_dict['Course'][m] = float(current_line[4])

			m = m + 1

	# truncate redundant lines:
	last_nonnan = np.where(~np.isnan(pstrack_dict['time_sec']))[0][-1] + 1		# + 1 because of python indexing
	for key in pstrack_dict.keys(): pstrack_dict[key] = pstrack_dict[key][:last_nonnan]

	# time to datetime:
	pstrack_dict['time'] = np.asarray([dt.datetime.utcfromtimestamp(tt) for tt in pstrack_dict['time_sec']])

	return pstrack_dict, attribute_info


def import_MOSAiC_Radiosondes_PS122_Level2_tab(filename):

	"""
	Imports level 2 radiosonde data launched from Polarstern
	during the MOSAiC campaign. Time will be given in seconds since 1970-01-01 00:00:00 UTC
	and datetime. Furthermore, the Integrated Water Vapour will be computed
	using the saturation water vapour pressure according to Hyland and Wexler 1983.

	Maturilli, Marion; Holdridge, Donna J; Dahlke, Sandro; Graeser, Jürgen;
	Sommerfeld, Anja; Jaiser, Ralf; Deckelmann, Holger; Schulz, Alexander 
	(2021): Initial radiosonde data from 2019-10 to 2020-09 during project 
	MOSAiC. Alfred Wegener Institute, Helmholtz Centre for Polar and Marine 
	Research, Bremerhaven, PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.928656 
	(DOI registration in progress)

	Parameters:
	-----------
	filename : str
		Filename + path of the Level 2 radiosonde data (.tab) downloaded from the DOI
		given above.
	"""

	n_sonde_prel = 3000		# just a preliminary assumption of the amount of radiosondes
	n_data_per_sonde = 12000	# assumption of max. time points per sonde
	reftime = dt.datetime(1970,1,1)
	# the radiosonde dict will be structured as follows:
	# rs_dict['0'] contains all data from the first radiosonde: rs_dict['0']['T'] contains temperature
	# rs_dict['1'] : second radiosonde, ...
	# this structure allows to have different time dimensions for each radiosonde
	rs_dict = dict()
	for k in range(n_sonde_prel):
		rs_dict[str(k)] = {'time': np.full((n_data_per_sonde,), reftime),		# datetime object
							'time_sec': np.full((n_data_per_sonde,), np.nan),	# in seconds since 1970-01-01 00:00:00 UTC
							'Latitude': np.full((n_data_per_sonde,), np.nan),	# in deg N
							'Longitude': np.full((n_data_per_sonde,), np.nan),	# in deg E
							'Altitude': np.full((n_data_per_sonde,), np.nan),	# in m
							'h_geom': np.full((n_data_per_sonde,), np.nan),		# geometric height in m
							'ETIM': np.full((n_data_per_sonde,), np.nan),		# elapsed time in seconds since sonde start
							'P': np.full((n_data_per_sonde,), np.nan),			# in hPa
							'T': np.full((n_data_per_sonde,), np.nan),			# in deg C
							'RH': np.full((n_data_per_sonde,), np.nan),			# in percent
							'wdir': np.full((n_data_per_sonde,), np.nan),		# in deg
							'wspeed': np.full((n_data_per_sonde,), np.nan)}		# in m s^-1


	f_handler = open(filename, 'r')

	# identify header size and save global attributes:
	attribute_info = list()
	for k, line in enumerate(f_handler):
		if line.strip().split("\t")[0] in ['Citation:', 'In:', 'Abstract:', 'Keyword(s):']:
			attribute_info.append(line.strip().split("\t"))	# split by tabs
		if line.strip() == "*/":
			break


	m = -1		# used as index to save the entries into rs_dict; will increase for each new radiosonde
	mm = 0		# runs though all time points of one radiosonde and is reset to 0 for each new radiosonde
	precursor_event = ''
	for k, line in enumerate(f_handler):
		if k > 0:		# skip header
			current_line = line.strip().split("\t")		# split by tabs
			current_event = current_line[0]			# marks the radiosonde launch

			if current_event != precursor_event:	# then a new sonde is found in the current_line
				m = m + 1
				mm = 0

			# convert time stamp to seconds since 1970-01-01 00:00:00 UTC:
			rs_dict[str(m)]['time'][mm] = dt.datetime.strptime(current_line[1], "%Y-%m-%dT%H:%M:%S")
			rs_dict[str(m)]['time_sec'][mm] = datetime_to_epochtime(rs_dict[str(m)]['time'][mm])

			# extract other info:
			try:
				rs_dict[str(m)]['Latitude'][mm] = float(current_line[2])
				rs_dict[str(m)]['Longitude'][mm] = float(current_line[3])
				rs_dict[str(m)]['Altitude'][mm] = float(current_line[4])
				rs_dict[str(m)]['h_geom'][mm] = float(current_line[5])
				rs_dict[str(m)]['ETIM'][mm] = float(current_line[6])
				rs_dict[str(m)]['P'][mm] = float(current_line[7])
				rs_dict[str(m)]['T'][mm] = float(current_line[8])
				rs_dict[str(m)]['RH'][mm] = float(current_line[9])
				rs_dict[str(m)]['wdir'][mm] = float(current_line[10])
				rs_dict[str(m)]['wspeed'][mm] = float(current_line[11])

			except ValueError:		# then at least one measurement is missing:
				for ix, cr in enumerate(current_line):
					if cr == '':
						current_line[ix] = 'nan'
				try:
					rs_dict[str(m)]['Latitude'][mm] = float(current_line[2])
					rs_dict[str(m)]['Longitude'][mm] = float(current_line[3])
					rs_dict[str(m)]['Altitude'][mm] = float(current_line[4])
					rs_dict[str(m)]['h_geom'][mm] = float(current_line[5])
					rs_dict[str(m)]['ETIM'][mm] = float(current_line[6])
					rs_dict[str(m)]['P'][mm] = float(current_line[7])
					rs_dict[str(m)]['T'][mm] = float(current_line[8])
					rs_dict[str(m)]['RH'][mm] = float(current_line[9])
					rs_dict[str(m)]['wdir'][mm] = float(current_line[10])
					rs_dict[str(m)]['wspeed'][mm] = float(current_line[11])

				except IndexError:		# GPS connection lost
					rs_dict[str(m)]['Latitude'][mm] = float('nan')
					rs_dict[str(m)]['Longitude'][mm] = float('nan')
					rs_dict[str(m)]['Altitude'][mm] = float(current_line[4])
					rs_dict[str(m)]['h_geom'][mm] = float('nan')
					rs_dict[str(m)]['ETIM'][mm] = float(current_line[6])
					rs_dict[str(m)]['P'][mm] = float(current_line[7])
					rs_dict[str(m)]['T'][mm] = float(current_line[8])
					rs_dict[str(m)]['RH'][mm] = float(current_line[9])
					rs_dict[str(m)]['wdir'][mm] = float('nan')
					rs_dict[str(m)]['wspeed'][mm] = float('nan')

			mm = mm + 1
			precursor_event = current_event

	# truncate redundantly initialised sondes:
	for k in range(m+1, n_sonde_prel): del rs_dict[str(k)]
	
	# finally truncate unneccessary time dimension for each sonde and compute IWV:
	for k in range(m+1):
		last_nonnan = np.where(~np.isnan(rs_dict[str(k)]['time_sec']))[0][-1] + 1		# + 1 because of python indexing
		for key in rs_dict[str(k)].keys(): rs_dict[str(k)][key] = rs_dict[str(k)][key][:last_nonnan]
		rs_dict[str(k)]['q'] = np.asarray([convert_rh_to_spechum(t+273.15, p*100, rh/100) 
								for t, p, rh in zip(rs_dict[str(k)]['T'], rs_dict[str(k)]['P'], rs_dict[str(k)]['RH'])])
		rs_dict[str(k)]['rho_v'] = np.asarray([convert_rh_to_abshum(t+273.15, rh/100) 
								for t, rh in zip(rs_dict[str(k)]['T'], rs_dict[str(k)]['RH'])])
		rs_dict[str(k)]['IWV'] = compute_IWV_q(rs_dict[str(k)]['q'], rs_dict[str(k)]['P']*100)
	
	return rs_dict, attribute_info


def import_PS_mastertrack(
	filename,
	keys='all'):

	"""
	Imports Polarstern master track data during MOSAiC published on PANGAEA. Time
	will be given in seconds since 1970-01-01 00:00:00 UTC and datetime. It also
	returns global attributes in the .tab file so that the information can be
	forwarded to the netcdf version of the master tracks.

	Parameters:
	-----------
	filename : str
		Filename + path of the Polarstern Track data (.nc).
	keys : list of str
		List of names of variables to be imported. 'all' will import all keys.
		Default: 'all'
	"""

	file_nc = nc.Dataset(filename)

	if keys == 'all':
		keys = file_nc.variables.keys()

	elif isinstance(keys, str) and (keys != 'all'):
		raise ValueError("Argument 'keys' must either be a string ('all') or a list of variable names.")

	ps_track_dict = dict()
	for key in keys:
		if not key in file_nc.variables.keys():
			raise KeyError("I have no memory of this key: '%s'. Key not found in file '%s'." %(key, filename))
		ps_track_dict[key] = np.asarray(file_nc.variables[key])

	return ps_track_dict


def import_synthetic_TBs(
	files):

	"""
	Imports simulated TBs (must be .nc) from a certain site (e.g. Ny Alesund radiosondes)
	saved in a format readable by the mwr_pro retrieval.
	Also converts time to sec since 1970-01-01 00:00:00 UTC.

	Parameters:
	-----------
	files : str or list of str
		Path in which all files are selected to be imported or list of files
			to be imported.
	"""

	if type(files) == str:
		files = sorted(glob.glob(files + "*.nc"))

	DS = xr.open_mfdataset(files, concat_dim='n_date', combine='nested',
								preprocess=syn_MWR_cut_useless_variables_TB)

	# convert time:
	time_dt = np.asarray([dt.datetime.strptime(str(tttt), "%Y%m%d%H") for tttt in DS.date.values])
	time = datetime_to_epochtime(time_dt)
	DS['time'] = xr.DataArray(time,	dims=['n_date'])

	# Cut unwanted dimensions in variables 'frequency' and 'elevation_angle':
	DS['frequency'] = DS.frequency[0,:]
	DS['elevation_angle'] = DS.elevation_angle[0]

	return DS


def import_mwr_pro_radiosondes(
	files,
	**kwargs):

	"""
	Imports radiosonde (or radiosonde-like) data (must be .nc) from a certain site 
	(e.g. Ny Alesund radiosondes) saved in a format readable by the mwr_pro retrieval. 
	Also converts time to sec since 1970-01-01 00:00:00 UTC and computes specific 
	and relative humidity from absolute humidity. 

	Parameters:
	-----------
	files : str or list of str
		Path in which all files are selected to be imported or list of files
			to be imported.

	**kwargs:
	with_lwp : bool
		If True, the sonde dict will also include liquid water path.
	"""

	if type(files) == str:
		files = sorted(glob.glob(files + "*.nc"))

	DS = xr.open_mfdataset(files, concat_dim='n_date', combine='nested',
								preprocess=syn_MWR_cut_useless_variables_RS)

	# convert time:
	time_dt = np.asarray([dt.datetime.strptime(str(tttt), "%Y%m%d%H") for tttt in DS.date.values])
	time = datetime_to_epochtime(time_dt)
	n_time = len(time)
	DS['time'] = xr.DataArray(time,	dims=['n_date'])

	# compute specific and relative humidity:
	spec_hum = convert_abshum_to_spechum(DS.atmosphere_temperature.values, DS.atmosphere_pressure.values, 
											DS.atmosphere_humidity.values)
	spec_hum_sfc = convert_abshum_to_spechum(DS.atmosphere_temperature_sfc.values, DS.atmosphere_pressure_sfc.values, 
											DS.atmosphere_humidity_sfc.values)
	DS['atmosphere_spec_humidity'] = xr.DataArray(spec_hum, dims=['n_date', 'n_height'])
	DS['atmosphere_spec_humidity_sfc'] = xr.DataArray(spec_hum_sfc, dims=['n_date'])

	rel_hum = convert_abshum_to_relhum(DS.atmosphere_temperature.values, DS.atmosphere_humidity.values)
	rel_hum_sfc = convert_abshum_to_relhum(DS.atmosphere_temperature_sfc.values, DS.atmosphere_humidity_sfc.values)
	DS['atmosphere_rel_humidity'] = xr.DataArray(rel_hum, dims=['n_date', 'n_height'])
	DS['atmosphere_rel_humidity_sfc'] = xr.DataArray(rel_hum_sfc, dims=['n_date'])


	# convert to dict:
	sonde_dict = {	'pres': DS.atmosphere_pressure.values,		# in Pa, time x height
					'temp': DS.atmosphere_temperature.values,	# in K
					'rh': DS.atmosphere_rel_humidity.values,	# between 0 and 1
					'height': DS.height_grid.values,			# in m; also time x height
					'rho_v': DS.atmosphere_humidity.values,		# in kg m^-3
					'q': DS.atmosphere_spec_humidity.values,	# in kg kg^-1
					'wspeed': np.zeros(DS.atmosphere_temperature.shape), 	# in m s^-1; is 0 because unknown
					'wdir': np.zeros(DS.atmosphere_temperature.shape), 	# in deg; is 0 because unknown
					'lat': np.repeat(DS.latitude, n_time),		# in deg N
					'lon': np.repeat(DS.longitude, n_time),		# in deg E
					'launch_time': time,						# in sec since 1970-01-01 00:00:00 UTC
					'iwv': DS.integrated_water_vapor.values}	# in kg m^-2

	if kwargs['with_lwp']:
		sonde_dict['lwp'] = DS.liquid_water_path.values			# in kg m^-2

	return sonde_dict
