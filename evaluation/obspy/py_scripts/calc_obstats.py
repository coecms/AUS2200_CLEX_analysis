"""
########################################################################
#
# Calculate statistics obs-forecast diffs from netcdf files (e.g. from interp_ob.py)
#    TIME averaged, STATION averaged and HOURLY averaged
#
# This can iterate over EITHER model basedatetimes (adjust_basedate=True)
#     OR forecast steps from a single basedatetime (adjust_basedate=False)
#
########################################################################
"""
######################################
# Get configuration
#####################################
import os,sys
sys.path=['', '/projects/access/apps/odbserver/0.16.2-ops-intel1903-ompi402/lib/python2.7/site-packages', '/g/data/access/varpy/build/varpy/2021.03.0_patched', '/g/data/dp9/pjs548/mosrs_utils/r5885_obspy', '/g/data/hh5/public/apps/miniconda3/envs/analysis3-22.10/lib/python39.zip', '/g/data/hh5/public/apps/miniconda3/envs/analysis3-22.10/lib/python3.9', '/g/data/hh5/public/apps/miniconda3/envs/analysis3-22.10/lib/python3.9/lib-dynload', '/g/data/hh5/public/apps/miniconda3/envs/analysis3-22.10/lib/python3.9/site-packages', '/g/data/hh5/public/apps/miniconda3/envs/analysis3-22.10/lib/python3.9/site-packages/cdat_info-8.2.1-py3.9.egg']
os.environ["PYTHONPATH"]="/g/data/dp9/pjs548/mosrs_utils/r5885_obspy/util_helpers"
from pathlib import Path
import obs_py as opy

# PKG_PATH = Path(__file__).parent.absolute()
PKG_PATH = "/g/data/tm70/dm5220/AUS2200_CLEX_analysis/evaluation/obspy/py_scripts"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def check_file (fname, debug=None):
       """
           check if file exists and can be read
       """
       if not debug is None:
         opy.logger.debug(f'fname  : {fname}')
         opy.logger.debug(f'isfile : {os.path.isfile(fname)}')
         opy.logger.debug(f'access : {os.access(fname,os.R_OK)}')

       return os.path.isfile(fname) and os.access(fname,os.R_OK)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_config_file ( f_default=None, f_run=None, debug=None ):
   """
      Returns the names for default configuration files and
         run specific files
   
      f_default = filename of default configuration file
      f_run     = filename for configuration of this particular run
     
      Returns at most 3 filenames in a list
         in the order they should be executed via execfile
         file_list[0]  f_default from script source code directory
         file_list[1]  f_default from current directory
         file_list[2]  f_run from current directory

      Usage:

         ini_f = get_config_file( 'file.def', 'file.cfg' )
         for f in ini_f:
           execfile(f)

   """

   file_list=[]
 
   f_def_src = os.path.join(PKG_PATH,f_default)
   if check_file(f_def_src, debug) :
      file_list.append(f_def_src)
   
   # check current directory for default settings
   #   will be assumed to override defaults
   
   if check_file( f_default, debug ):
      file_list.append(f_default )
   
   # check current directory for run-time specific settings
   
   if (f_run is not None) and check_file( f_run, debug ):
      file_list.append(f_run )
   
   if len(file_list) < 1:
      print(f'  WARNING >> No config files found ')
      print(f'     f_def_src = [ {f_def_src} ]')
      print(f'     f_default = [ {f_default} ]')
      print(f'     f_run     = [ {f_run} ]')
   
   return file_list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ini_f = get_config_file( 'calc_obstats.def', 'calc_obstats.cfg' )

for f in ini_f:
  print('Executing commands from file ',f)
  exec(compile(open(f, "rb").read(), f, 'exec'))

######################################
# Main declarations 
######################################

# import the modules we need

import numpy as np
import cdms2 as cdms
from datetime import datetime as dt
from datetime import timedelta as delt
import pandas as pd
import time
import json
from util_helpers import get_station_dict, get_state_region, GenEncoder

######################################
# Start main program 
######################################

# default settings

stats_hdr = save_input_info()
stats_hdr['Version'] = '1.0.0'

model_fields=[model_field]
for i in range(1,num_exp):
   model_fields.append(model_field+'_'+str(i))

stime = time.time()

stn_dict, nstn = get_station_dict( 'wmo_idi', upd_station_dict )
if not(state_name is None) and (domain is None) :
   domain  = get_state_region( state_name )
   stn_names = state_name
else:
   stn_names = 'Aus'

################# read file and variables ###############################

# Set up list of keys to read
#   keys to identify each obs
#   required field (model & obs)
#   lat, lon for generating station plots

read_keys = obs_sort_keys[:]
read_keys.extend([ 'lat', 'lon', model_fields[0], obs_field])
if len(aux_fields) > 0:
     read_keys.extend(aux_fields)

slat = domain[0]
nlat = domain[1]
wlon = domain[2]
elon = domain[3]

if lat_range is None:
   lat_range = [slat, nlat]
if lon_range is None:
   lon_range = [wlon, elon]

base_dts=[]
file_dts=[]
if adjust_basedate:
  print(dates[0],date_fmt)
  # date0 is the first basedate of experiment0
  f_date0  = dt.strptime( dates[0], date_fmt )
  end_date = dt.strptime( date_last, date_fmt )
else:
  # date0 is the time of the observations
  # date1 is the first basedate of experiment1
  #PJS fc_steplast needs to be defined in the config file
  #   only needs to be defined for experiment0 as
  #   assume same number of steps in experiment1 as experiment0
  f_date0  = delt( seconds=int(fc_steps[0])*60*60)
  end_date = delt( seconds=int(fc_steplast)*60*60)
  for n in range(num_exp):
    base_dts.append(dt.strptime( dates[n], date_fmt ))
    file_dts.append(dt.strptime( dates[n], date_fmt ))
  #file_dt1 = dt.strptime( date1, date_fmt ) #unneeded?
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*
#* Allows for parallelism if required
#*
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_data_par( args ):
   """
   Reads in data from files and performs some preliminary
     processing - dropping missing data and any subsetting
   """
   iproc   = args['iproc']
   sr_name = 'plot_obstats.read_data_par (' +str(iproc)+ ') >> '

   fnames = args['fnames']
   num_exp    = args['num_exp']
   #fname1 = args['fname1']
   start_time = args['start_time']
   end_time   = args['end_time']
   obs_fields=['' for i in range(num_exp)]
   mdl_fields=['' for i in range(num_exp)]
   obs_fields[0] = args['obs_field']
   mdl_fields[0] = args['mdl_fields'][0]
   for n in range(1,num_exp):
     obs_fields[n] = args['obs_field']+'_'+str(n)
     mdl_fields[n] = args['mdl_fields'][n]
   obs_offset = args['obs_offset']
   obs_scale  = args['obs_scale']
   fc_offsets = args['fc_offsets']
   #fc1_offset = args['fc1_offset']
   
   ret_val = { 'missing_data':False, 'iproc':iproc }

   in_file = opy.obfile()
   print(sr_name+' Opening '+fnames[0])
   try:
      in_file.open(fnames[0])#, source='obstats')
   except :
      print(sr_name,fnames[0],' missing, skipping this file')
      ret_val['missing_data'] = True
      return ret_val

   obs_df0, obs_units  = in_file.read( obs_sort_keys, \
                          list_keys=read_keys, \
                          time_range=[start_time,end_time],\
                          lat_range=lat_range, lon_range=lon_range, \
                          id_list=stn_id, debug=trace )

   ret_val['field_units'] = obs_units[obs_field]
   in_file.close()

   #Reset the obs_index to be 0-N
   obs_df0['obs_index']=np.arange(len(obs_df0))
   obs_df0.set_index('obs_index',drop=False,inplace=True)

   if num_exp > 1:
     for n in range(1,num_exp):
       in_file = opy.obfile()
       print(f'Opening {fnames[n]}')
       try:
         in_file.open(fnames[n]) #, source='obstats')
       except :
         # skip this date
         print(f'{fnames[n]} missing, skipping this file')
         ret_val['missing_data'] = True
         return ret_val
       
       obs_df1, obs_units  = in_file.read( obs_sort_keys, \
                             list_keys=read_keys, \
                             time_range=[start_time,end_time],\
                             lat_range=lat_range, lon_range=lon_range, \
                             id_list=stn_id, debug=trace )
       in_file.close()
       #Reset the obs_index to be 0-N
       obs_df1['obs_index']=np.arange(len(obs_df1))
       obs_df1.set_index('obs_index',drop=False,inplace=True)

       obs_dfx = obs_df0.join( obs_df1, how='outer', rsuffix='_'+str(n))
       obs_df0 = obs_dfx
   else:
      obs_dfx = obs_df0

   obs_df = obs_dfx

   if not(obs_offset is None):
     for n in range(num_exp):
       obs_df[obs_fields[n]] = obs_df[obs_fields[n]] + obs_offset
     #if two_expts:
     #   obs_df[obs_fields[1]] = obs_df[obs_fields[1]] + obs_offset
   
   if not(obs_scale is None):
     for n in range(num_exp):
       obs_df[obs_fields[n]] = obs_df[obs_fields[n]] * obs_scale
     #if two_expts:
     #   obs_df[obs_field1] = obs_df[obs_field1] * obs_scale
   
   for n in range(num_exp):
     if not(fc_offsets[n] is None):
       obs_df[mdl_fields[n]] = obs_df[mdl_fields[n]] + fc_offsets[n]
   #if not(fc1_offset is None):
   #  obs_df[mdl_field1] = obs_df[mdl_field1] + fc1_offset
   
   if not(df_subset is None):
      obs_df = df_subset( obs_df )


   obs_df.dropna(thresh=10*num_exp,inplace=True)
   #SJR create multi_index
   dt_arr=[dt.strptime('{:}{:06d}'.format(int(d),int(t)),'%Y%m%d%H%M%S') for d,t in zip(obs_df['date'],obs_df['time'])]
   obs_df['dt']=dt_arr
   obs_df['id']=obs_df['Station_Identifier']
   obs_df.set_index(obs_sort_keys,inplace=True)
   obs_df['lat']=obs_df['latitude']
   obs_df['lon']=obs_df['longitude']

   ret_val['obs_df'] = obs_df

   stn_list = obs_df.index.get_level_values('id')
   stn_count={}
   for s in stn_list:
      stn_count[s] = obs_df.xs(s,level='id')[obs_field].count()
   
   ret_val['stn_count'] = stn_count
   ret_val['stn_list']  = stn_list

   return ret_val

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f_date   = f_date0
add_date = delt( seconds=date_inc*60*60 )
obs_df   = None

fnames_dt=[ [] for i in range(num_exp) ]

#fname1_dt  = []
start_time = []
end_time   = []
print("Looping dates for",obs_field,"and",model_field)
#######################################################
# Set up filenames etc. for each call to read_data_par
#######################################################
opy.logger.info(f' using dates (start, date,increment) : {f_date} {end_date} {add_date}')
while f_date <= end_date:

   # set up time for actual observation
   tobs_start = f_date - delt( minutes=diff_mins )
   tobs_end = f_date + delt( minutes=diff_mins-1 )
   if adjust_basedate:
      # Keep forecast step constant & adjust basedate
      fname_fcstep = fc_steps[0]
      file_dt = f_date - delt(hours=int(fc_steps[0]))   
   else:
      # Keep basedatetime constant & adjust forecast step
      tobs_start = base_dts[0] + tobs_start
      tobs_end = base_dts[0] + tobs_end
      fname_fcstep = str(int(f_date.total_seconds()/60/60))
      file_dt = file_dts[0]

   start_time.append( int( dt.strftime(tobs_start,'%Y%m%d%H%M%S') ) )
   end_time.append( int( dt.strftime(tobs_end,'%Y%m%d%H%M%S') ) )

   # Read in data
   # exec statement required if using functions such as zfill or fmt in filenames
   try:
       exec('fname_dtfmt='+in_fnames[0])
   except:
       fname_dtfmt = in_fnames[0]

   fnames_dt[0].append( dt.strftime( file_dt, fname_dtfmt ) )

   if num_exp > 1:
     for n in range(1,num_exp):
       if adjust_basedate:
         # Keep forecast step constant & adjust basedate
         fname_fcstep = fc_steps[n]
         file_dt = f_date - delt(hours=int(fc_steps[n]))   
       else:
         # Keep basedatetime constant & adjust forecast step
         #   allow for possible different basetime in exp1
         #     tobs_start+delt(diff_mins) gives middle of window
         fctimes[n] = tobs_start + delt(minutes=diff_mins) - file_dts[n]
         date_last = str(int(fctimes[n].total_seconds()/60/60))
         file_dt = file_dts[n]

      # exec statement required if using functions such as zfill or fmt in filenames
       try:
           exec('fname_dtfmt='+in_fnames[n])
       except:
           fname_dtfmt = in_fnames[n]
       fnames_dt[n].append( dt.strftime( file_dt, fname_dtfmt ) )
   #else:
   #   fnames_dt.append( None )
   
   f_date = f_date + add_date


#######################################################
# Read and process data
#######################################################
iproc = list(range(len(fnames_dt[0])))
args  = {}
stn_count = {}
stn_list  = []

for ip,ts,te in zip( iproc,start_time,end_time ):
   
   ff=[x[ip] for x in fnames_dt]
   
   args['iproc']  = ip
   args['fnames'] = ff
   args['start_time'] = ts
   args['end_time']   = te
   args['obs_field']  = obs_field
   args['mdl_fields'] = model_fields
   args['obs_offset'] = obs_offset
   args['obs_scale']  = obs_scale
   args['num_exp']  = num_exp
   args['fc_offsets'] = fc_offsets

   ret_val = read_data_par( args )
   if ret_val['missing_data']:
       continue

   field_units = ret_val['field_units'] 

   # Accumulate data from experiment

   if isinstance( obs_df, pd.DataFrame) :
      obs_df = obs_df.append( ret_val['obs_df'] )
   else:
      obs_df = ret_val['obs_df']

   stn_proc = ret_val['stn_count']
   for sk, si in stn_proc.items() :
      if sk in stn_count :
          stn_count[sk] = stn_count[sk] + si
      else: 
          stn_count[sk] = si

   stn_proc = ret_val['stn_list']
   for sid in stn_proc :
      if not sid in stn_list:
          stn_list.append(sid)

opy.logger.info(f'Time to read data = {time.time()-stime}')
stime = time.time()
   
#############################################
# Construct master DataFrame
#############################################

opy.logger.info(f'Time to merge data = {time.time()-stime}')

stime = time.time()

# Delete stations with less than minrep obs
stn_drop = []
for s in stn_list:
   if stn_count[s] < min_rep:
      stn_drop.append(s)

if len(stn_drop) > 0:
    obs_df = obs_df.drop(stn_drop,level='id')

stn_list = obs_df.index.get_level_values('id')

opy.logger.info(f'No of stations processed = {len(stn_count)}')
opy.logger.info(f'Time to delete {len(stn_drop)} under-reporting stns = {time.time()-stime}')
stime = time.time()

if len(obs_df) == 0:
   raise RunTimeError('Empty data frame, ending now')

#Fix pressure units. Ensure in hPa
if obs_field == 'obs_station_pressure': 
    dobs=obs_df[obs_field]
    obs_df[obs_field] = np.where( np.logical_and(dobs>10000,dobs<1000000),dobs/100.,dobs)
    for n in range(num_exp):
      dobs=obs_df[model_fields[n]]
      obs_df[model_field[n]] = np.where( np.logical_and(dobs>10000,dobs<1000000),dobs/100.,dobs)


obs_df['diff_0'] = obs_df[model_fields[0]]-obs_df[obs_field]

opy.logger.info(f'Time to calculate model-obs = {time.time()-stime}')
stime = time.time()

# Make sure wind direction difference in range [-180, 180]
if obs_field == 'obs_wind_direction' :
    dobs = obs_df['diff_0']
    obs_df['diff_0'] = np.where( dobs < -180., dobs+360., dobs )
    obs_df['diff_0'] = np.where( dobs > 180. , dobs-360., dobs )

    opy.logger.info(f'Time to adjust obs wind direction = {time.time()-stime}')
stime = time.time()

obs_df['AbsErr_0'] = obs_df['diff_0'].abs()
obs_df['SqErr_0'] = obs_df['diff_0']*obs_df['diff_0']
#if two_expts:
if num_exp > 1:
  for n in range(1,num_exp):
    i=str(n)
    obs_df['diff_'+i] = obs_df[model_fields[n]]-obs_df[obs_field]
    # Make sure wind direction difference in range [-180, 180]
    if obs_field == 'obs_wind_direction' :
        dobs = obs_df['diff_'+i]
        obs_df['diff_'+i] = np.where( dobs < -180., dobs+360., dobs )
        obs_df['diff_'+i] = np.where( dobs > 180. , dobs-360., dobs )
    
    obs_df['AbsErr_'+i] = obs_df['diff_'+i].abs()
    obs_df['SqErr_'+i] = obs_df['diff_'+i]*obs_df['diff_'+i]

opy.logger.info(f'Time to adjust wind direction differences = {time.time()-stime}')
stime = time.time()

# Set up output file name
# allow for fname_out to be defined in config file as 
#  EITHER
#       simple string
#  OR
#       a procedure call to be evaluated 

try:
   exec("tmp_fname_out="+fname_out)
except:
   pass
else:
   fname_out = tmp_fname_out


if station_stats:
   #############################################
   # station by station time averaged stats
   #############################################

   # stats should be calculated where all data present - so drop Nan
   obs_data  = obs_df.dropna().groupby(level='id')
   opy.logger.info(f'Time to drop missing station data = {time.time()-stime}')
   stime = time.time()

   obs_stats = calc_stats( obs_data, num_exp, stn_dict[stn_names] )

   opy.logger.info('Time to calculate station-based stats = {time.time()-stime}')
   stime = time.time()

   # generate list of stats
   # Unsure if this is redundant or a dodgey way to sort obs_stats.keys
   # stats_list = list(obs_stats.keys())
   # stats_list.sort()
   obs_stats['Header'] = {}
   for hk, hi in stats_hdr.items():
      obs_stats['Header'][hk] = hi
   obs_stats['Aux Info'] = dict(lat=obs_data['lat'].mean().values, \
                                lon=obs_data['lon'].mean().values, \
                                units=field_units )

   stn_out = open( fname_out+'_stn.json','w')
   stn_out.write( json.dumps(obs_stats,cls=GenEncoder) )
   stn_out.close()

#############################################
# Domain averaged time series stats
#############################################

def nearest_hour( obs_dt ):
   half_hr = delt(seconds=30*60)
   one_hr  = delt(seconds=60*60)
   prev_hr = dt(obs_dt.year, obs_dt.month, obs_dt.day, obs_dt.hour)
   next_hr = prev_hr + one_hr
   if (obs_dt-prev_hr) < half_hr :
       return prev_hr
   else :
       return next_hr

if tseries_stats:
   # stats should be calculated where all data present - so drop Nan
   obs_data   = obs_df.dropna().groupby(nearest_hour,level='dt')
   opy.logger.info(f'Time to drop missing time-based data = {time.time()-stime}')
   stime = time.time()
   
   obs_stats  = calc_stats( obs_data, num_exp )

   opy.logger.info(f'Time to calculate time-based stats = {time.time()-stime}')
   stime = time.time()
   
   # Unsure if this is redundant or a dodgey way to sort obs_stats.keys
   # stats_list = list(obs_stats.keys())
   # stats_list.sort()
   obs_stats['Header'] = {}
   for hk, hi in stats_hdr.items():
      obs_stats['Header'][hk] = hi
   obs_stats['Aux Info'] = dict( units=field_units )
       
   date_out = open( fname_out+'_date.json','w')
   date_out.write( json.dumps(obs_stats,cls=GenEncoder) )
   date_out.close()


opy.logger.info(f' Time to close files = {time.time()-stime}')
stime = time.time()

#############################################
# Domain & Date averaged stats for each hour
#############################################

def obs_hr( obs_dt ):
   return '{:02d}'.format(obs_dt.hour)

if hourly_stats:
   # stats should be calculated where all data present - so drop Nan
   obs_data   = obs_df.dropna().groupby(obs_hr,level='dt')
   opy.logger.info(f'Time to drop missing hourly-based data = {time.time()-stime}')
   stime = time.time()
   
   obs_stats  = calc_stats( obs_data, num_exp )
   
   opy.logger.info(f'Time to calculate hourly-based stats = {time.time()-stime}')
   stime = time.time()
   # Unsure if this is redundant or a dodgey way to sort obs_stats.keys
   # stats_list = list(obs_stats.keys())
   # stats_list.sort()
   obs_stats['Header'] = {}
   for hk, hi in stats_hdr.items():
      obs_stats['Header'][hk] = hi
   obs_stats['Aux Info'] = dict( units=field_units )
   
   hour_out = open( fname_out+'_hour.json','w')
   hour_out.write( json.dumps(obs_stats,cls=GenEncoder) )
   hour_out.close()
