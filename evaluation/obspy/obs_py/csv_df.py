from datetime import datetime as dt
from datetime import timedelta
import pandas as pd
import numpy as np
import csv
from thermo import es, dewpt, C2K

from util_helpers import miss_val
from .defns import dtf_obs, logger
from .encode_id import encode_id
from .local_datetime_utc import local_datetime_utc
from .local_datetime_utc_date import local_datetime_utc_date
from .local_datetime_utc_time import local_datetime_utc_time
from .ffddd_to_uv import ffddd_to_u, ffddd_to_v
from .uv_to_ffddd import uv_to_ddd, uv_to_ff

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def csv_df( infile, stn_dict, nw, ob_info, miss_vals=miss_val() ) :
    """
      Read in CSV file and generate pandas dataframe
      infile    : input CSV file
      stn_dict  : dictionary of station info from obs_csv_stn_info
      nw        : instance of ob_network class
      ob_info   : instance of obs_var_info class
    """

    obs_d = nw.info['obs_dict']
    index_col = nw.info['csv_axis']
    df = pd.read_csv(infile, index_col=index_col)

    # rename rest of columns to match obs_py variables
    #   - needs to be done before start adding new columns
    new_col = {}
    col_list = list(df.keys())
    for c in col_list:
       if (c in obs_d) and (obs_d[c]['var'] != c):
           new_col[c] = obs_d[c]['var']
           logger.info('renaming column from {0} to {1}\n'.format(c,new_col[c]))
           print('renaming column from {0} to {1}'.format(c,new_col[c]))

    if len(new_col.keys()) > 0:
        df.rename(columns=new_col,inplace=True)

    # prune unwanted columns
    if 'nc_skip' in nw.info:
      for v in nw.info['nc_skip']:
        if v in df: del( df[v] )

    ncol = len(df.columns)

    # generate BoM encoded station id
    id_axis = nw.info['csv_axis'][1]
    id_type = ob_info.var_dict()['station_identifier']['type']
    # should only be one item in the following
    for kk, ii in nw.info['obs_dict'][id_axis]['transform'].items():
       stn_id = [id for id in df.index.get_level_values(level=id_axis)]
       stn_id = [ encode_id(ii, id) for id in np.array(stn_id,dtype=id_type) ]
       df.insert(ncol,kk,stn_id)
    ncol += 1

    # add netcdf axis variables
    rec_indx = np.arange(len(stn_id))
    df.insert(ncol, ob_info.nc_axis(), rec_indx)
    ncol += 1

    # fields to extracted from station dictionary rather than reported data
    #    e.g. lat, lon, instrument height, ...
    stn0 = list(stn_dict.keys())[0]
    for k in nw.info['station_data']:
       if k not in df:
           kv = nw.info['inv_dict'][k]['var']
           stn_data = [ stn_dict[ss][kv] for ss in stn_id ]
           df.insert(ncol,k,stn_data)
           ncol += 1
      
    # observed variables that need transformations
    #    assume there is another key in nw_info that provides
    #       any extra info
    # obs_dict is indexed by names in CSV file
    #    so use inv_dict to get the CSV version of 'datetime'
    #    format etc. is ob_network dependent, so need to go 
    #       to actual entry describing datetime in the CSV file
    dt_axis = nw.info['csv_axis'][0]

    # If CSV file has datetime variable
    #    assume transformed into date and time

    if 'datetime' in nw.info['inv_dict']:
        dt_name = nw.info['inv_dict']['datetime']['var']
    else:
        dt_name = nw.info['inv_dict']['date']['var']

    dt_vals = [dt.strptime(dd,nw.info['obs_dict'][dt_name]['format'])
                   for dd in df.index.get_level_values(level=dt_axis) ]
    for k_od,i_od in obs_d.items():
       if k_od in df: continue
       if 'transform' in i_od:
         for kv, iv in i_od['transform'].items():
            if iv['func'] == 'local_datetime_utc_date':
               otz = [ stn_dict[ss]['tz'] for ss in stn_id ]
               stn_data = local_datetime_utc_date(dt_vals, otz, iv['units'])
            elif iv['func'] == 'local_datetime_utc_time':
               otz = [ stn_dict[ss]['tz'] for ss in stn_id ]
               stn_data = local_datetime_utc_date(dt_vals, otz, iv['units'])
            elif iv['func'] == 'add_vert_temp_diff':
               stn_data = df[i_od['var']] + df[ iv['ref'] ]
            elif iv['func'] == 'ffddd_to_u':
               stn_data = ffddd_to_u(df[iv['speed']].values,
                                     df[iv['dir']].values)
            elif iv['func'] == 'ffddd_to_v':
               stn_data = ffddd_to_v(df[iv['speed']].values,
                                     df[iv['dir']].values)
            elif iv['func'] == 'uv_to_ff':
               stn_data = uv_to_ff(df[iv['u']].values,
                                    df[iv['v']].values)
            elif iv['func'] == 'uv_to_ddd':
               stn_data = uv_to_ddd(df[iv['u']].values,
                                    df[iv['v']].values)
            elif iv['func'] == 'dewpt_from_rh':
               vp = es(df[iv['t']].values+C2K)
               stn_data = dewpt(vp*df[iv['rh']].values/100.)
            elif iv['func'] == 'encode_id':
               continue # already done above
            else:
               logger.error('Unknown transform name: {0}\n'.format(iv['func']))
               raise ValueError
            df.insert(ncol,kv,stn_data)
            ncol += 1
 
    for k_od,i_od in obs_d.items():
       if (i_od['var'] in df) and ('transform' in i_od):
          if i_od['var'] not in i_od['transform']:
             del df[i_od['var']]
             ncol -= 1

    # replace NaN with missing value

    for k,i in df.items():
       if isinstance(i[0],int):
          df[k].fillna(miss_vals.imiss(),inplace=True)
       else:
          # leave floast as np.nan
          continue

 
    return df
