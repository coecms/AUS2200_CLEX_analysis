from zoneinfo import ZoneInfo
from datetime import datetime as dt
from timezonefinder import TimezoneFinder
import csv

from util_helpers import miss_val

from .encode_id import encode_id
from .ob_convert import ob_convert
from .defns import logger

tzf = TimezoneFinder()
xmiss = miss_val().xmiss()
imiss = miss_val().imiss()

def obs_csv_stn_info(stn_file, nw, ob_info):
    """
    Set up dictionary of station meta-data, read from CSV file
       stn_file: csv file of station data
       nw: instance of ob_network class
       ob_info: instance of obfile_info class var_dict
    """

    stn_info = {}
    obs_d = nw.info['obs_dict']

    # find csv names for stn_index and station_name
    for k,i in obs_d.items():
       if 'var' in i:
          iv = i['var']
       if iv == 'station_identifier' : csv_stn_index = k
       if iv == 'local_id' :
            csv_stn_index = k
            id_key = list(i['transform'].keys())[0]
            encode_info = i['transform'][id_key]
       if iv == 'station_name' : csv_stn_name = k

    # read CSV file and unpack into dictionary
    with open(stn_file) as csvfile:
        for row in csv.DictReader(csvfile):
           stn_id = encode_id( encode_info, row[csv_stn_index] )
           stn_name = row[csv_stn_name]
           stn_info[stn_id] = {}

           for k,val in row.items():
               if (k in obs_d) and (obs_d[k]['var'] in ob_info) :
                 otype = ob_info[obs_d[k]['var']]['type']

                 if otype[:3] == 'np.' :
                   stn_info[stn_id][k] = ob_convert(otype, val)
                 else:
                   stn_info[stn_id][k] = val
   
           for k,i in nw.info['station_data'].items():
               if (k not in stn_info) and (i is not None) :
                   stn_info[stn_id][k] = i['val']

           stn_info[stn_name] = stn_info[stn_id]

    stn_list = stn_info.keys()
    for stn_id in stn_list:
        if nw.info['tzone'] == 'find':
            stn_info[stn_id]['tz'] = tzf.timezone_at(lng=stn_info[stn_id]['lon'], 
                                                    lat=stn_info[stn_id]['lat'] )
        else:
            stn_info[stn_id]['tz'] = ZoneInfo(nw.info['tzone'])
        stn_info[stn_name]['tz'] = stn_info[stn_id]['tz']

                                         
    return stn_info
