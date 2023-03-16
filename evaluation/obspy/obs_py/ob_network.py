# Manages meta data and various helper functions for 
#   processing various 3rd party obs networks

import os
import yaml
from zoneinfo import ZoneInfo
from pathlib import Path

from .defns import PKG_PATH, logger

# index of varno can be found at: https://code.metoffice.gov.uk/doc/ops/ops-2022.12.0/doc/varno.hh

varno_translate_def={'obs_screen_pressure' : 107,
                 'obs_mslp' : 108,
                 'obs_screen_temperature' : 39,
		 'obs_screen_dew_point' : 40,
		 'obs_screen_relative_humidity' : 58,
		 'obs_10m_wind_speed' : 112,
		 'obs_10m_wind_direction' : 111,
		 'obs_10m_windgust' : 261,
		 'obs_10m_zonal_wind' : 41,
		 'obs_10m_meridional_wind' : 42,
                 'obs_height' : 1,
                 'obs_temperature' : 2,
		 'obs_dew_point' : 59,
		 'obs_specific_humidity' : 7,
		 'obs_wind_speed' : 112,
		 'obs_wind_direction' : 111,
		 'obs_zonal_wind' : 3,
		 'obs_meridional_wind' : 4,
		 }

# obs types are in VarMod_ObsIO/Var_WriteObsSensNRL.f90		 
odb_obs_type_def = {'synop':[10100,10101,10102], 
                    'mobsyn':[10800,11600], 
                    'metar':[11100,11101,11102], 
                    'ship':[10200,10201,10202], 
                    'buoy':[10300,11700,10310,11702], 
                    'amdar':[30100],
                    'airep':[30200],
                    'tamdar':[30300],
                    'mode-s':[30500],
                    'temp':[50100, 50101],
                    'tempship':[50102],
                    'tempmob':[50103],
                    'pilot':[50200, 50201],
                    'pilotship':[50202],
                    'pilotmob':[50203],
                    'dropsonde':[50300],
                    'winpro':[50400],
                    'sondebufr':[50500,50501]
                    }
odb_obs_type_def['surface'] = odb_obs_type_def['synop']  \
         + odb_obs_type_def['mobsyn']  \
         + odb_obs_type_def['metar'] + odb_obs_type_def['ship'] \
         + odb_obs_type_def['buoy']

odb_obs_type_def['temp_all'] = odb_obs_type_def['temp']  \
         + odb_obs_type_def['tempship']  \
         + odb_obs_type_def['tempmob']

odb_obs_type_def['pilot_all'] = odb_obs_type_def['pilot']  \
         + odb_obs_type_def['pilotship']  \
         + odb_obs_type_def['pilotmob']

odb_obs_type_def['upper_all'] = odb_obs_type_def['temp_all']  \
         + odb_obs_type_def['pilot_all']  \
         + odb_obs_type_def['winpro']  \
         + odb_obs_type_def['sondebufr']

odb_obs_type_def['sonde_all'] = odb_obs_type_def['temp_all']  \
         + odb_obs_type_def['pilot_all']  \
         + odb_obs_type_def['sondebufr']

odb_obs_type_def['aircraft'] = odb_obs_type_def['amdar']  \
         + odb_obs_type_def['airep'] + odb_obs_type_def['tamdar'] \
         + odb_obs_type_def['mode-s']

odb_columns = {'surface': [ \
                   'andate', 'antime', 'ops_obstype', 'date', 'time', 'statid','ident', 'lat',
                   'lon', 'stalt', 'initial_obsvalue', 'obsvalue', 'varno', 'seqno',
                   'ops_report_flags', 'ops_datum_flags'], \
               'upper': [ \
                   'andate', 'antime', 'ops_obstype', 'date', 'time', 'statid', 'lat',
                   'lon', 'stalt', 'initial_obsvalue', 'obsvalue', 'varno', 'seqno',
                   'ops_report_flags', 'ops_datum_flags',
                   'vertco_type', 'vertco_reference_1', 'vertco_reference_2'], \
               'aircraft': [ \
                   'andate', 'antime', 'ops_obstype', 'date', 'time', 'statid','lat',
                   'lon', 'stalt', 'initial_obsvalue', 'obsvalue', 'varno', 'seqno',
                   'ops_report_flags', 'ops_datum_flags',
                   'vertco_type', 'vertco_reference_1', 'flight_phase'] \
              } 

odb_regionfilter_def='(lon>110) & (lon<160) & (lat>-50) & (lat<-5)'

class ob_network( object ):
    """
        Class for information about observation networks
    """
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # class constructor
    #   initialization of variables

    def __init__(self, name='odb_surface', network_info=None, yaml_user_cfg=None, 
                       varno_translate=varno_translate_def,
                       odb_obs_type=odb_obs_type_def,
                       odb_columns=odb_columns['surface'],
                       odb_regionfilter=odb_regionfilter_def):

        # sets up meta-data and helper functions for a network
        # name - type of network, 

        yaml_base_cfg = os.path.join(PKG_PATH,'ob_network.yml')
        self.yaml_base_cfg = yaml_base_cfg
        print('base ob_network definition from: ',yaml_base_cfg)
        with open(yaml_base_cfg,'r') as stream:
            network_yaml = yaml.safe_load(stream)

        if yaml_user_cfg is not None:
          yaml_user_cfg = Path(yaml_user_cfg).absolute()
          print('user ob_network definition from: ',yaml_user_cfg)
          with open(yaml_user_cfg,'r') as stream:
              network_yaml_user = yaml.safe_load(stream)
          for k,i in network_yaml_user.items():
              network_yaml[k] = i.copy()
        self.yaml_user_cfg = yaml_user_cfg

        if network_info is not None:
           for k,i in network_info.items():
              network_yaml['network_info'][k] = i.copy()

        network_info = network_yaml['network_info'].copy()
        network_info['conditions over the land'] = network_info['cotl']
        network_info['landscape south australia'] = network_info['lsa']
        network_info['ga aws'] = network_info['geoscience']

        try:
            self.info = network_info[name]
        except:
            logger.error('unknown network name: '+name)
            logger.error('valid name: '+str(network_info.keys()) )

        self.name = name
        self.set_varno_translate(varno_translate)
        self.set_odb_obs_type(odb_obs_type)
        self.set_odb_columns(odb_columns)
        self.set_odb_regionfilter(odb_regionfilter)

        trans_vars = {}
        if 'inv_dict' in self.info:
           if self.info['inv_dict'] is None: self.info['inv_dict'] = {}
           for k,i in self.info['obs_dict'].items():
              if not isinstance(i,dict) : continue
              if 'transform' in i:
                for kt in i['transform']:
                    self.info['inv_dict'][kt] = {'units': i['transform'][kt]['units'],
                                                 'var': k}
              elif 'units' in i:
                    self.info['inv_dict'][i['var']] = {'units': i['units'], 
                                                       'var': k}

        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def __dir__(self):
 
        list_dir = ['__init__', '__dir__' \
                    '__odb_columns', '__odb_obs_type', \
                    '__odb_regionfilter', \
                    '__varno_translate', \
                    'odb_columns', 'odb_obs_type', \
                    'odb_regionfilter', \
                    'set_info', \
                    'set_odb_columns', 'set_odb_obs_type', \
                    'set_odb_regionfilter', \
                    'set_varno_translate',  \
                    'varno_translate',  \
                   ]
        return list_dir
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_info( nw_dict ):

        for k,i in nw_dict.items():
            self.info[k] = i
           
        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def varno_translate(self):

        return self.__varno_translate

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_varno_translate(self, val):

        self.__varno_translate = val
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def odb_obs_type(self):

        return self.__odb_obs_type

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def odb_columns(self):

        return self.__odb_columns

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def odb_regionfilter(self):

        return self.__odb_regionfilter

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_odb_obs_type(self, val):

        self.__odb_obs_type = val
        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_odb_columns(self, val):

        self.__odb_columns = val
        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_odb_regionfilter(self, val):

        self.__odb_regionfilter = val
        return

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

