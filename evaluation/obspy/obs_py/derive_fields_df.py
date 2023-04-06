import numpy as np
import thermo

from .defns import field_dep, logger
from .ffddd_to_uv import ffddd_to_u, ffddd_to_v
from .uv_to_ffddd import uv_to_ddd, uv_to_ff
from .find_cvrt_keys import find_cvrt_keys

def derive_fields_df(df, pfx, gfile_info):
   """
   This is for 'simple' variable transformations - not height corrections etc.
      df: data frame
      pfx: one of 'obs_','bkg_','adj_'
   Will end up with versions of each variable listed in defns.field_dep

   Expect the different cases to be dealt with by different calls, as this 
      simplifies the logic a bit when dealing with different conversions
      (e.g. which dew_point or relative_humidty calculations to use - depends on which
            other variables are present)

   """

   # set up list of dataframe variables for this prefix type
   df_keys = []
   for k in df.keys():
     if k[:4] == pfx:
        df_keys.append(k)

   # determine which variables can be derived
   cvrt_info = {}
   var_list, cvrt_info, ncalls = find_cvrt_keys(df_keys, cvrt_info, field_dep)

   # cvrt_info is indexed by df variable to be created
   #   'l_indx' gives index of list of variables to use in transformation
   #   'var' gives name of basic variable to be generated (without prefix)
   # find_cvrt_keys also checks that ALL dependencies are available in df
   #    - as long as created in the order of var_list
   logger.debug(f"var_list: {var_list}")
   logger.debug(f"starting df keys: {df.keys()}")

   for v in var_list:

      if v in cvrt_info:
         cvrt_ok = False
         civ = cvrt_info[v] 
         logger.debug(f"deriving variable: {v}, {civ.keys()}")
         logger.debug(f"cvrt_info var: {civ['var']} l_indx: {civ['l_indx']}")

         if civ['var'] == 'dew_point':
            if civ['l_indx'] == 0:
               # td from t, rh 
               t = df[civ['input'][0]].values
               rh = df[civ['input'][1]].values
               vp = thermo.es(t)*rh/100.   # saturation vapour pressure
               df[v] = thermo.dewpt(vp)
               cvrt_ok = True

            elif civ['l_indx'] == 1:
               # td from p, spec_hum
               p = df[civ['input'][0]].values
               sphu = df[civ['input'][1]].values
               mr = thermo.mix_rto_from_specificHumidity(sphu)
               vp = thermo.water_vap_pres(p, mr)
               df[v] = thermo.dewpt(vp)
               cvrt_ok = True
   
         elif civ['var'] == 'relative_humidity':
            if civ['l_indx'] == 0:
               # rh from p, t & td
               p = df[civ['input'][0]].values
               t = df[civ['input'][1]].values
               td = df[civ['input'][2]].values
               df[v] = thermo.RH_from_dewpt(p, td, t)
               cvrt_ok = True

            elif civ['l_indx'] == 1:
               # rh from p, t, spec_hum 
               p = df[civ['input'][0]].values
               t = df[civ['input'][1]].values
               sphu = df[civ['input'][2]].values
               mr = thermo.mix_rto_from_specificHumidity(sphu)
               df[v] = thermo.relativeHumidity(p, t, mr)
               cvrt_ok = True
   
         elif civ['var'] == 'specific_humidity':
            if civ['l_indx'] == 0:
               # spechum from p, & td
               p = df[civ['input'][0]].values
               td = df[civ['input'][1]].values
               mr = thermo.mix_rto_from_dewpt(p,td)
               df[v] = thermo.specificHumidity(mr)
               cvrt_ok = True

            elif civ['l_indx'] == 1:
               # spechum from p, t & rh
               p = df[civ['input'][0]].values
               t = df[civ['input'][1]].values
               rh = df[civ['input'][2]].values
               mr = thermo.mix_rto_from_RH(p, t, rh)
               df[v] = thermo.specificHumidity(mr)
               cvrt_ok = True
   
         elif civ['var'] == 'zonal_wind':
            if civ['l_indx'] == 0:
               # zonal wind from speed & direction
               ws = df[civ['input'][0]].values
               wd = df[civ['input'][1]].values
               df[v] = ffddd_to_u(ws,wd)
               cvrt_ok = True
   
         elif civ['var'] == 'meridional_wind':
            if civ['l_indx'] == 0:
               # zonal wind from speed & direction
               ws = df[civ['input'][0]].values
               wd = df[civ['input'][1]].values
               df[v] = ffddd_to_v(ws,wd)
               cvrt_ok = True
   
         elif civ['var'] == 'wind_speed':
            if civ['l_indx'] == 0:
               # wind speed from u & v
               wu = df[civ['input'][0]].values
               wv = df[civ['input'][1]].values
               df[v] = np.sqrt((wu*wu) + (wv*wv))
               cvrt_ok = True
   
         elif civ['var'] == 'wind_direction':
            if civ['l_indx'] == 0:
               # wind direction from u & v
               wu = df[civ['input'][0]].values
               wv = df[civ['input'][1]].values
               df[v] = uv_to_ddd(wu,wv)
               cvrt_ok = True
   
         if not cvrt_ok:
            raise ValueError
            logger.error(f"Unknown transformation option {cvrt_info[v]['l_indx']} for {v}")

   return df
