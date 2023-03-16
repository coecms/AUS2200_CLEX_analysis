import pandas as pd
import numpy as np
from datetime import datetime as dt

from .nc_df import nc_df

# This is a deprecated class to provide an interface to old
#   stats routines
# once we make the transition/interface to METplus should not
#   be needed

class obfile_deprecated( object ):

    def __init__(self):
       self.dt_fmt = '%Y%m%d%H%M%S'
       self.date_scale = 100*100*100
       return

    def read( self, obs_sort_keys=None, list_keys=None,
              time_range=None, lon_range=None, lat_range=None, \
              id_list=None, debug=None ):
       """
       sort_keys etc. are redundant: nc_df handles all that

       
       """
       df, self.file_info, nw_info  = nc_df(self.fname)
       self.nw_info = nw_info.info
       self.set_dt_fmt( self.file_info.var_dict()['date']['attributes']['units']+ \
                        self.file_info.var_dict()['time']['attributes']['units'] )
       self.set_date_scale( self.file_info.var_dict()['time']['attributes']['units'] )

       if list_keys is not None:
           pfx = ['adj','obs','bkg']
           for k in df.keys():
              if (k[:3] in pfx) and (k not in list_keys):
                  df = df.drop(columns=[k])
    
       odt = [ (d*self.date_scale)+t for d,t in zip(df['date'],df['time']) ]

       if obs_sort_keys is not None:
           if 'id' in obs_sort_keys:
                df['id'] = df['Station_Identifier']
           if 'dt' in obs_sort_keys:
                ddt = np.array([dt.strptime(str(t),self.dt_fmt) for t in odt])
                df['dt'] = pd.Series(ddt, df.index )
   
       if time_range is not None:
          mask1 = np.greater_equal(odt, time_range[0])
          mask2 = np.less_equal(odt, time_range[1])
          df = df[ np.logical_and(mask1, mask2) ]
       if lon_range is not None:
          df = df[ np.logical_and(df['longitude'] >= lon_range[0],  \
                                  df['longitude'] <= lon_range[1]) ]

       if lat_range is not None:
          df = df[ np.logical_and(df['latitude'] >= lat_range[0],  \
                                  df['latitude'] <= lat_range[1]) ]
    
       if (id_list is not None) and (len(id_list) > 0):
          ids = [s_id in id_list for s_id in df['Station_Identifier'] ]
          df = df[ ids ]
 
       units = {}
       for k,i in self.file_info.var_dict().items() :
          if isinstance(i,dict):
             if 'attributes' in i:
                 units[k] = i['attributes']['units']

       return df, units           
  
    def open(self, fname):
       """
       Just return filename (fname), as open
          will be done in nc_df later
       """
       self.fname = fname
       return
  
    def close(self):
       """
       Just return nothing, as close is nc_df 
       """
       return

    def set_dt_fmt(self, dt_fmt):
       self.dt_fmt = dt_fmt
       return

    def set_date_scale(self, date_scale):
       if isinstance(date_scale,str):
           # convert from time format (assume %H%M%S or subset)
           self.date_scale = np.power(10,len(date_scale))
       else:
           self.date_scale = date_scale
       return
