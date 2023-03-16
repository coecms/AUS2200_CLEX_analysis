import numpy as np
from .defns import ht_corr_dep, logger
from .adjust_T_by_dz import adjust_T_by_dz
from .adjust_w10_by_dz import adjust_w10_by_dz
from .adjust_RH_by_dz import adjust_RH_by_dz
from .adjust_p_HY import adjust_p_HY

def height_correct_df(df, pfx):

   # set up list of dataframe variables for this prefix type
   df_keys = []
   for k in df.keys():
     if k[:4] == pfx:
        df_keys.append(k)

     if pfx == 'bkg_':
         if ('bkg_surface_height' in df.keys()) \
              and ('obs_surface_height' in df.keys())  :

            dz = df['obs_surface_height'].values - df['bkg_surface_height'].values 

            if ( ('obs_screen_temperature' in df) and \
                 ('bkg_screen_pressure' in df) and \
                 ('obs_screen_dew_point' in df) ):
               t  = df['obs_screen_temperature'].values
               p  = df['bkg_screen_pressure'].values
               td = df['obs_screen_dew_point'].values
               df['adj_screen_pressure'] = adjust_p_HY(t, p, dz, td)

            if 'bkg_screen_temperature' in df:
               t = df['bkg_screen_temperature'].values
               df['adj_screen_temperature'] = adjust_T_by_dz(t, dz)

            if 'bkg_10m_zonal_wind' in df:
               u10 = df['bkg_10m_zonal_wind'].values
               df['adj_10m_zonal_wind'] = adjust_w10_by_dz(u10, v10=None, dz=dz)

            if 'bkg_10m_meridional_wind' in df:
               u10 = df['bkg_10m_meridional_wind'].values
               df['adj_10m_meridional_wind'] = adjust_w10_by_dz(u10, v10=None, dz=dz)

            if 'bkg_10m_windgust' in df:
               u10 = df['bkg_10m_windgust'].values
               df['adj_10m_windgust'] = adjust_w10_by_dz(u10, v10=None, dz=dz)

            if 'bkg_screen_relative_humidity' in df:
               rh = df['bkg_screen_relative_humidity'].values
               df['adj_screen_relative_humidity'] = adjust_RH_by_dz(rh, dz)
               # Should be safe to assume all of relhum, spec.hum and mix.ratio all defined by
               #    call to derive_fields_df
               # dewpoint & specific humidty will be caught by later call to derive_fields

     use_10m = False
     for k in df.keys():
         try:
            use_10m = use_10m or ( k.index('_10m_') > 0 )
         except:
            pass
     logger.debug(f"use_10m: {use_10m}")
     logger.debug(f"df.keys: {df.keys()}")
     logger.debug(f"z_lev: {'z_lev' in df}")
     logger.debug(f"z_lev_u: {'z_lev_u' in df}")
     
     if (use_10m):
        if ( (pfx+'screen_temperature' in df) and \
             (pfx+'screen_pressure' in df) and \
             (pfx+'screen_dew_point' in df) ):
           dz = np.array( [10.]*len(df) )
           df[pfx+'10m_pressure'] = adjust_p_HY(df[pfx+'screen_temperature'].values,
                                                    df[pfx+'screen_pressure'].values,dz, 
                                                    df[pfx+'screen_dew_point'].values)
        if ('z_lev' in df) and ('z_lev_u' in df):
           if pfx+'screen_temperature' in df:
              dz = df['z_lev_u'].values - df['z_lev'].values
              df[pfx+'10m_temperature'] = adjust_T_by_dz(df[pfx+'screen_temperature'].values, dz)

           if ('z_lev_u' in df) and \
              ('bkg_roughness_length' in df) :
                 ruff = df['bkg_roughness_length'].values
                 lev2m = df['z_lev'].values
                 lev10m = df['z_lev_u'].values
                 wscale = np.log(lev2m/ruff) / np.log(lev10m/ruff)

                 # screen winds from 10m are bkg
                 if (pfx+'screen_zonal_wind' not in df) and (pfx+'10m_zonal_wind' in df):
                    df[pfx+'screen_zonal_wind'] = df[pfx+'10m_zonal_wind'].values * wscale

                 if (pfx+'screen_meridional_wind' not in df) and (pfx+'10m_meridional_wind' in df):
                    df[pfx+'screen_meridional_wind'] = df[pfx+'10m_meridional_wind'].values * wscale

                 if (pfx+'screen_wind_speed' not in df) and (pfx+'10m_wind_speed' in df):
                    df[pfx+'screen_wind_speed'] = df[pfx+'10m_wind_speed'].values * wscale
       
                 # 10m winds from screen winds are bias corrected
                 if pfx[:4] == 'obs_':
                    pfxb = 'bco'+pfx[3:]
                    if (pfxb+'screen_zonal_wind' not in df) and (pfx+'screen_zonal_wind' in df):
                       df[pfxb+'screen_zonal_wind'] = df[pfx+'screen_zonal_wind'].values / wscale

                    if (pfxb+'screen_meridional_wind' not in df) and (pfx+'screen_meridional_wind' in df):
                       df[pfxb+'screen_meridional_wind'] = df[pfx+'screen_meridional_wind'].values / wscale

                    if (pfxb+'screen_wind_speed' not in df) and (pfx+'screen_wind_speed' in df):
                       df[pfxb+'screen_wind_speed'] = df[pfx+'screen_wind_speed'].values / wscale

     return df
