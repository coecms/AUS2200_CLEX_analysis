#keys=wind_speed t2m td_2m merid_wind zonal_wind sfcp gales high_t2m low_wind mid_wind" 
import scipy.constants as scc

obs_field={'wind_speed':'obs_wind_speed',\
    't2m':'obs_screen_temperature',\
    'td_2m':'obs_dew_point',\
    'merid_wind':'obs_meridional_wind',\
    'zonal_wind':'obs_zonal_wind',\
    'sfcp':'obs_station_pressure',\
    'gales':'obs_wind_speed',\
    'high_t2m':'obs_screen_temperature',\
    'low_wind':'obs_wind_speed',\
    'mid_wind':'obs_wind_speed',\
    'sphm':'obs_specific_humidity',\
    'rh2m':'obs_relative_humidity'}

model_field={'wind_speed':'adj_wind_speed',\
    't2m':'adj_screen_temperature',\
    'td_2m':'adj_dew_point',\
    'merid_wind':'adj_meridional_wind',\
    'zonal_wind':'adj_zonal_wind',\
    'sfcp':'adjHY_station_pressure',\
    'gales':'adj_wind_speed',\
    'high_t2m':'adj_screen_temperature',\
    'low_wind':'adj_wind_speed',\
    'mid_wind':'adj_wind_speed',\
    'sphm':'adj_specific_humidity',\
    'rh2m':'bkg_relative_humidity'}

header={'wind_speed':'10m wind speed',\
    't2m':'2m temperature',\
    'td_2m':'2m dewpoint',\
    'merid_wind':'10m v-wind speed',\
    'zonal_wind':'10m u-wind speed',\
    'sfcp':'Sfc Pressure',\
    'gales':'10m wind speed > 35 kts',\
    'high_t2m':'2m temperature > 35kts',\
    'low_wind':'10m wind speed < 5kts',\
    'mid_wind':'10m wind speed (moderate)',\
    'sphm':'2m specific humidity',\
    'rh2m':'2m relative humidity'}

fc_offset={'wind_speed':None,\
    't2m':scc.zero_Celsius,\
    'td_2m':scc.zero_Celsius,\
    'merid_wind':None,\
    'zonal_wind':None,\
    'sfcp':None,\
    'gales':None,\
    'high_t2m':scc.zero_Celsius,\
    'low_wind':None,\
    'mid_wind':None,\
    'sphm':None,\
    'rh2m':None}
