import xarray as xr
import xgcm as xgcm
import numpy as np
'''
cloud analysis library: 
contact: Fadhlil R Muhammad, Corey Robinson, & Yi Huang
'''


def ctt_calc_icesnow(qis,qg,qc,qr,t,
                    latname='lat',lonname='lon',
                    vertname = 'model_theta_level_number',
                    thres = 0.00001,
                    celsius = True
                    ):

    ''' 
        def ctt_calc_icesnow: calculates the CTT height based on qice+snow,qgraupel,
        qcloud,qrain, and temperature. All must have the same dimension size and name. 
        
        input variables:
        qis: Q ICE+SNOW
        qg:  Q GRAUPEL
        qc:  Q CLOUD
        qr:  Q RAIN
        T :  TEMPERATURE
        
        latname  : name of latitude (.e.g. lat) ;  default = lat
        lonname  : name of longitude (.e.g. lon);  default = lon
        vertname : name of vertical coordinate  ;  defailt = model_theta_level_number
    
        thres    : threshold of the mixing ratio where CTT is defined; default = 0.00001 (0.01g/kg)
        celsius  : convert to celsius; default = True
      '''
    
    sum_hydro = qis + qg + qc + qr

    qsum = sum_hydro.rename('qsum')
    t    = t.rename('t')

    ds = xr.merge([qsum,t], compat='override')
    ds

    grid = xgcm.Grid(
        ds,
        coords={
            lonname: {'center': lonname},
            latname: {'center': latname},
            vertname: {'center': vertname}
        },
        periodic=False
    )

    thres_levels = np.array([thres])

    ctt = grid.transform(
        ds['t'],
        vertname,
        thres_levels,
        target_data=ds['qsum'],
        method='linear'
    )
    
    if celsius is True:
        ctt = ctt-273.

    return(ctt)

def ctt_calc(qi,qs,qg,qc,qr,t,                   #ice,snow,graupel,cloud,rain,temp
            latname='lat',lonname='lon',
            vertname = 'model_theta_level_number',
            thres = 0.00001,
            celsius = True
            ):
    ''' 
        def ctt_calc: calculates the CTT height based on qice,qsnow,qgraupel,
        qcloud,qrain, and temperature. All must have the same dimension size and name. 
        
        input variables:
        qi:  Q ICE
        qs:  Q SNOW
        qg:  Q GRAUPEL
        qc:  Q CLOUD
        qr:  Q RAIN
        T :  TEMPERATURE
        
        latname  : name of latitude (.e.g. lat) ;  default = lat
        lonname  : name of longitude (.e.g. lon);  default = lon
        vertname : name of vertical coordinate  ;  defailt = model_theta_level_number
    
        thres    : threshold of the mixing ratio where CTT is defined; default = 0.00001 (0.01g/kg)
        celsius  : convert to celsius; default = True
      '''
    sum_hydro = qi + qs + qg + qc + qr

    qsum = sum_hydro.rename('qsum')
    t    = t.rename('t')

    ds = xr.merge([qsum,t], compat='override')
    ds

    grid = xgcm.Grid(
        ds,
        coords={
            lonname: {'center': lonname},
            latname: {'center': latname},
            vertname: {'center': vertname}
        },
        periodic=False
    )

    thres_levels = np.array([thres])

    ctt = grid.transform(
        ds['t'],
        vertname,
        thres_levels,
        target_data=ds['qsum'],
        method='linear'
    )
    
    if celsius is True:
        ctt = ctt-273.

    return(ctt)