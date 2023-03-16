import numpy as np

from .local_datetime_utc import local_datetime_utc
from .local_datetime_utc_date import local_datetime_utc_date
from .local_datetime_utc_time import local_datetime_utc_time
from .encode_id import encode_id
from thermo import mix_rto_from_dewpt, specificHumidity

# For conversion routines they will all be called with the same number of arguments
#    - so just add some dummy arguments d1, d2, ....
#       and make them optional for backward compatibility

# A dummy return to do nothing
def return_value(x, d1=None):
    return x

def recast_str8(x, d1=None):
    if not isinstance(x,str):
       x = np.array([id.strip() for id in x])
    return x.astype('<U8')

def spechum_from_dewpt(p, td):
   mix_rto = mix_rto_from_dewpt(p,td)
   return specificHumidity(mix_rto)

convert_fn = { 'local_datetime_utc': local_datetime_utc, 
               'local_datetime_utc_date': local_datetime_utc_date, 
               'local_datetime_utc_time': local_datetime_utc_time, 
               'spechum_from_dewpt': spechum_from_dewpt, 
               'return_value': return_value, 
               'encode_id': encode_id,
               '<U8': recast_str8,
               'np.float16': np.float16, 
               'np.float32': np.float32, 
               'np.float64': np.float64,
               'np.int16': np.int16, 
               'np.int32': np.int32, 
               'np.int64': np.int64 }

def ob_convert(fn_name, *args, **kwargs):
    """
    Link between entries in ob_network.yml that are not in 
    obfile_info.yml and functions to convert to obfile_info
    variables

    When Python 3.10 is common use
     match fn_name:
     case req_f:
           return req_f(kwargs)
    """
    return convert_fn[fn_name](*args, **kwargs)


