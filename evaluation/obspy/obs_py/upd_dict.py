import unittest

def upd_dict(main_d, upd_item):
   # update an item in a dictionary, preserving all existing keywords
   # original item may or may not exist
   #     and may or may not be a dictionary

   for k,i in upd_item.items():
      if k not in main_d:
         # adding new item to dictionary
         main_d[k] = i
      elif isinstance(main_d[k],dict):
         # updating a dictionary
         # cannot just replace main_d[k] with i
         upd_dict(main_d[k], i)
      else:
         main_d[k] = i
   return main_d

class test_upd_dict(unittest.TestCase):
    def test_rtn(self):

       main_d = {'time_type': 'timestamp', 'vtime_var': 'time',
               'surface_height': {'ndim': 2, 'fname': 'topog.nc', 'fvar_name': 'topog', 'fc_type': 'an', 'lev_type': 'sfc'},
               'lsm': {'ndim': 2, 'fname': 'lnd_mask.nc', 'fvar_name': 'lnd_mask', 'fc_type': 'an', 'lev_type': 'sfc'},
               'screen_temperature': {'fname': 'temp_scrn.nc', 'fvar_name': 'temp_scrn'},
               '10m_zonal_wind': {'fname':'uwnd10m.nc', 'fvar_name': 'uwnd10m'},
               '10m_meridional_wind': {'fname':'vwnd10m.nc', 'fvar_name': 'vwnd10m'},
               '10m_windgust': {'fname': 'wndgust10m.nc', 'fvar_name': 'wndgust10m'}
               }

       upd =  { '10m_zonal_wind': {'fvar_name': 'u10'},
               '10m_meridional_wind': {'fvar_name': 'v10'} }

       new_d =  upd_dict(main_d, upd)
       ans_u10 = { 'fname':'uwnd10m.nc', 'fvar_name': 'u10'}
       ans_v10 = { 'fname':'vwnd10m.nc', 'fvar_name': 'v10'}

       self.assertEqual( new_d['10m_zonal_wind'], ans_u10 )
       self.assertEqual( new_d['10m_meridional_wind'], ans_v10 )

if __name__ == '__main__':
   unittest.main()
