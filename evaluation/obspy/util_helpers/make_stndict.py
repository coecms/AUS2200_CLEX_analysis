#####
#####
#####  DEPRECATED AND TO BE DELETED
#####      - just available for transition
#####
#####
from .get_config_file import get_config_file, check_file

ftp_dir  = 'anon/gen/fwo/'
stn_file_def = 'IDY02126.dat'

i_stnno  = 0   # index of station number column in stn_file
i_name   = 1   #  -""-    station name     ---- "" -----
i_avname = 2   #          aviation name 
i_lat    = 3   #          latitude
i_lon    = 4   #          longitude
i_height = 5   #          station height
i_wmoid  = 6   #          WMO number
i_state  = 7   #          state

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_line( stn_l ):
   """
   puts a line of info from the station dictionary
     and into a simple dictionary
   """
   data = {}

   # get rid of trailing linefeeds, '"'s and split by commas
   stn_l = stn_l[:-1].replace('"','').split(',')

   # station name - get rid of "AWS"
   sname = stn_l[i_name].replace('AWS','').split()

   # station name - get rid of "AIRPORT" & recombine
   sname2 = sname
   if len(sname) > 1:
      add_name =  None
      if sname[-1][:4] == 'AIRP':
         add_name = '_airport'
      elif sname[-1][:4] == 'AERO':
         add_name = '_aerodrome'
      elif sname[-1] == 'AW' :
         add_name = '_aws'
      if not add_name is None:
         sname2 = sname[:-1]
         sname2[-1] = sname2[-1]+add_name

   str = ''
   for s in sname2 :
       str = str + s

   sname = str
   sname = sname.replace('(','_')
   sname = sname.replace(')','')
   
   data['bom_id'] = stn_l[i_stnno]
   data['name']   = sname.lower().replace(' ','')
   data['wmo_id'] = stn_l[i_wmoid]
   data['wmo_idi'] = int(stn_l[i_wmoid])
   data['avn_id'] = stn_l[i_avname].lower()
   data['lat']    = float(stn_l[i_lat])
   data['lon']    = float(stn_l[i_lon])
   data['height'] = float(stn_l[i_height])
   data['state']  = stn_l[i_state].strip()
   return data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_station_dict( dkey, update=None  ):
   """
   sets up the station dictionary (and number of stations)
   dkey = dictionary key, one of ['name', 'bom_id', 'wmo']
           'name'   = keys are station names
           'bom_id' = keys are Bureau ids
           'wmo_id' = keys are WMO ids
           'wmo_idi' = keys are WMO ids (integers)
           'avn_id' = keys are Aviation ids
   """
   stn_file = stn_file_def

   if update:
      import os
   
      # get reference station dictionary from ftp server
      
      ftp_cmd_file = 'ftp_session'
      f1 = open(ftp_cmd_file,'w')
      f1.write('cd '+ftp_dir+"\n")
      f1.write('get '+stn_file+"\n")
      f1.write('bye\n')
      f1.close()
   
      os.system('ftp ftp.bom.gov.au < '+ftp_cmd_file)
      os.system('dos2unix '+stn_file)
   else:
      if not check_file(stn_file) :
           stn_file_list = get_config_file(stn_file)
           stn_file = stn_file_list[-1]
           print('Using station file ',stn_file)
      

   # Now read & process file

   f = open( stn_file,'r')
   stn_dict = { 'Aus':{} }
   nstn = int( f.readline() )

   while 1:
        stn_info = f.readline()

        if not stn_info:
           break

        stn      = parse_line( stn_info )
        state    = stn['state']

        if state not in stn_dict:
             stn_dict[state] = {}

        stn_dict[state][stn[dkey]]  = stn.copy()
        stn_dict['Aus'][stn[dkey]]  = stn.copy()

   f.close()

   return stn_dict,nstn
