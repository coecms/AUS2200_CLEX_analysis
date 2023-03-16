from datetime import datetime as dt
from datetime import timedelta as delt
from glob import glob
import os
from pathlib import Path
import numpy as np
import yaml
import cdms2 as cdms

from .defns import dtf_grid, logger
from .get_nc_time import get_nc_time
from .upd_dict import upd_dict

PKG_PATH = Path(__file__).parent.absolute()

var_delim_def = ['{{','}}'] # Default delimiters for internal variables in file names

class archive(object):
   """
    A class for holding information about archived grid file formats
    - there are a number of interdepencies. This allows setting of
      values and maintaining the interdependcies

    Contains methods for checking base datetimes, valid datetimes and other meta data

    Different forecast validity times & units for time (hour, minutes etc.) handled by
      time_type = date_time  use valid_date & valid_time
                        OR timestamp  POSIX std. of time since-01-JAN-1970
                        OR time       time since basetime

   """

   def __init__(self, dir_type='mars_ncdf.0', 
                       root=None, dtfmt=None,
                       da_win=None, fcst_freq=None,
                       fc_type=None, lev_type=None,
                       cfg='archive.yml', info_items=None):
      """
      Class Constructor
         __XX are the variable values
         XX are functions to return the values

      dir_type: style for generating directory names
          mars_ncdf.0 = Pre 20200721 NCI archive
          mars_ncdf.1 = Post 20200721 NCI archive (no valid_date, base_time etc.)
          iris_ncdf.0 = output from IRIS

      root: top of archive directory - prefix for filenames

      dtfmt: dt_fmt class instance

      da_win: assimilation window (minutes)

      fcst_freq: forecast frequency (minutes)
                 (first versions of ACCESS_C3 had da_win=60, fcst_freq=360, 
                   i.e, forecasts every 6hours, but hourly DA)

      cfg : user yaml config file (default is archive.yml in working directory)

      info_items: run time modifications or mods too small to bother about
                    putting into config file
      """

      yaml_base_cfg = os.path.join(PKG_PATH,'archive.yml')
      print('Archive meta-data read from ',yaml_base_cfg)
      with open(yaml_base_cfg,'r') as stream:
         arch_yml = yaml.safe_load(stream)
    
      if os.access(cfg,os.R_OK):
        with open(cfg,'r') as stream:
           user_yml = yaml.safe_load(stream)
        print('User defined archive meta-data read from ',cfg)
        arch_yml = upd_dict( arch_yml, user_yml)

      if root is None:
        root=os.getenv('NC_ARCH')
        # "if root" will evaluate false for None, '', False etc.
        if not root: root = './' 

      # None ==> use default value
      #  (if have these values as default settings, then get
      #     reset to None when ask for default)
      if da_win is None and fcst_freq is None:
         da_win = 60
         fsct_freq = 360
      if da_win is None and fcst_freq is not None: da_win = fcst_freq
      if da_win is not None and fcst_freq is None: fcst_freq = da_win

      self.set_use_end_dt( arch_yml['use_end_dt'] )
      self.set_da_win( da_win )
      self.set_fcst_freq( fcst_freq )
      self.set_dir_type( dir_type )
      self.set_arch_root( root )
      self.set_file_info( file_dict=arch_yml[self.dir_type()], 
                          file_dict_items=info_items )
      self.set_var_delimiter( var_delim_def[0], var_delim_def[1])


      # defaults from YAML file
      # fc_type: type of forecast (an, fc, fcmm etc.)
      #     an for analysis (forecast time=0)
      #     fc for forecast (forecast time > 0, hourly data)
      #     fcmm for forecast (forecast time > 0, sub-hourly data)
      if fc_type is None:
           fc_type = arch_yml[self.dir_type()]['fc_type']
           # ... shouldn't really happen
           # default if not set in YAML file
           if fc_type is None: fc_type = 'fcmm'

      # lev_type: type of level (slv, spec, ml, pl etc.)
      #      spec for (10min) single level data,
      #      slv  for (hourly) single level data
      #      ml   for (hourly) model level data
      #      pl   for (hourly) pressure level data
      if lev_type is None:
           lev_type = arch_yml[self.dir_type()]['lev_type']
           # default if not set in YAML file
           # ... shouldn't really happen
           if lev_type is None: lev_type = 'sfc'

      pert_id = arch_yml[self.dir_type()]['pert_id']

      self.set_lev_type( lev_type )
      self.set_fc_type( fc_type )
      self.set_pert_id( pert_id )

      if dtfmt is None:
         self.__dt_fmt = dtf_grid
      else:
         self.__dt_fmt = dtfmt
      return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def __dir__(self):
       list_dir = [ \
                   '__dir__',  \
                   '__init__', \
                   '__da_win', '__fcst_freq', '__ic_info', '__bdt_info', \
                   '_dir_type','_lev_type', '_arch_root', '_fc_type', \
                   '_file_info', '_pert_id', \
                   '_use_end_dt' \
                   'bdt_info', \
                   'da_win', 'da_time_from_dt', \
                   'expand_vars', \
                   'fcst_freq', \
                   'get_var_delimiter', \
                   'ic_info', \
                   'set_da_win', \
                   'set_dir_type','dir_type', \
                   'set_fcst_freq', \
                   'set_lev_type','lev_type', \
                   'set_arch_root','arch_root', \
                   'set_fc_type','fc_type', \
                   'set_file_info','file_info', 'file_attr', \
                   'set_use_end_dt','use_end_dt' \
                   'set_var_delimiter', \
                   'valid_dt','base_dt', \
                   'grid_file'  \
                  ]
       return list_dir

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   def file_name( self, var ):

   # returns full path-name that can be passed to datetime.strftime e.g. in grid_file(...)
   # for a given variable name (var)

      if 'fc_type' in self.file_info()[var]:
          # some variables only exist for certain forecast types
          #  i.e. topography and land sea mask only available under 'an'
          fc_type = self.file_info()[var]['fc_type']
      else: 
          fc_type = self.fc_type()

      if 'lev_type' in self.file_info()[var]:
          # some variables only exist for certain forecast types
          #  i.e. topography and land sea mask only available under 'slv'
          lev_type = self.file_info()[var]['lev_type']
      else: 
          lev_type = self.lev_type()
   
      fname = self.file_info()[var]['fname']
      fvar_name = self.file_info()[var]['fvar_name']
      logger.debug(f"var: {var},  fname: {fname}")
      if fname[0] == '/':
        fpath = fname
      else:
        if self.dir_type() == 'mars_ncdf.0' : 
          fpath = os.path.join( self.arch_root(),'%Y%m%d/%H00/',fc_type,
                               lev_type,fname)

        elif self.dir_type() == 'mars_ncdf.1' : 
          fpath = os.path.join( self.arch_root(),'%Y%m%d/%H00/',fc_type,
                               lev_type,fname)

        elif self.dir_type() == 'iris_ncdf.0' : 
          fpath = os.path.join( self.arch_root(),'%Y%m%dT%H00Z/nc/',
                               lev_type,
                            '-'.join([fname,lev_type,'%Y%m%dT%H00Z.nc']) )

        elif self.dir_type() == 'era5_ncdf.0' : 
          ens_id = os.getenv('ENS_ID')
          # if ens_id will evaluate false for None, '', False etc.
          if ens_id :
              if self.pert_id() is not None:
                 ens_id = ''.join([self.pert_id(),ens_id])
              fpath = os.path.join( self.arch_root(),'%Y/%m/%d',
                             '.'.join([fname,ens_id,lev_type,'%Y%m%d%H.nc']) )
          else:
              fpath = os.path.join( self.arch_root(),'%Y/%m/%d',
                             '.'.join([fname,lev_type,'%Y%m%d%H.nc']) )

        elif self.dir_type() == 'barra2_ncdf.0' : 
          ens_id = os.getenv('ENS_ID')
          # if ens_id will evaluate false for None, '', False etc.
          if ens_id :
              if self.pert_id() is not None:
                 ens_id = ''.join([self.pert_id(),ens_id])
#             fpath = os.path.join( self.arch_root(),'%Y%m%dT%H00Z',ens_id,lev_type,fname)
              fpath = os.path.join( self.arch_root(),'%Y/%m/%Y%m%dT%H00Z',ens_id,lev_type,fname)
          else:
              logger.debug(f"lev_type: {lev_type},   fname: {fname}")
#             fpath = os.path.join( self.arch_root(),'%Y%m%dT%H00Z',lev_type,fname)
              fpath = os.path.join( self.arch_root(),'%Y/%m/%Y%m%dT%H00Z',lev_type,fname)

        elif self.dir_type() == 'barra1_ncdf.0' : 
              fpath = os.path.join(self.arch_root(),fc_type,lev_type,fvar_name,'%Y/%m',fname)

        else:
            fpath = os.path.join( self.arch_root(),fname)

      logger.info(f"fpath: {fpath}")
      # not sure what os.path.expandvars will do with '{{ var }}stuff'
      fpath = self.expandvars(fpath, var)
      logger.info(f"file_info expanded fpath: {fpath}")
      fpath = os.path.expandvars(fpath)
      logger.info(f"environment expanded fpath: {fpath}")

      return fpath
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def dir_type(self):
       return self._dir_type
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_dir_type(self, str):
       self._dir_type = str
       return
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def lev_type(self):
       return self._lev_type
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_lev_type(self, str):
       self._lev_type = str
       return
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def arch_root(self):
       return self._arch_root

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_arch_root(self, str):
       self._arch_root = str
       return
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def fc_type(self):
       return self._fc_type

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_fc_type(self, str):
       self._fc_type = str
       return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def pert_id(self):
       return self._pert_id

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_pert_id(self, str):
       self._pert_id = str
       return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def file_attr(self, key=None, var=None):
       # Returns attribute (key) for variable (var) 
       #    returns default value if var=None

       if var is None:
          if key in self.file_info():
              return self.file_info()[key]
          else:
              logger.error(f'Unknown default key: {key}')
              raise ValueError

       elif key in self.file_info()[var]:
          return self.file_info()[var][key]
       elif key in self.file_info():
          return self.file_info()[key]
       else:
          logger.error(f'Unknown (key,var) pair:[{key}] [{var}]')
          raise ValueError

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def file_info(self):
      return self._file_info

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_file_info(self, file_dict=None,
                           file_dict_items=None):
      # can be used to set up the whole dictionary of meta data
      # using file_dict
      #  and/or
      # update individual items using file_dict_items

      if file_dict is not None:
         self._file_info = file_dict

      if file_dict_items is not None:
        self._file_info = upd_dict( self._file_info, file_dict_items)

      return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_da_win(self, da_win):
      day_sec = 24*60*60
      self._da_win = da_win

      self.__ic_info = {}
      # set up start times in seconds
      #      da_win already in minutes
      # assume first analysis time spans days
      for h in np.arange(0,24,round(da_win/60)):
          ic_mid = h*60*60
          ic_start = ic_mid - 30*da_win 
          if ic_start < 0: ic_start += day_sec
          ic_end = ic_mid + 30*da_win
          if ic_end >= day_sec: ic_end -= day_sec
          if self.use_end_dt():
             ic_end += 1.e-3
      # Catch first interval, assumed to span the change of day
          self.__ic_info[h] = {'start_sec':ic_start, 
                               'mid_sec':ic_mid,
                               'ref_hr':h,
                               'end_sec':ic_end}
      return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_fcst_freq(self, fcst_freq):
      day_sec = 24*60*60
      self._fcst_freq = fcst_freq

      self.__bdt_info = {}
      # set up start times in seconds
      #      fcst_freq already in minutes
      # assume first analysis time spans days
      for h in np.arange(0,24,round(fcst_freq/60)):
          ic_mid = h*60*60
          ic_start = ic_mid - 30*fcst_freq 
          if ic_start < 0: ic_start += day_sec
          ic_end = ic_mid + 30*fcst_freq
          if ic_end >= day_sec: ic_end -= day_sec
          if self.use_end_dt():
             ic_end += 1.e-3
          self.__bdt_info[h] = {'start_sec':ic_start, 
                               'mid_sec':ic_mid,
                               'ref_hr':h,
                               'end_sec':ic_end}
      return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def da_win(self):
      return self.__da_win

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def fcst_freq(self):
      return self.__fcst_freq

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def da_time_from_dt(self, ref_dt, use_basedt=False):
      if use_basedt:
         ic_info = self.bdt_info()
      else:
         ic_info = self.ic_info()
      ref_sec = (ref_dt.hour*60 + ref_dt.minute)*60 + ref_dt.second
      da_time = None
      

      # Catch first interval, assumed to span the change of day
      if ref_sec >= ic_info[0]['start_sec']:
         # analysis time is "tomorrow"
         da_time = (ref_dt+delt(days=1)).replace(hour=ic_info[0]['ref_hr'])
      elif ref_sec < ic_info[0]['end_sec'] :
         da_time = ref_dt.replace(hour=ic_info[0]['ref_hr'])
      else:
        for h in ic_info.keys():
          if h == 0 : continue   # just in case order of keys gets mangled
          if ref_sec >= ic_info[h]['start_sec'] and \
             ref_sec  < ic_info[h]['end_sec']:
                da_time = ref_dt.replace(hour=ic_info[h]['ref_hr'])
                break
               
      if da_time is not None:
         da_time = da_time.replace(minute=0,second=0,microsecond=0)
      else:
         logger.error(f"Could not find window containing datetime: {ref_dt}")
         if use_basedt:
           logger.error(f"check archive.bdt_info()")
         else:
           logger.error(f"check archive.ic_info()")
         raise ValueError
      
      return da_time

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def ic_info(self):
      return self.__ic_info

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def bdt_info(self):
      return self.__bdt_info

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def valid_dt(self, f, var_name):
       """
       Returns valid datetimes for file containing variable 'var_name'
       f : cdms file handle
       var_name : variable being extracted from gridfile
       """

       if 'no_date_check' in self.file_info()[var_name] : return None

       vt_ncname = self.file_attr('vtime_var',var_name)
       dt_type   = self.file_attr('time_type',var_name)

       if dt_type == 'date_time' :
          dt_ncname = self.file_attr('vdate_var',var_name)
          fvtime    = f(vt_ncname)
          fvtime    = np.array(fvtime, dtype=np.int64)
          fvdate    = f(dt_ncname)
          fvdate    = np.array(fvdate, dtype=np.int64)
          fval_dt   = self.__dt_fmt.int_to_dt(fvdate, fvtime)

       elif dt_type == 'timestamp':
          vt_ref_fmt = self.file_attr('vtime_ref_format', var_name)
          fval_dt = get_nc_time(f[vt_ncname], vt_ref_fmt)

       return fval_dt

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def base_dt(self, f, var_name):
       """
       Returns base datetimes for file containing variable 'var_name'
       f : cdms file handle
       var_name : variable being extracted from gridfile
       """

       if 'no_date_check' in self.file_info()[var_name] : return None

       dt_type   = self.file_attr('time_type', var_name)
       bt_ncname = self.file_attr('btime_var',var_name)

       if dt_type == 'date_time' :
          dt_ncname = self.file_attr('bdate_var', var_name)
          fbtime    = f(bt_ncname)
          fbtime    = np.array(fbtime, dtype=np.int64)
          fbdate    = f(dt_ncname)
          fbdate    = np.array(fbdate, dtype=np.int64)
          fbase_dt  = self.__dt_fmt.int_to_dt(fbdate, fbtime)

       elif dt_type == 'timestamp' :
          bt_ref_fmt = self.file_attr('btime_ref_format', var_name)
          vt_ncname = self.file_attr('vtime_var', var_name)
    
          if bt_ncname == vt_ncname:
             source='units'
          else:
             source = None
          fbase_dt = get_nc_time( f[bt_ncname], bt_ref_fmt, source)

       return fbase_dt

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def grid_file(self, var_name, fname_dt, valid_dt=None, base_dt=None, rw='r'):
      """
       Returns handle for file containing requested grid (file_handle), 
          name of file (fname), and
          name of grid (gname)
   `
       Name generated from grid_file in config file
         var_name : name of required gridfield (grid_file key)
         fname_dt : datetime value to be editted into generic filename
                        - generally basetime, but does not need to be
         valid_dt : list of datetimes required in file
                      only check first and last entries
                      range of values in file can be wider, but must contain valid_dt
         base_dt  : as for valid_date but for base datetimes
         rw       : access method, i.e. 'r', 'r+'
      """

      fname = self.file_name(var_name)
      logger.debug(f"var_name: {var_name},   fname: {fname}")

      if fname_dt is not None :
         fname = dt.strftime(fname_dt,fname)

      fname_glob = glob(fname)
      if len(fname_glob) > 1:
          logger.error(f'File name is ambiguous: {fname} resolved to {fname_glob}')
          raise ValueError
      elif len(fname_glob) == 0:
          logger.error(f'File not found: {fname}')
          raise ValueError
      fname = fname_glob[0]
      
      logger.debug(f'grid file: {fname}')
      if not os.path.exists(fname):
         logger.error(f'file {fname} does not exist')
         raise ValueError

      file_handle = cdms.open(fname)

      if valid_dt is not None:
         fval_dt = self.valid_dt(file_handle, var_name)
         if fval_dt[0] > valid_dt[0] :
             logger.error(f'valid datetime mismatch in {fname}')
             logger.error(f'expected first datetime {valid_dt[0]}')
             logger.error(f'first datetime found {fval_dt[0]}')
             raise ValueError
         if fval_dt[-1] < valid_dt[-1] :
             logger.error(f'valid datetime mismatch in {fname}')
             logger.error(f'expected last datetime {valid_dt[-1]}')
             logger.error(f'last datetime found {fval_dt[-1]}')
             raise ValueError

      if base_dt is not None:
         bval_dt = self.base_dt(file_handle, var_name)
         if bval_dt[0] > base_dt[0] :
             logger.error(f'base datetime mismatch in {fname}')
             logger.error(f'expected first datetime {base_dt[0]}')
             logger.error(f'first datetime found {bval_dt[0]}')
             raise ValueError
         if bval_dt[-1] < base_dt[-1] :
             logger.error(f'base datetime mismatch in {fname}')
             logger.error(f'expected last datetime {base_dt[-1]}')
             logger.error(f'last datetime found {bval_dt[-1]}')
             raise ValueError

      return file_handle,fname,self.file_info()[var_name]['fvar_name']

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_use_end_dt(self, flag):
      self._use_end_dt = flag
      return
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def use_end_dt(self):
      return self._use_end_dt 
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_var_delimiter(self, start_str, end_str):
      """
        Set delimiting strings to identify internal variables in filenames
      """
      self.start_str = start_str
      self.end_str = end_str
      return
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def get_var_delimiter(self):
      """
        Return delimiting strings to identify internal variables in filenames
      """
      return self.start_str, self.end_str
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def expandvars(self, fpath, field):
      """
        expand internal variables in pathname
        variable names surrounded by specific character strings (var_delimiter)
 
        variable names are read from a dictionaary (fdict)
      """

      start_var, end_var = self.get_var_delimiter()
      if start_var in fpath:
          iv0 = fpath.index(start_var)
          iv1 = fpath[iv0:].index(end_var)
          f0 = fpath[:iv0]         # filepath before variable to be replaced
          f1 = fpath[iv0+iv1+2:]   # filepath after variable to be replaced
          var = fpath[iv0+2:iv0+iv1].strip()   # variable to be replaced

          if var in self.file_info()[field]:
              var = self.file_info()[field][var]
          elif var in self.file_info():
              var = self.file_info()[var]
          else:
              logger.debug(f"Keys for field {field}: {self.file_info()[field].keys()}")
              logger.debug(f"Keys for file_info: {self.file_info().keys()}")
              raise ValueError

          # Reconstruct fiepath
          new_path = f0 + var + f1
          # Check if any more variables to replace
          if start_var in new_path: new_path = self.expandvars(new_path, field)
      else:
          new_path = fpath
          
      return new_path
