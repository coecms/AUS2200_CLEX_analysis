#####
#####
#####  DEPRECATED AND TO BE DELETED
#####      - just available for transition
#####
#####

mod_name='get_config_file'
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def check_file (fname, debug=None):
       """
           check if file exists and can be read
       """
       sr_name = mod_name+'.check_file >> ' 
       if not debug is None:
          if debug > 0:
              print(sr_name+'fname  : ',fname)
              print(sr_name+'isfile : ',os.path.isfile(fname))
              print(sr_name+'access : ',os.access(fname,os.R_OK))

       return os.path.isfile(fname) and os.access(fname,os.R_OK)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_config_file ( f_default=None, f_run=None, debug=None ):
   """
      Returns the names for default configuration files and
         run specific files
   
      f_default = filename of default configuration file
      f_run     = filename for configuration of this particular run
     
      Returns at most 3 filenames in a list
         in the order they should be executed via execfile
         file_list[0]  f_default from PYTHONPATH
         file_list[1]  f_default from current directory
         file_list[2]  f_run from current directory

      Usage:
         import get_config_file as ini

         ini_f = ini.get_config_file( 'file.def', 'file.cfg' )
         for f in ini_f:
           execfile(f)

   """

   sr_name = mod_name+'get_config_file'
   py_path        = os.getenv('PYTHONPATH').split(':')
   default_file_l = [ d+'/'+f_default for d in py_path ]

   file_list=[]

   # check PYTHONPATH for default files
   for f in default_file_l:
       if check_file(f, debug) :
          file_list.append(f)
          break
   
   # check current directory for default settings
   #   will be assumed to override PYTHONPATH defaults
   
   if check_file( f_default, debug ):
      file_list.append(f_default )
   
   # check current directory for run-time specific settings
   
   if not(f_run is None) and check_file( f_run, debug ):
      file_list.append(f_run )
   
   if len(file_list) < 1:
      print(sr_name,'  WARNING >> No files found ')
      print('     f_default = [',f_default,']')
      print('     f_run     = [',f_run,']')
   
   return file_list

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
