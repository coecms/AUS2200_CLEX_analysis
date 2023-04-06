import numpy as np
from datetime import datetime as dt

class dt_fmt( object ):
   """
   A class for holding information about date time formats
    - there are a number of interdepencies. This allows setting of
      values and maintaining the interdependcies

   Contains methods for converting between numpy arrays of ints and strings
   """

   def __init__( self, date_fmt='%Y%m%d', time_fmt='%H%M%S' ):

       """
       Class Constructor
       _XX are the variable values
       XX are functions to return the values

      
       Need to know what are date and time formats, so only allow
         specification of date & time and not datetime
       (i.e. cannot back out date and time formats from datetime format
       """

       self.set_vals(date_fmt=date_fmt, time_fmt=time_fmt)
       return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def __dir__(self):
       list_dir = [ \
                   '__dir__',  \
                   '__init__', \
                   '_date_fmt','_dt_fmt', \
                   '_time_fmt', \
                   '_scale_yr,', \
                   'date_fmt','dt_fmt','dt_to_int', 'scale_yr', \
                   'int_to_dt', \
                   'set_vals', \
                   'time_fmt' \
                  ]
       return list_dir

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_vals(self, date_fmt=None, time_fmt=None):
       """
       Sets parameters, preserving relationships
       """
       if not(date_fmt is None):
         self._date_fmt = date_fmt

       if not(time_fmt is None):
         self._time_fmt = time_fmt
 
       self._dt_fmt = self.date_fmt() + self.time_fmt()

       self._scale_yr = np.power( 10,len( self.time_fmt() ) )
       return
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def date_fmt(self):
       return self._date_fmt
      
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def time_fmt(self):
       return self._time_fmt
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def dt_fmt(self):
       return self._dt_fmt

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def scale_yr(self):
       return self._scale_yr
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def int_to_dt(self, int_date, int_time ):
       """
       Convert integer arrays of dates (int_date) and times (int_time)
         to datetime objects
       """
   
       len_t = len( self._time_fmt )

       dt_array = [ dt.strptime( str(d)+str(t).zfill(len_t), self.dt_fmt() )  \
                     for d,t in zip(int_date,int_time) ]
       dt_array =  np.array( dt_array )

       return dt_array

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def dt_to_int(self, date, dt_fmt=None):
       """
       Convert datetime arrays to a date (or time)
       """
    
       if dt_fmt is None:
          dt_fmt = self.dt_fmt()
   
       t_array = [dt.strftime(d,dt_fmt) for d in date]
       t_array =  np.array( t_array )
       return t_array


