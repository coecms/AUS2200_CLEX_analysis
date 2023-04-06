class miss_val(object):
   """
   A class for holding information about missing data
    - there are a number of interdepencies:
          missing values
          values used for testing if missing
          type of test to use to detect missing values
      This allows setting of values and maintaining the interdependcies

   """

   def __init__( self, xmiss=1.e+17, xmiss_scale=0.999, test_val=None, imiss=None, cmiss=None ):

       """
       Class Constructor
       _XX are the variable values
       XX are functions to return the values
    
       xmiss : missing data indicator
       xmiss_scale : scaling of xmiss to test value (xmiss_test)
       imiss : integer missing data indicator
       cmiss : string missing data indicator

       """

       self.set_vals(xmiss=xmiss, xmiss_scale=xmiss_scale, test_val=test_val, 
                     imiss=imiss, cmiss=cmiss)
       return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def __dir__(self):
      list_dir=[ \
                 '__dir__', \
                 '__init__', \
                 '_cmiss', '_imiss', '_xmiss','_xmiss_scale','_test_val', \
                 'cmiss', 'imiss', \
                 'xmiss','xmiss_scale','test_val' \
                 'set_vals', \
               ]
      return list_dir

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def set_vals( self, xmiss=None, xmiss_scale=None, test_val=None, 
                       imiss=None, cmiss=None ):

      if xmiss is not None:
          self._xmiss = xmiss
          # self._test_val will be updated later if needed
   
      if xmiss_scale is not None:
          self._xmiss_scale = xmiss_scale
          # self._test_val will be updated later if needed
   
      if test_val is not None:
          # Assume want to force in test_val, so skip xmiss_scale
          self._test_val = test_val
          self._xmiss_scale = None
      else:
          self._test_val = self.xmiss()*self.xmiss_scale()

      if imiss is not None:
          self._imiss = imiss
      else:
          self._imiss = round(self.xmiss())
    
      if cmiss is not None:
          self._cmiss = cmiss
      else:
          self._cmiss = '--------'
    
      return

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def xmiss_scale(self):
      return self._xmiss_scale

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   def test_val(self):
      return self._test_val

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def xmiss(self):
      return self._xmiss

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def imiss(self):
      return self._imiss

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   def cmiss(self):
      return self._cmiss
