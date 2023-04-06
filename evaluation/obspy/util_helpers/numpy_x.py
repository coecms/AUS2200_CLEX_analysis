#
# Some additional routines to make numpy more usable
#  (or rather more like F90 etc.)
#
import numpy as np

def is_numpy_ndarray(x):
  # checks if array or masked array

  return isinstance(x, np.ndarray) 

def is_numpy_maarray(x):
  # checks if array or masked array

  return isinstance(x, np.ma.core.MaskedArray)

def is_numpy_array(x):
  # checks if array or masked array

  return is_numpy_ndarray(x) or is_numpy_maarray(x)

def max_loc(x):
  # For numpy arrays use x.argmax() and np.unravel_index(xmax, x.shape)
  m = max(x)
  ind_max = [i for i,j in enumerate(x) if j == m ]
  return ind_max

def min_loc(x):
  # For numpy arrays use x.argmin() and np.unravel_index(xmin, x.shape)
  m = min(x)
  ind_min = [i for i,j in enumerate(x) if j == m ]
  return ind_min

def false_loc(x):
  return [i for  i,j in enumerate(x) if not j ]

def IsEqualWithin(x, y, tol):
  xy_diff = abs(x-y)
  if xy_diff.max() < tol:
     return {'flag':True, 'val':None, 'loc':None}
  else:
     return {'flag':False, 'val':xy_diff.max(), 'loc':max_loc(xy_diff) }

def IsTrue(x):
   if all(x):
     return {'flag':True, 'val':None, 'loc':None}
   else:
     return {'flag':False, 'val':None, 'loc':false_loc(x) }

def vmax(x, y):
   val = np.where( x > y, x, y )
   return val

def vmin(x, y):
   val = np.where( x < y, x, y )
   return val

def nint(x, vtype=np.int32):
   val = np.array( np.rint(x), dtype=vtype )
   return val

def int(x,vtype=np.int32):
   val = x.astype(vtype)
   return val

def xfloat(x, vtype=np.float64):
   val = np.array( x, dtype=vtype )
   return val

def lin_intrp( xval, x, y):
   wt = (x[1]-xval)/(x[1]-x[0])
   yval = y[1] + wt*(y[0]-y[1])
   return yval

def lin_intrpx( xval, x, y):
   wt = (x[1]-xval)/(x[1]-x[0])
   print('wt, xval=  ',wt,xval)
   print('x0, x1=  ',x[0],x[1])
   yval = y[1] + wt*(y[0]-y[1])
   print('y0, y1=  ',y[0],y[1])
   return yval

