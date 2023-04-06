from numpy import radians
from numpy.ma import cos, sin

from .defns import logger

def ffddd_to_uv(ff,ddd):

     """
     Convert wind speed and direction to wind components.
       ff: wind speed
       dd: wind direction (angle between the northern axis and direction where wind blows from, in degrees)

      Output: u,v wind components in the natural coordinate system.
     """

     phi = radians(270.-ddd)

     u = ff*cos(phi)
     v = ff*sin(phi)
     return u,v
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ffddd_to_u(ff,ddd):
     u,v = ffddd_to_uv(ff,ddd)
     return u
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ffddd_to_v(ff,ddd):
     u,v = ffddd_to_uv(ff,ddd)
     return v
