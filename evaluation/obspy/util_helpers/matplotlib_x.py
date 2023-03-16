# A few routines to assist with using matplotlib

import matplotlib.pyplot as plt
from numpy import floor, ceil, log10, arange, abs

mod_name = "matplot_extras."

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_new_plot_area( fig, iplot=0, ncol=2, nrow=2, header="", font_size=12,
                        left=0.10, bottom=0.15, right=0.95, top=0.90,
                        wspace=0.05, hspace=None,horizontalalignment='center',
                        projection=None,title_vspace=0.05, title_hspace=0.05,
                        title_sep=',', title_split_len=150):
    """
    set up a new plotting area
    """

    new_page = (iplot == ncol*nrow) or ( iplot == 0 )
    
    if new_page:
        iplot = 1
        fig   = plt.figure()
        if len(header) > title_split_len:
            hd_split = header.split(',')
            hd2 = title_sep.join(hd_split[1:])

            fig.subplots_adjust(left=left, bottom=bottom, right=right, \
                                top=top-title_vspace,
                                wspace=wspace, hspace=hspace)
            fig.text(0.5-title_hspace, 1-title_vspace, \
                     hd_split[0], horizontalalignment=horizontalalignment, \
                     fontsize=font_size)
            fig.text(0.5-title_hspace, 1-2*title_vspace,   \
                     hd2, horizontalalignment=horizontalalignment, \
                     fontsize=font_size)
        else:
            fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                                wspace=wspace, hspace=hspace)
            fig.text(0.5-title_hspace, 1-title_vspace, \
                     header, horizontalalignment=horizontalalignment, \
                     fontsize=font_size)

    else:
        iplot +=1

    ax = plt.subplot(nrow,ncol,iplot,projection=projection)

    end_page = (iplot == ncol*nrow)

    return fig, ax, iplot, end_page

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def scale_val(x, val=None):
   """
   Provide reasonable max & min scaling values
     -leading digit is 1, 2 or 5
   """

   ax = abs(x)
   if x==0:
     sgn=1
   else:
     sgn = ax/x

   if val is None:
      val = 'max'
   else:
      val = val.lower()

   if ax > 1.e-8:
      i = float( 10**int(log10(ax)) )
   else:
      i = 0

   if ax < 1.:
      i = i/10.

   if val == 'min':
      #loop through again, but change sign to allow for min to now be like a max
      if ax >= 9*i :
         return sgn*9*i, i
      elif ax >= 8*i :
         return sgn*8*i, i
      elif ax >= 5*i:
         return sgn*5*i, i/2.
      else:
         return sgn*i, i/10.
      
   else:
      if i >= ax:
         return sgn*i, i/10.
      elif 2*i >= ax:
         return sgn*2*i, i/2.
      elif 5*i >= ax:
         return sgn*5*i, i
      else:
         return sgn*10*i, i

def scale_2val(mny,mxy):
    """
    Get nice axis values
    """
    sr_name = mod_name+"scale_2val >> "
    rng=mxy-mny
    try:
      nsig=-floor(log10(rng))
    except:
      print(sr_name,"mxy,mny,rng = ",mxy,mny,rng)
      nsig=1
    ymin=floor(mny*10**nsig)/10**nsig
    ymax=ceil(mxy*10**nsig)/10**nsig
    
    if abs(ymax-ymin) < 1.e-10 :
        ymin = 0.9*ymax

    yinc=float(10**int(log10(ymax-ymin)))
    if yinc>rng:
       yinc=yinc/10
    elif (2*yinc)>rng:
       yinc=yinc/2
    
    return ymin,ymax,yinc

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_state_region( state_name ):
    """
    Returns regular lat,lon
    domains that cover various Australian states
    """
    sr_name = mod_name+"get_state_region >>"
    # Northern territory
    if state_name[:2] in [ 'No', 'NO', 'NT', 'nt', 'no'] :
      domain = [-30., -9., 125., 141.]
 
    # Adelaide
    elif state_name[:2]  in ['AD', 'Ad', 'ad'] :
      domain = [-39.75, -30.3, 131.7, 142.25 ]
 
    # Australia
    elif state_name[:1]  in ['A', 'a'] :
      domain = [-45., -9., 108., 155.]
 
    # Brisbane
    elif state_name[:1]  in ['B', 'b'] :
      domain = [-31.25, -21.8, 147.7, 156.25]
 
    # NSW
    elif state_name[:1]  in ['N', 'n'] :
      domain = [-38., -26., 139., 155.]
 
    # Perth
    elif state_name[:1]  in ['P', 'p'] :
      domain = [-37.25, -27.8, 111.7, 120.25]
 
    # Queensland
    elif state_name[:1]  in ['Q', 'q'] :
      domain = [-32., -9., 135., 155.]
 
    # Sydney
    elif state_name[:2]  in ['Sy', 'SY', 'sy'] :
      domain = [-38.25, -29.8, 146.7, 155.25]
 
    # SA
    elif state_name[:1]  in ['S', 's'] :
      domain = [-42., -23., 125., 145.]
 
    # Tasmania
    elif state_name[:1]  in ['T', 't'] :
      domain = [-44., -38., 143., 150.]
 
    # Victoria
    elif state_name[:1]  in ['V', 'v'] :
      domain = [-42., -33., 139., 151.]
 
    # WA
    elif state_name[:1]  in ['W', 'w'] :
      domain = [-37., -12., 110., 135.]
    else:
      print(sr_name,'Unknown state = ',state_name)
      sys.exit()
 

    return domain

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def  get_latlon_lst( domain, grid_res, grid_prec=0.0001 ):
   """
   Returns vectors of latitudes and longitudes
     with a given range & spacing
   These can be passed to ax.set_xticks and ax.set_yticks
        grid_prec: precision of lat/lon values
   """

   grid_res = grid_res / grid_prec
   lat_min = domain[0] / grid_prec
   lat_max = domain[1] / grid_prec
   lon_min = domain[2] / grid_prec
   lon_max = domain[3] / grid_prec

   iln1   = ( int(lon_min-0.001) / int(grid_res) + 1 ) * int(grid_res)
   ilnn   = ( int(lon_max-0.001) / int(grid_res) + 1 ) * int(grid_res)
   ilt1   = ( int(lat_min-0.001) / int(grid_res)     ) * int(grid_res)
   iltn   = ( int(lat_max-0.001) / int(grid_res) + 1 ) * int(grid_res)
   
   lonLst = arange( float(iln1), ilnn, grid_res) * grid_prec
   latLst = arange( float(ilt1), iltn, grid_res) * grid_prec

   return latLst, lonLst

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
