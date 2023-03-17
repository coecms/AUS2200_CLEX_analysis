#!/usr/bin/env python

def plot_timeseries(st_id,
                    var_name,
                    date_start,
                    date_end):
    def setup(ax):
        """Set up common parameters for the Axes in the example."""
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(60*60*24))
        ax.xaxis.set_ticks_position('bottom')
        # ax.set_xlabel("s")
        ax.set_ylabel(f"{var_name} "+"["+data[f"obs_{var_name}"].units+"]")

    startdate=dt.strptime(date_start,"%Y%m%dT%H%M")
    enddate=dt.strptime(date_end,"%Y%m%dT%H%M")
    filenames = [fname for fname in glob.glob("/scratch/public/pjs548/aus2200_v2/obs_int/*.nc") if (((dt.strptime(os.path.split(fname)[-1][-17:-3],"%Y%m%dT%H%MZ")-startdate).total_seconds()>=0) and ((dt.strptime(os.path.split(fname)[-1][-17:-3],"%Y%m%dT%H%MZ")-enddate).total_seconds()<=0))]
    for _id in st_id:
        plt.figure()
        obs_data=[]
        model_data=[]
        time=[]
        for fname in filenames:
            time0=(dt.strptime(os.path.split(fname)[-1][-17:-3],"%Y%m%dT%H%MZ")-startdate).total_seconds()
            data=xr.open_dataset(fname)
            time=np.append(time,[time0 + x for x in data.time[data.station_identifier.values==_id]])
            obs_data=np.append(obs_data,data[f"obs_{var_name}"][data.station_identifier.values==_id])
            model_data=np.append(model_data,data[f"adj_{var_name}"][data.station_identifier.values==_id])
        idx=np.argsort(time)
        time = time[idx]
        # obs_data = obs_data[idx]
        # model_data = model_data[idx]
        ax=plt.axes()
        # ax.plot(time,obs_data,label=f"Observations")
        # ax.plot(time,model_data,label=f"Aus2200")
        ax.plot(obs_data,label=f"Observations")
        ax.plot(model_data,label=f"Aus2200")
        setup(ax)
        plt.title(f"Station id: {_id}")
        plt.legend()
    plt.show()
    
description='''Plot the timeseries of observed data and modelled data 
    (interpolated over station's heights) for a specific station id
    and for a specific time window.'''

import argparse

parser = argparse.ArgumentParser(description=description, allow_abbrev=False)
parser.add_argument('--id', '--station_id', dest='st_id', required=True, nargs='+', type=str,
                    help='Station ID.')
parser.add_argument('-v','--var', '--var_name', dest='var_name', required=True, type=str,
                    help='Variable name.')
parser.add_argument('-s','--start', '--start_date', dest='date_start', required=True, type=str,
                    help="Start date in the format 'YYYYmmddTHHMM', where YYYY is the year, "
                        "mm is the month number, dd is the day number, 'T' is the letter 'T' "
                        "used as a date/time divider, HH is the hours, MM is the minutes.")
parser.add_argument('-e','--end', '--end_date', dest='date_end', required=True, type=str,
                    help="End date in the format 'YYYYmmddTHHMM', where YYYY is the year, "
                        "mm is the month number, dd is the day number, 'T' is the letter 'T' "
                        "used as a date/time divider, HH is the hours, MM is the minutes.")

args = parser.parse_args()
st_id=args.st_id
var_name=args.var_name
date_start=args.date_start
date_end=args.date_end

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from datetime import datetime as dt
from matplotlib import ticker
import re

plot_timeseries(st_id,var_name,date_start,date_end)
