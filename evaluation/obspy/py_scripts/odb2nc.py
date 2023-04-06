import numpy as np
from netCDF4 import Dataset
import argparse
import sys

import obs_py as opy

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('-i',"--input",help="input file")
    parser.add_argument('-o',"--output",help="output file")
    parser.add_argument('-O',"--obs_type",help="group of obstypes read from ODB", \
                        default='surface')
    parser.add_argument('-n',"--network",help="network name, default: odb_surface", \
                        default='odb_surface')
    parser.add_argument('-r',"--region",help=f"odb_regionfilter default:{opy.odb_regionfilter_def} ", \
                        default=opy.odb_regionfilter_def)
    parser.add_argument('-v',"--verbose",help=f"logger level, default:{opy.logger.level}", \
                        default=opy.logger.level)
    
    args=parser.parse_args()
    opy.logger.setLevel(args.verbose)

    fob_info = opy.obs_var_info()
    netwk = opy.ob_network(args.network, odb_columns=opy.odb_columns[args.obs_type])
    odb_filter = opy.make_odb_typefilter(netwk.odb_obs_type()[args.obs_type], args.region)
    df = opy.odb2_df(args.input, netwk, fob_info, odb_filter, obs_type=args.obs_type)
    opy.df_nc(df, args.output, netwk, fob_info, obs_type=args.obs_type)
