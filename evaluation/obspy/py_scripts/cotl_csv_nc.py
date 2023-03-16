import argparse
import obs_py as opy

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--input", \
                    default='cotl.csv', \
                    help="input .csv file (default: cotl.csv)" )
parser.add_argument('-o','--output', \
                    default='obs.nc', \
                    help="output netcdf file name (default: obs.nc)" )
parser.add_argument('-s',"--stnlist", \
                    default='stns.csv', \
                    help="input file list of station names (default: stns.csv)")
parser.add_argument('-n','--network', \
                    default='cotl', \
                    help="name of network (default: cotl)")
parser.add_argument('-v',"--verbose",help=f"logger level, default:{opy.logger.level}", \
                        default=opy.logger.level)


args = parser.parse_args()
opy.logger.setLevel(args.verbose)

fob_info = opy.obs_var_info()
netwk = opy.ob_network(args.network)
stn_dict = opy.obs_csv_stn_info(args.stnlist, netwk, fob_info.var_dict())
df = opy.csv_df( args.input, stn_dict, netwk, fob_info)
df = opy.df_std_units(df, netwk.info, fob_info.var_dict())
opy.df_nc( df, args.output, netwk, fob_info)
