#!/usr/bin/python2.7

import argparse
import datetime as dt
import os
from ConfigParser import SafeConfigParser

local_dir=os.getcwd()

parser = argparse.ArgumentParser()

parser.add_argument("--gsflow", dest="gsflow_bin",
                    required=True, help="Full path of gsflow binary")
parser.add_argument("--dem", dest="dem",
                    required=True, help="Name of dem file (without extension)")
parser.add_argument("--start", dest="start_date",
                    required=True, help="start date (YYYY-mm-dd)")
parser.add_argument("--end", dest="end_date",
                    required=True, help="end date (YYYY-mm-dd)")
parser.add_argument("--name", dest="proj_name", default="No Name Project",
                    required=False, help="Name of project")
parser.add_argument("--dir", dest="gsflow_dir", default=local_dir + "/GSFLOW/",
                    required=False, help="gsflow project directory")

args = parser.parse_args()

try:
    start_date = dt.datetime.strptime(args.start_date, "%Y-%m-%d")
    end_date = dt.datetime.strptime(args.end_date, "%Y-%m-%d")
except:
    raise
                    
settings_ini = local_dir + "/settings.ini"

cfg_fid = open(settings_ini, "w+")

Config = SafeConfigParser()

control_dir = args.gsflow_dir + "/control/"
PRMSinput_dir = args.gsflow_dir + "/inputs/PRMS/"
MODFLOWinput_dir = args.gsflow_dir + "/inputs/MODFLOW_NWT/"
PRMSoutput_dir = args.gsflow_dir + "/outputs/PRMS/"

Config.add_section("settings")
Config.set("settings", "proj_name", args.proj_name)
Config.set("settings", "local_dir", local_dir)
Config.set("settings", "gsflow_exe", args.gsflow_bin)
Config.set("settings", "control_dir", control_dir)
Config.set("settings", "PRMSinput_dir", PRMSinput_dir)
Config.set("settings", "MODFLOWinput_dir", MODFLOWinput_dir)
Config.set("settings", "PRMSoutput_dir", PRMSoutput_dir)
Config.set("settings", "DEM", args.dem)
Config.set("settings", "start_date", args.start_date)
Config.set("settings", "end_date", args.end_date)

Config.write(cfg_fid)

cfg_fid.close()

for dirname in [control_dir, PRMSinput_dir, MODFLOWinput_dir, PRMSoutput_dir]:
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

open(local_dir + "/python_scripts/MODFLOW_scripts/__init__.py", "w+").close()
