#!/usr/local/anaconda3/envs/chil_3_8/bin/python
#20200916 - metsensors_ncas_v3 has had all tabs replaced by spaces, but is still working in v2.7
# =========================================================================
# GENERATE DAYS OF NETCDF FOR METSENSORS 
# Produces NetCDF data from up to a month of Halo lidar processed files 
#
# Author:        Judith Jeffery, RAL 
# History:        
# Version:	 0.1
# Last modified: 09/08/17 
# =========================================================================

import getopt, sys
#This makes it possible to produce plots without an X window, e.g. from a cronjob
#Needs to be in here, not just in any function that's imported
import matplotlib as mpl
mpl.use('Agg')
import metsensors_ncas

import datetime
from datetime import date, datetime
from pylab import date2num

print("Python version in code = ", sys.version)
print("Matplotlib version:")
mpl.__version__

# -------------------------
# Explain function usage
# Needs to be before main()
# -------------------------
def usage():

    print("\n")
    print("Command line structure")
    print("./generate_days_netcdf_metsensors.py -s yyyymmdd (of start day) -e yyyymmdd (of finish day) -x sensor choice")
    print('x = instrument number')
    print('1 = Pluvio weighing raingauge')
    print('2 = Campbell Scientific receive cabin datalogger 1 - met. sensors and raingauges')
    print('3 = Chilbolton raingauges from format5 files (before 20200728)')
    print('4 = Chilbolton broadband radiometers from format5 files')
    print('5 = Campbell Scientific Sparsholt datalogger - raingauges')
    print('6 = Chilbolton disdrometer')
    print('7 = Sparsholt disdrometer')
    print('8 = Campbell Scientific PWS100 present weather sensor')
    print('9 = Chilbolton broadband radiometers from Campbell datalogger files')
    print('10 = Chilbolton CNR4 net flux radiometer')
    print('11 = Chilbolton Metek USA-1 sonic anemometer')
    print (' ')



# ------------------------------------------------------------------------
#def main():
# ------------------------------------------------------------------------

if len(sys.argv) <= 1:  #no arguments provided
    usage()
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:e:x:")
    print (' ')
    print('Example -s 20200105 -e 20200106 -x 3')
    print (' ')
    print('x = instrument number')
    print('1 = Pluvio weighing raingauge')
    print('2 = Campbell Scientific receive cabin datalogger 1 - met. sensors and raingauges')
    print('3 = Chilbolton raingauges from format5 files (before 20200728)')
    print('4 = Chilbolton broadband radiometers from format5 files')
    print('5 = Campbell Scientific Sparsholt datalogger - raingauges')
    print('6 = Chilbolton disdrometer')
    print('7 = Sparsholt disdrometer')
    print('8 = Campbell Scientific PWS100 present weather sensor')
    print('9 = Chilbolton broadband radiometers from Campbell datalogger files')
    print('10 = Chilbolton CNR4 net flux radiometer')
    print('11 = Chilbolton Metek USA-1 sonic anemometer')
    print (' ')
except getopt.GetoptError as err:	#Changed from v2
    # print help information and exit:
    print(str(err)) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)


# -------------------------
# Set default date as today
# -------------------------
today = date.today()
year  = today.year
month = today.month
day   = today.day
epoch_offset = 719163

# --------------------
# Command line parsing
# --------------------
for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-s", "--startday"):
            #startday = int(a)
            startday = str(a)
        elif o in ("-e", "--endday"):
            #endday = int(a)
            endday = str(a)
        elif o in ("-x", "--sensor"):
            #endday = int(a)
            sensor = int(a)
        else:
            assert False, "unhandled option"

print("start day: ", startday)
print("end day: ", endday)
print("sensor: ", sensor)

#startnum = int(date2num(datetime.datetime(int(startday[0:4]),int(startday[4:6]),int(startday[6:8]),0,0,0)))
#endnum = int(date2num(datetime.datetime(int(endday[0:4]),int(endday[4:6]),int(endday[6:8]),0,0,0)))
startnum = int(date2num(datetime(int(startday[0:4]),int(startday[4:6]),int(startday[6:8]),0,0,0)))
endnum = int(date2num(datetime(int(endday[0:4]),int(endday[4:6]),int(endday[6:8]),0,0,0)))
print('Julian start, end day = ',startnum,endnum) 
#Get today's 'nday' value
str_today = datetime.strftime(datetime.now(), '%Y%m%d')     #Today's date as a string
num_today = int(date2num(datetime(int(str_today[0:4]),int(str_today[4:6]),int(str_today[6:8]),0,0,0)))
print('String, number today = ',str_today, num_today)

#Depending on version of Python, date2num can either give day wrt 01/01/0000 (earlier versions) or 01/01/1970 (later versions)
#This causes startnum, endnum to either be ~738000 or ~19000
#Other parts of code are based on 01/01/0000 so to accommodate different versions of Python, add to value if it's less than 700000

if startnum < 100000:
    startnum = startnum + epoch_offset
    endnum = endnum + epoch_offset
    num_today = num_today + epoch_offset


for nday in range(startnum,endnum+1):

    datevals = metsensors_ncas.generate_netcdf_common(nday)
    print('\n')
    print('Date being generated = ',datevals[3])
    print('\n')

    if sensor == 1:
        metsensors_ncas.generate_netcdf_pluvio(nday)
    elif sensor == 2:
        metsensors_ncas.generate_netcdf_met(nday,num_today,sensor)
    elif sensor == 3:
        metsensors_ncas.generate_netcdf_rain_f5(nday)
    elif sensor == 4:
        metsensors_ncas.generate_netcdf_bbrad_f5(nday)
    elif sensor == 5:
        #metsensors_ncas_trial.generate_netcdf_sparsholt_raingauge(nday,num_today)
        metsensors_ncas.generate_netcdf_met(nday,num_today,sensor)
    elif sensor == 6 or sensor == 7:
        metsensors_ncas.generate_netcdf_disdro_f5(nday,sensor)
    elif sensor == 8:
        metsensors_ncas.generate_netcdf_pws100(nday)
    elif sensor == 9:
        metsensors_ncas.generate_netcdf_met(nday,num_today,sensor)
    elif sensor == 10:
        metsensors_ncas.generate_netcdf_cnr4_netflux(nday)
    elif sensor == 11:
        metsensors_ncas.generate_netcdf_sonic(nday)
    else:
        print("Sensor option does not exist")



