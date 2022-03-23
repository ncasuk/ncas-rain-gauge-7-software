
# =========================================================================
# MET SENSORS PYTHON MODULE
# Tools to produces NetCDF data from a day of time series met data files
#
# Author:        Judith Jeffery, RAL 
# History:       Based on Python code written for Leosphere lidar by Chris Walden 
# Version:	 1.0
# Last modified: 09/08/17
# =========================================================================
module_version = 1.0

# -------------------------------------------------------------------------
# Import required tools
# -------------------------------------------------------------------------
import numpy as np
import numpy.ma as ma
import os, re, sys, getopt, shutil, zipfile, string, pwd
import netCDF4 as nc4

import datetime
import time
import calendar
import scipy.signal
from pylab import *
import module_data_object_python3
import module_distrometer_format5
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.dates as mdates
from matplotlib.ticker import LogFormatter
from matplotlib import colors
from matplotlib import rcParams

rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12

import pyart
import cmocean
import cftime
import colorcet

# ----------------------
# Define useful function
# ----------------------
def in_interval(seq,xmin,xmax):
    for i, x in enumerate(seq):
        if x>=xmin and x<xmax:
            yield i

# ------------------------------------------------------------------------
# Pluvio raingauge netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_pluvio(nday):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon

    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['pldcrg_ch','pldcrg_ch']	#Need a file for each dataset read, but can be same one. May want them separate for different corrections 
    nvar = np.size(chids)
    n = 0 	#Number of data points

    # -------------------------
    # Define various parameters
    # -------------------------

    missing_value = -1.0E+20
    epoch_offset = 719163

    # ----------------------
    # Initialize data arrays
    # ----------------------
    vals        = missing_value * np.ones((10000,nvar))
    timesecs    = np.zeros((10000))
    input_missing =  np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_moe_home2/range/pluvio/"
    path_out = "/data/amof-netCDF/ncas-rain-gauge-4/ncas-rain-gauge-4_cao_"
    data_version = "_v1.0"	#Changed temporarily so as to show this was generated with v3 python
    graph_path_out = "/data/amof-netCDF/graphs/ncas-rain-gauge-4" + data_version + "/" + str(year) + '/' + strmon + '/'
    data_product = "_precipitation"

    template_file_path = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-rain-gauge-4_metadata-template.yaml"      #YAML template

    #------------Setting up YAML data object----------------------
    #Need export PYTHONPATH=/home/jla/python/global_modules
    handler = module_data_object_python3.Handler()
    data_object_id = handler.load_template(template_file_path)
    print('Data_object_id = ',data_object_id)
    #handler.show_requirements_for_template(data_object_id)

    # ---------------------------------------------------------------
    # Loop to process raw data from current and next day's data files
    # Since first data point in file is at 00:00:00, first data point
    # contains data entirely from the previous day.
    # Hence, need to open file from next day too.
    # ---------------------------------------------------------------

    for day_incr in range(2):

        nday_file = nday + day_incr


        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        #print(nday_file, datevals)


        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "pldc_" + datestring_now + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print('Raw files = ',infiles)
            nfiles = len(infiles)
            if nfiles == 0:
                print('No files on this day')

        else:
            nfiles = 0
            print('No directory found')

        # ----------------------
        # Process data files
        # ----------------------
        print('Starting to process raw data files')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            for z in range(9):      #9 header lines, may eventually want some of them
                line = f.readline()

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #print fdata
                fdatatmp = line.split()
                #print len(line), line
                startnum = date2num(datetime.datetime(int(year_now),int(month_now),int(day_now),int(line[6:8]),int(line[9:11]),int(line[12:14])))
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum > nday and startnum <= (nday + 1):
                    timesecs[n] = 3600.0*float(line[6:8]) + 60.0*float(line[9:11]) + float(line[12:14])
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0
                    if len(line) > 45:	#There's data in the line beyond the quantities we're writing to a file
                        fdata=fdatatmp[1].split(';')
                        tempstr = fdata[0]
                        char_flag = 0
                        for m in range(len(tempstr)):
                            char_test = tempstr[m].isalpha()
                            if char_test:	#Alphabet character found
                                char_flag = char_flag + 1
                        if tempstr[0] == '+' and char_flag == 0:	#Lines with errors don't tend to start records with + 
                            vals[n,0] = float(fdata[0])
                        else:
                            vals[n,0] = missing_value
                            input_missing[n,0] = 1
                        tempstr = fdata[1]
                        char_flag = 0
                        for m in range(len(tempstr)):
                            char_test = tempstr[m].isalpha()
                            if char_test:	#Alphabet character found
                                char_flag = char_flag + 1
                        if fdata[1][0] == '+' and char_flag == 0: 
                            vals[n,1] = float(fdata[1])
                        else:
                            vals[n,1] = missing_value
                            input_missing[n,1] = 1
                    else:	#No relevant data in line or data missing
                        input_missing[n,0] = 1
                        input_missing[n,1] = 1

                    n += 1

            f.close()

        print('No. values from raw data files = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        input_missing = input_missing[0:n,:]
        timesecs   = timesecs[0:n]

        #print(timesecs[0:9])
        #print(vals[0:9,0])
        #print(vals[0:9,1])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]	#Assumes evenly spaced, may want to improve this

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        #Make any qc_flag values which were missing at input == 2. These won't have been read from the correction file
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)
        #qc_flag = qc_flag + input_missing	#Add value of 1 from any values that were set as missing when reading data
        #Can't use valid_min_max from function as won't have the input_missing values
        #valid_min_max = corr_details[1]


        #print('qc_flag shape = ',qc_flag.shape)
        #print('valid_min_max = ',valid_min_max)

        #vals_qc = np.where(qc_flag <> 2, vals, 0)
        #vals_qc = np.where(qc_flag <> 3, vals_qc, missing_value)
        #vals = np.where(qc_flag <> 3, vals, missing_value)
        vals = np.where(qc_flag != 2, vals, missing_value)


        #print(qc_flag[0:49,0])


    #--------------------------------------------
    # Sorting out day/time information
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #lengths_of_dimensions = time_details[8]
        #substitutions = time_details[9]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])
        #print(first_last_datetime)

        #print('vals shape = ',vals.shape)

        # -------------------------------
        # Open new netCDF file for output
        # -------------------------------
        #cfarr_head = 'amof-pluvio_chilbolton_'
        cfarr_head = 'ncas-rain-gauge-4_cao_'
        #out_file   = os.path.join(path_out,cfarr_head+datestring+'.nc')
        out_file   = path_out+datestring+data_product+data_version+'.nc'
        print('Opening new NetCDF file ' + out_file)

        lengths_of_dimensions = {"time": n}

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec": np.amax(sc_arr), "min_rr": valid_min_max[0,0], "max_rr": valid_min_max[1,0], "min_thick": valid_min_max[0,1], "max_thick": valid_min_max[1,1], "uname": user_id, "codename": __file__, "machinename": computer_id}

        data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
        data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
        data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]      #Need to include fraction of day
        data_object["variables"]["year"]["values"][:] = yr_arr[:]
        data_object["variables"]["month"]["values"][:] = mn_arr[:]
        data_object["variables"]["day"]["values"][:] = dy_arr[:]
        data_object["variables"]["hour"]["values"][:] = hr_arr[:]
        data_object["variables"]["minute"]["values"][:] = mi_arr[:]
        data_object["variables"]["second"]["values"][:] = sc_arr[:]
        data_object["variables"]["rainfall_rate"]["values"][:] = vals[:,0]
        data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,1]
        data_object["variables"]["qc_flag_rainfall_rate"]["values"][:] = qc_flag[:,0]
        data_object["variables"]["qc_flag_thickness_of_rainfall_amount"]["values"][:] = qc_flag[:,1]
        exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
        print(exit_code)



        variables_to_plot = ["rainfall_rate", "thickness_of_rainfall_amount"]
        generate_netcdf_graphs(out_file, graph_path_out, variables_to_plot, datestring)

        print('Changing permissions of file ',out_file)
        #oscommand = "chmod g+w " + out_file
        oscommand = "chmod 775 " + out_file
        os.system(oscommand) 


# ------------------------------------------------------------------------
# Chilbolton Campbell datalogger netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_met(nday,num_today,sensor):


    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]

    # -------------------------
    # Define various parameters
    # -------------------------
    daily_file_read = 0	#Flag showing whether we are reading data from an ongoing logger file (= 0) or an already-written daily file (= 1)
    daily_file_exist = 0	#Flag showing whether a daily file already exists(= 1)
    #format_threshold1 = 737635   #Date when added raingauges to CR1000X = 20200729
    #format_threshold2 = 737867   #Date when real data from gauge near flux compound (009) = 20210318
    #if nday < format_threshold1:
        #nn_max = 3	#Just met sensors, no raingauges, before 20200729
    #elif nday >= format_threshold1 and nday < format_threshold2:
        #nn_max = 8      #Met sensors, raingauge and multiple raingauges at Chilbolton
    #else:
        #nn_max = 9

 
    missing_value = -1.0E+20
    varmin2 = 0	#Dummy values, not used if 1 variable in file
    varmax2 = 0

    # -------------------------
    # Define various parameters
    # -------------------------
    A = -4.56582625e-7
    B =  8.97922514e-5
    C = -6.95640241e-3
    D =  2.74163515e-1
    E = -6.23224724
    F = 66.1824897

    A1=0.0010295
    B1=2.391e-4
    G1=1.568e-7
    Tk = 273.15
    SBconst = 5.67e-8

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_grape_loggernet/"
    #path_out_9 = "/data/amof-netCDF/diagnostics/ncas-nocorr-temperature-rh-1/ncas-nocorr-temperature-rh-1_cao_"
    if sensor == 2:
        daily_path_out = "/data/range/daily_met/cr1000x_rxcabin_1/"
        daily_file = 'CR1000XSeries_Chilbolton_Rxcabinmet1_'
        infiles = "CR1000XSeries_Chilbolton_Rxcabinmet1.dat"
        #infiles = "CR1000XSeries_Chilbolton_Rxcabinmet1_20210804.dat"
        logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_1.txt'
        if nday < 737825:	#20210204, the day after an extra channel was added to the datafiles for rg009
            logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_1_20200729.txt'
        if nday < 737636:	#2020730, the day after extra channels were added to the datafiles for rg001, rg006, rg008, rg004
            logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_1_20191024.txt'
        print('nday, logger_details = ', nday, logger_details)
    if sensor == 5:
        daily_path_out = "/data/range/daily_met/cr1000x_sparsholt_1/"
        daily_file = 'CR1000XSeries_Sparsholt1_'
        infiles = "CR1000XSeries_Sparsholt_dc_tb.dat"
        #infiles = "CR1000XSeries_Sparsholt_dc_tb_20210113.dat"
        #infiles = "CR1000XSeries_Sparsholt_dc_tb_20200623.dat"
        logger_details = '/data/netCDF/corrections/sparsholt_1.txt'
    if sensor == 9:
        daily_path_out = "/data/range/daily_met/cr1000x_rxcabin_2/"
        daily_file = 'CR1000XSeries_Chilbolton_Rxcabinmet2_'
        infiles = "CR1000XSeries_Chilbolton2_Rxcabinmet2.dat"
        logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_2.txt'
        if nday < 738167:	#20220112, the day the CM21 and CG4 were reconnected after cal. (no change to datalogger file columns) 
            logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_2_20211015.txt'
        if nday < 738078:       #20211015, the day after an extra channel was added to the datafiles for rg009
            logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_2_20211007.txt'
            infiles = "CR1000XSeries_Chilbolton2_Table1.dat"

    daily_file   = os.path.join(daily_path_out,daily_file+datestring+'.dat')

    data_version = "_v1.0"

    path_in_file = path_in + infiles

    # ---------------------------------------------------
    # Process data files
    # Read the cumulative file from the datalogger first
    # ---------------------------------------------------
    print('Starting to read raw datalogger files')

    #Check if a daily file already exists. If it does, we want to avoid overwriting it by opening it for writing
    if os.path.isfile(daily_file):	#A daily file exists
        daily_file_exist = 1


    #file_vals = read_cr1000x_chilbolton1(path_in_file, daily_file, daily_file_read, daily_file_exist, nday, num_today,nvar)
    file_vals = read_cr1000x_general(path_in_file, daily_file, daily_file_read, daily_file_exist, logger_details, nday, num_today)
    #Returns local n, timesecs, vals, chosen_chids, chosen_varnames, chosen_qcflagnames, ordered_chosen_files, sorted_indicies
    n = file_vals[0]
    timesecs = file_vals[1]
    vals = file_vals[2]
    chids = file_vals[3]
    varnames = file_vals[4]
    qcflagnames = file_vals[5]
    files_to_write = file_vals[6]
    #sorted_indicies = file_vals[7]
    chosen_files = file_vals[7]

    print('No. values from current datalogger file  = ', n)

    # ---------------------------------------------------
    # If you don't find any data from that day in the cumulative file, 
    # look for a daily file instead.
    # This case is most likely if you are processing older data which have
    # been removed from the cumulative file to reduce its size.
    # In this case, look for a previously written daily file
    # and don't try to write a daily file again. 
    # ---------------------------------------------------

    if n == 0: #No data found so far
        if os.path.isfile(daily_file):	#A daily file exists
            print('Reading from a daily file')
            daily_file_read = 1
            #file_vals = read_cr1000x_chilbolton1(daily_file, daily_file, daily_file_read, daily_file_exist, nday, num_today, nvar)
            file_vals = read_cr1000x_general(daily_file, daily_file, daily_file_read, daily_file_exist, logger_details, nday, num_today)
            n = file_vals[0]
            timesecs = file_vals[1]
            vals = file_vals[2]
            chids = file_vals[3]
            varnames = file_vals[4]
            qcflagnames = file_vals[5]
            files_to_write = file_vals[6]
            #sorted_indicies = file_vals[7]
            chosen_files = file_vals[7]

            print('No. values from previously written daily datalogger file = ', n)


    if n > 0:
        nvar = len(chids)
        nn_max = len(files_to_write)
        input_missing =  np.zeros((10000,nvar))
        valid_min_max = np.zeros((2,(nvar+1)))	#Making it nvar+1 is a fudge to make substitutions line work for writing the diagnostics T/RH file
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]
        #print('In main, sorted_indicies = ',sorted_indicies)

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        #print(timesecs[0:9])
        #print(vals[0:9,0])
        #print(vals[0:9,1])
        #print(vals[0:9,2])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #Something here to calculate body temperature and CG4 irradiance

        if sensor == 9:
            for m in range(len(varnames)):
                if varnames[m].find('downwelling_longwave_flux_in_air') >= 0:
                    m_longwave = m      #Column of vals which contains CG4 downwelling flux which needs to be corrected
            for m in range(len(varnames)):
                if varnames[m].find('body_temperature') >= 0 and chids[m].find('pyrCG4_ch') >= 0:
                    raw_v = vals[:,m]   #Column with voltage divider values that are used to derive temperature
                    R_CG4 = 10.0*(raw_v/(1-raw_v))      #Resistance in kohms
                    for k in range(n):
                        if R_CG4[k] <= 100.0:   #
                            vals[k,m] = (A*R_CG4[k]**5 + B*R_CG4[k]**4 + C*R_CG4[k]**3 + D*R_CG4[k]**2 + E*R_CG4[k] + F) + Tk
                            corr_CG4 = SBconst * vals[k,m] ** 4
                            vals[k,m_longwave] = vals[k,m_longwave] + corr_CG4
                        else:
                            vals[k,m_longwave] = missing_value  #This it should be missing as can't give it a value? Set a value for QC flag too?
                            vals[k,m] = missing_value   #This it should be missing as can't give it a value? Set a value for QC flag too?
                if varnames[m].find('body_temperature') >= 0 and chids[m].find('pyrCP1_T_ch') >= 0:
                    raw_v = vals[:,m]   #Column with voltage divider values that are used to derive temperature
                    R_CHP1 = 10000.0*(raw_v/(1-raw_v))  #Resistance in ohms
                    for k in range(n):
                        if R_CHP1[k] <= 100000.0:
                            lnR_CHP1 = np.log(R_CHP1);
                            vals[k,m] = 1.0/(A1 + (B1 * lnR_CHP1[k]) + (G1 * lnR_CHP1[k] ** 3))
                        else:
                            vals[k,m] = missing_value   #This it should be missing as can't give it a value? Set a value for QC flag too?

        #----------------------------------------------------------------------
        # Generate template and output filenames from input logger file details
        #----------------------------------------------------------------------
        #template_root = '/home/jla/python/metsensors/'
        template_root = '/home/chilbolton_software/python/ncas_python/metsensors/'
        template_suffix = '_metadata-template.yaml'
        tfp = [template_root + s + template_suffix for s in files_to_write] 
        #print('tfp = ', tfp)

        path_out_root = '/data/amof-netCDF/'
        path_out = [path_out_root+ s + '/' + s + '_cao_' for s in files_to_write] 
        if sensor == 5:
            path_out = [path_out_root+ s + '/' + s + '_cao-sparsholt_' for s in files_to_write] 
        print('path_out = ', path_out)
        graph_path_out_root = '/data/amof-netCDF/graphs/'
        graph_path_out = [graph_path_out_root+ s + data_version + '/'+ str(year) + '/' + strmon + '/' for s in files_to_write] 
        #print('graph_path_out = ', graph_path_out)


        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        #qc_flag = qc_flag + input_missing	#Add value of 1 from any values that were set as missing when reading data

        for k in range(nvar):
            col_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = col_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        #valid_min_max = corr_details[1]
        #print('qc_flag shape = ',qc_flag.shape)
        #print('valid_min_max = ',valid_min_max)
        #print('qc_flag shape = ',qc_flag.shape)

        raw_vals = vals	#Keep a copy of uncorrected values
        #Set to missing value if correction file showed BADDATA
        vals = np.where(qc_flag != 2, vals, missing_value)


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])
        #print(first_last_datetime)

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]


        for nn in range(nn_max):	#Each data file
            match_counter = 0
            for mm in range(nvar):      #Each variable

                if chosen_files[mm] == files_to_write[nn]:
                    if match_counter == 0:
                        vals_arr = [mm]
                    else:
                        vals_arr.append(mm)
                    match_counter = match_counter + 1

            print('vals_arr = ',vals_arr)
#            vals_arr = np.asarray(vals_arr, dtype = np.int)
#            print('vals_arr = ',vals_arr)
            how_many = len(vals_arr)
            #if how_many == 1:
            chids_arr = [chids[vals_arr[0]]]
            varnames_arr = [varnames[vals_arr[0]]]
            qcnames_arr = [qcflagnames[vals_arr[0]]]
            if how_many > 1:
                for k in range(how_many - 1):
                    chids_arr.append(chids[vals_arr[k+1]])
                    varnames_arr.append(varnames[vals_arr[k+1]])
                    qcnames_arr.append(qcflagnames[vals_arr[k+1]])

            print('chids_arr = ',chids_arr)
            print('varnames_arr = ',varnames_arr)
            print('qcnames_arr = ',qcnames_arr)

            sub_vals = np.squeeze(vals[:,vals_arr])
            print('sub_vals.shape = ',sub_vals.shape)
            print(sub_vals[0:9])
            sub_qcs = np.squeeze(qc_flag[:,vals_arr])
            #print('sub_qcs.shape = ',sub_qcs.shape)
            #print(sub_qcs[0:9])
            sub_min_max = valid_min_max[:,vals_arr]
            print('sub_min_max = ',sub_min_max)
            #print('sub_min_max.shape = ',sub_min_max.shape)

            if how_many == 1:
                varmin1 = sub_min_max[0] 
                varmax1 = sub_min_max[1]
                print('varmin1, varmax1 = ',varmin1, varmax1)
            else:
                varmin1 = sub_min_max[0,0] 
                varmax1 = sub_min_max[1,0]
                varmin2 = sub_min_max[0,1] 
                varmax2 = sub_min_max[1,1] 
                print('varmin1, varmax1, varmin2, varmax2 = ',varmin1, varmax1, varmin2, varmax2)


            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "valid_min_1": varmin1, "valid_max_1": varmax1, "valid_min_2": varmin2, "valid_max_2": varmax2, "uname": user_id, "codename": __file__, "machinename": computer_id}


            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print('Data_object_id = ',data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------

            #Derive data_product value, will need others or more robust method
            if 'rain-gauge' in files_to_write[nn]:
                data_product = "_precipitation"
            elif 'radiometer' in files_to_write[nn]:
                data_product = "_radiation"
            else:
                data_product = "_surface-met"

            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]


            for m in range(how_many):
                if how_many == 1:
                    print('varnames_arr, qcnames_arr = ', varnames_arr, qcnames_arr)
                    data_object["variables"][varnames_arr[m]]["values"][:] = sub_vals[:]
                    data_object["variables"][qcnames_arr[m]]["values"][:] = sub_qcs[:]
                else:
                    data_object["variables"][varnames_arr[m]]["values"][:] = sub_vals[:,m]
                    data_object["variables"][qcnames_arr[m]]["values"][:] = sub_qcs[:,m]

            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print('Data object handler exit code = ',exit_code)
            print('Changing permissions of file ',out_file)
            oscommand = "chgrp netcdf " + out_file 
            os.system(oscommand)
            #oscommand = "chmod g+w " + out_file
            oscommand = "chmod 775 " + out_file
            os.system(oscommand) 

            generate_netcdf_graphs(out_file, graph_path_out[nn], varnames_arr, datestring)


            #chids = ["oatnew_ch", "rhnew_ch"]
            #data_object["variables"]["air_temperature"]["values"][:] = raw_vals[:,0]
            #data_object["variables"]["relative_humidity"]["values"][:] = raw_vals[:,1]
            #data_object["variables"]["qc_flag_air_temperature"]["values"][:] = qc_flag[:,0]
            #data_object["variables"]["qc_flag_relative_humidity"]["values"][:] = qc_flag[:,1]

            #if nn < 8:	#Don't plot T/RH file with no corrections to data - it's just for quality control
                #generate_netcdf_graphs(out_file, graph_path_out[nn], varnames_arr, datestring)



# ------------------------------------------------------------------------
# Read data from Chilbolton1 Campbell datalogger
# ------------------------------------------------------------------------
def read_cr1000x_general(in_file, write_file, daily_file_read, daily_file_exist, logger_details, nday, num_today):

    # -------------------------
    # Define various parameters
    # -------------------------
    missing_value = -1.0E+20
    epoch_offset = 719163
    n = 0	#Times
    timesecs    = np.zeros((10000))
    #vals        = missing_value*np.ones((10000,nvar))
    write_daily_ok = 0	#Variable which shows whether we write a daily .dat file, 1 if we do write one
    wd_col = -999	#Error value for column which has wind direction data. Needs to be identified as needs extra filtering


    # ----------------------------
    # Read datalogger file details
    # ----------------------------

    #logger_details = '/data/netCDF/corrections/chilbolton_rxcabin_1.txt'
    chosen_chids = [' ']	#Dummy first element so that we can use the append command to add to the array
    chosen_files = [' ']	#Dummy first element so that we can use the append command to add to the array
    chosen_varnames = [' ']	#Dummy first element so that we can use the append command to add to the array
    chosen_qcflagnames = [' ']	#Dummy first element so that we can use the append command to add to the array
    f_logger = open(logger_details, 'r') 
    line = f_logger.readline()
    fdata = line.split(':')
    n_cols = int(fdata[1])
    on_off = np.zeros(n_cols)	#Flag showing whether each column in the data file is on or off. Some may have a fault or no instrument connected
    logger_cals = np.zeros((2,n_cols))
    line = f_logger.readline()
    line = line.strip('\n')
    fdata = line.split(':')
    first_column = int(fdata[1]) - 1
    for m in range(n_cols):
        line = f_logger.readline()
        line = line.strip('\n')
        fdata = line.split(',')
        print(fdata)
        if m == 0:
            chids = [fdata[0]]
        else:
            chids.append(fdata[0])
        if fdata[0] == 'wd_ch':
            wd_col = m
            #print('Wind direction column == ', wd_col)
        on_off[m] = int(fdata[1])
        if on_off[m] > 0:
            chosen_chids.append(fdata[0])
            chosen_files.append(fdata[4])
            chosen_varnames.append(fdata[5])
            chosen_qcflagnames.append(fdata[6])
        logger_cals[0,m] = float(fdata[2])
        logger_cals[1,m] = float(fdata[3])
    #print('n_cols, on_off = ', n_cols, on_off) 
    #print('logger_cals = ',logger_cals)

    f_logger.close()
    print('chids from text file = ',chids)

    chosen_chids.pop(0)
    chosen_files.pop(0)
    chosen_varnames.pop(0)
    chosen_qcflagnames.pop(0)
    print('chosen_chids where on_off == 1 = ',chosen_chids)
    print('chosen_varnames where on_off == 1 = ',chosen_varnames)
    #print('chosen_qcflagnames where on_off == 1 = ',chosen_qcflagnames)
    print('chosen_files where on_off == 1 = ',chosen_files)
    unique_chosen_files, indicies, inverses = np.unique(chosen_files, return_index=True, return_inverse=True)
    #sorted_indicies = np.sort(indicies)
    arg_sorted_indicies = np.argsort(indicies)
    #print('unique_chosen_files = ',unique_chosen_files)
    print('indicies = ',indicies)
    #print('sorted_indicies = ',sorted_indicies)
    print('arg_sorted_indicies = ',arg_sorted_indicies)
    ordered_chosen_files = np.array(unique_chosen_files)[arg_sorted_indicies] 
    #unique_chosen_files = chosen_files[sorted_indicies]
    print('ordered_chosen_files = ',ordered_chosen_files)
    nvar_local = len(chosen_chids)
    vals        = missing_value*np.ones((10000,nvar_local))

    f = open(in_file, 'r')
    print("In general read function, num_today, nday = ", num_today, nday)
    #In order to check whether we should write the CR1000X data line by line to a daily file check that
    #you're not reading data from a daily file, that the daily file doesn't already exist
    #that you're processing yesterday's data and that the time of day is later than 01UT
    # (so that the data have been mirrored to /data/range/mirror_grape_loggernet)
    # In addition, cronjobs to produce yesterday's data should be after ~ 01:15
    if daily_file_read == 0 and daily_file_exist == 0 and (num_today - nday) >= 1 and int(str(datetime.datetime.now())[11:13]) >= 1:
        write_daily_ok = 1
        fd = open(write_file, 'w')

    for z in range(4):      #4 header lines, may eventually want some of them
        line = f.readline()
        if daily_file_read == 0 and daily_file_exist == 0:	#Keep the header lines to write out later unless you're reading a daily file
            if z == 0:
                header_list = [line]
            else:
                header_list.append(line)

    # ---------------------------
    # Repeated data for each time
    # ---------------------------

    while True:

        line = f.readline()
        if not line: break
        column_count = 0	#Counter for each column in that line that is on
        #startnum = int(date2num(datetime.datetime(int(line[1:5]),int(line[6:8]),int(line[9:11]),0,0,0)))
        startnum = date2num(datetime.datetime(int(line[1:5]),int(line[6:8]),int(line[9:11]),int(line[12:14]),int(line[15:17]),int(line[18:20])))
        if startnum < 100000:
            startnum = startnum + epoch_offset
        if startnum > nday and startnum <= (nday + 1):
            #print(line)
            #if n == 0 and daily_file_read == 0 and daily_file_exist == 0:	#Write header to file if this is first line for today's date
            if n == 0 and write_daily_ok == 1:	#Write header to file if this is first line for today's date
                for z in range(4):
                    fd.write(header_list[z])	    

            timesecs[n] = 3600.0*float(line[12:14]) + 60.0*float(line[15:17]) + float(line[18:20])
            if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                timesecs[n] = timesecs[n] + 86400.0

            #print('timesecs = ',timesecs[n]) 
            fdata = line.split(',')
            #print(fdata)
            #If a data value is missing, the string has 5 characters, "NAN". Not all values will create such a long string, so just look for the first N
            for m in range(n_cols):	#All columns in raw datalogger file, whether on or off
                if on_off[m] > 0:	#Column is on
                    #print(fdata[m + first_column][1:2])
                    if fdata[m + first_column][1:2] != "N" and fdata[m + first_column][1:2] != "I":
                        temp_val = float(fdata[m + first_column])
                        if m != wd_col:
                            vals[n, column_count] = logger_cals[0, m] * temp_val + logger_cals[1, m]
                        else:
                            if np.absolute(temp_val) <= 5000.0:
                                temp_val = logger_cals[0, m] * temp_val + logger_cals[1, m]
                                if  temp_val >= 360.0:
                                    vals[n, column_count] = temp_val - 360.0
                                elif temp_val < 0.0:
                                    vals[n, column_count] = temp_val + 360.0
                                else:
                                    vals[n, column_count] = temp_val

                        #print('n,column_count,temp_val = ',n, column_count, temp_val)

                    column_count = column_count + 1
  
            #Write line to daily file
            if write_daily_ok == 1:
                fd.write(line)
            n += 1

    f.close()
    if write_daily_ok == 1:
        fd.close()	#Writes daily files
        print('Changing permissions of daily file ', write_file)
        oscommand = "chgrp netcdf " + write_file 
        os.system(oscommand)
        #oscommand = "chmod g+w " + out_file
        oscommand = "chmod 775 " + write_file
        os.system(oscommand) 
    
    print('No. values in read function = ', n)

    return n, timesecs, vals, chosen_chids, chosen_varnames, chosen_qcflagnames, ordered_chosen_files, chosen_files


# ------------------------------------------------------------------------
# Raingauge from chan* file netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_rain_f5(nday):

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['rg001dc_ch', 'rg006dc_ch', 'rg008dc_ch', 'rg004tb_ch']	#rain-gauge-1, rain-gauge-2, rain-gauge-3, rain-gauge-5
    nvar = np.size(chids)
    print('nvar = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    #Calibration factors
    if nday < 737599:	#Before 23/06/20
        rgcal = np.asarray([0.0033, 0.0033, 0.00189, 0.2]) 
    elif nday == 737599:	#On 23/06/20 - rg001, rg006 were calibrated previous day
        rgcal = np.asarray([0.00331, 0.00329, 0.00189, 0.2]) 
    else:	#After 23/06/20
        rgcal = np.asarray([0.00331, 0.00329, 0.00186, 0.2]) 

    missing_value = -1.0E+20
    epoch_offset = 719163

    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nvar))
    input_missing    = np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))


    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_marvin_home2/ranged/analog/"
    path_out_1 = "/data/amof-netCDF/ncas-rain-gauge-1/ncas-rain-gauge-1_cao_"
    path_out_2 = "/data/amof-netCDF/ncas-rain-gauge-2/ncas-rain-gauge-2_cao_"
    path_out_3 = "/data/amof-netCDF/ncas-rain-gauge-3/ncas-rain-gauge-3_cao_"
    path_out_4 = "/data/amof-netCDF/ncas-rain-gauge-5/ncas-rain-gauge-5_cao_"

    template_file_path_1 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-rain-gauge-1_ulink_metadata-template.yaml"	#YAML template
    template_file_path_2 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-rain-gauge-2_ulink_metadata-template.yaml"	#YAML template
    template_file_path_3 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-rain-gauge-3_ulink_metadata-template.yaml"	#YAML template
    template_file_path_4 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-rain-gauge-5_ulink_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1, template_file_path_2, template_file_path_3, template_file_path_4]
    path_out = [path_out_1, path_out_2, path_out_3, path_out_4]
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-rain-gauge-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_2 = "/data/amof-netCDF/graphs/ncas-rain-gauge-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_3 = "/data/amof-netCDF/graphs/ncas-rain-gauge-3" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_4 = "/data/amof-netCDF/graphs/ncas-rain-gauge-5" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1, graph_path_out_2, graph_path_out_3, graph_path_out_4]
    data_product = "_precipitation"

    # -------------------------------------------------------------------------
    # Loop to process raw data from previous, current and next day's data files
    # Data for the current day could be found in each of these format5 files.
    # contains data entirely from the previous day.
    # Hence, need to open files from all 3 days.
    # -------------------------------------------------------------------------

    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to process files')
    n = 0       #Times

    for day_incr in range(3):

        nday_file = nday + day_incr - 1
        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        print(nday_file, datevals)

        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "chan" + datestring_now[2:] + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print(infiles)
            nfiles = len(infiles)
            if nfiles == 0:
                print('No files on this day')

        else:
            nfiles = 0
            print('No directory found')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            for z in range(17):      #4 header lines, may eventually want some of them
                line = f.readline()

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #May need to change to handle change of year. Or maybe easier to just lose 1 point once a year?! Stems from not having year in the f5 file
                startnum = date2num(datetime.datetime(int(year_now),int(line[0:2]),int(line[3:5]),int(line[6:8]),int(line[9:11]),int(line[12:14])))
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum > nday and startnum <= (nday + 1):

                #startnum = int(date2num(datetime.datetime(int(datestring[0:4]),int(line[0:2]),int(line[3:5]),0,0,0)))
                #if startnum == nday:

                    timesecs[n] = 3600.0*float(line[6:8]) + 60.0*float(line[9:11]) + float(line[12:14])
                    if n == 0:
                        print('First point = ', nday, timesecs[0], startnum, datestring_now)
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0

                    fdata = line.split()	#Need to not specify ' ' or it treats sucessive delimiters as individuals
                    #vals[n,1] = stdrgcal * int(fdata[1])	#rg006dc_ch
                    #vals[n,2] = lowrgcal * int(fdata[2])	#rg008dc_ch
                    #vals[n,3] = tbrgcal * int(fdata[4])	#rg004tb_ch
                    #vals[n,0] = stdrgcal * int(fdata[15])	#rg001dc_ch

                    vals[n,1] = rgcal[1] * int(fdata[1])	#rg006dc_ch
                    vals[n,2] = rgcal[2] * int(fdata[2])	#rg008dc_ch
                    vals[n,3] = rgcal[3] * int(fdata[4])	#rg004tb_ch
                    vals[n,0] = rgcal[0] * int(fdata[15])	#rg001dc_ch
                    n += 1

            f.close()

        print('no. values  = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])
        print(vals[0:9,2])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs,vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 3
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        valid_min_max = corr_details[1]
        print('qc_flag shape = ',qc_flag.shape)

        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)	#Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 2, vals, missing_value)


        print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values 


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        #substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],sc_arr[0]), "end_time": datetime.datetime(year,month,day,hr_arr[n-1],mi_arr[n-1],sc_arr[n-1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "T_min": np.amin(vals[:,0]), "T_max": np.amax(vals[:,0]), "RH_min": np.amin(vals[:,1]), "RH_max": np.amax(vals[:,1]), "P_min": np.amin(vals[:,2]), "P_max": np.amax(vals[:,2]), "WS_min": np.amin(vals[:,3]), "WS_max": np.amax(vals[:,3]), "WD_min": np.amin(vals[:,4]), "WD_max": np.amax(vals[:,4]), "min_thick": np.amin(vals[:,5]), "max_thick": np.amax(vals[:,5]), "uname": user_id, "codename": __file__, "machinename": computer_id}


        for nn in range(nvar):	#Each data file

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "min_thick": valid_min_max[0,nn], "max_thick": valid_min_max[1,nn], "uname": user_id, "codename": __file__, "machinename": computer_id}

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]


            #chids = ["rg001dc_ch"]
            data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,nn] 
            data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,nn]
            variables_to_plot = ["thickness_of_rainfall_amount"]

            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print(exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)


# ------------------------------------------------------------------------
# RD-80 impact disdrometer chds* or spds* file netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_disdro_f5(nday,sensor):


    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    num_today = int(mpl.dates.date2num(datetime.datetime(int(datestring[0:4]),int(datestring[4:6]),int(datestring[6:8]),0,0,0)))
    print("In disdro, num_today = ", num_today)
    if num_today < 100000:
        epoch_offset = 719163 
    else:
        epoch_offset = 0
    #print "datestring = ",datestring
    #print(datevals)
    if sensor == 6:
        chids = ['disdrom_ch']  #Chilbolton disdrometer
        tempstr = 'chds' + datestring[2:] + '.00'	#Root of filename to search for in directory
    else:
        chids = ['disdrom_sp']	#Sparsholt disdrometer 
        tempstr = 'spds' + datestring[2:] + '.00'
    nvar = np.size(chids)
    print('No. variables to read = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    missing_value = -1.0E+20
    nbins = 127
    disdro_A=5.0265E-03 #m2
    #Minimum diameters for each bin
    rd80_thresholds=[0.313, 0.318, 0.324, 0.331,
                    0.337, 0.343, 0.35, 0.357, 0.364, 0.371,
                    0.379, 0.387, 0.396, 0.405, 0.414, 0.423,
                    0.433, 0.443, 0.453, 0.464, 0.474, 0.484,
                    0.495, 0.505, 0.515, 0.526, 0.537, 0.547,
                    0.558, 0.57, 0.583, 0.596, 0.611, 0.628,
                    0.644, 0.662, 0.679, 0.696, 0.715, 0.735,
                    0.754, 0.771, 0.787, 0.806, 0.827, 0.845,
                    0.862, 0.879, 0.895, 0.912, 0.928, 0.944,
                    0.96, 0.978, 0.999, 1.024, 1.051, 1.08,
                    1.11, 1.14, 1.171, 1.202, 1.232, 1.262,
                    1.289, 1.318, 1.346, 1.374, 1.402, 1.429,
                    1.456, 1.483, 1.509, 1.533, 1.558, 1.582,
                    1.606, 1.631, 1.657, 1.683, 1.715, 1.748,
                    1.793, 1.841, 1.897, 1.955, 2.013, 2.077,
                    2.139, 2.2, 2.262, 2.321, 2.381, 2.441,
                    2.499, 2.558, 2.616, 2.672, 2.727, 2.781,
                    2.836, 2.893, 2.949, 3.011, 3.08, 3.155,
                    3.23, 3.306, 3.385, 3.466, 3.545, 3.625,
                    3.704, 3.784, 3.864, 3.945, 4.028, 4.127,
                    4.231, 4.34, 4.456, 4.573, 4.686, 4.801,
                    4.915, 5.03, 5.145]



    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((30000))
    vals        = missing_value*np.ones((30000,nbins))	#Could/should expand to 10000,nbins,nvar
    input_missing    = np.zeros((30000,nvar))
    #valid_min_max = np.zeros((2))
    mean_diam    = np.zeros((127))
    mean_vol     = np.zeros((127))
    ds_bin_width = np.zeros((127))



    # ---------------------------------------------------
    # Calculate mean drop diameter, volume and bin widths 
    # ---------------------------------------------------
    for n in range(126):
        mean_diam[n]=np.mean([rd80_thresholds[n],rd80_thresholds[n+1]])
        mean_vol[n]=0.5236*np.mean([rd80_thresholds[n]**3,rd80_thresholds[n+1]**3])
        ds_bin_width[n]=rd80_thresholds[n+1]-rd80_thresholds[n]
    mean_diam[126] = 5.203      #arbitrary value based on spacing of previous bins
    mean_vol[126] = 74.0
    ds_bin_width[126] = ds_bin_width[125]       #arbitrary value based on spacing of previous bins

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    data_version = "_v1.0"
    if sensor == 6:
        path_in = "/data/range/mirror_marvin_home2/ranged/distrom/"
        if nday < 737916:	#20210506
            path_out_1 = "/data/amof-netCDF/ncas-disdrometer-1/ncas-disdrometer-1_cao_"
            template_file_path_1 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-disdrometer-1_metadata-template.yaml"	#YAML template
            graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-disdrometer-1" + data_version + "/" + str(year) + '/' + strmon + '/'
        else:
            path_out_1 = "/data/amof-netCDF/ncas-disdrometer-2/ncas-disdrometer-2_cao_"
            template_file_path_1 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-disdrometer-2-chilbolton_metadata-template.yaml"	#YAML template
            graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-disdrometer-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    else:	#Should only be 6 or 7 if you get to this function
        path_in = "/data/sparsholt/mirror_rhubarb_home2/ranged/distrom/"
        path_out_1 = "/data/amof-netCDF/ncas-disdrometer-2-sparsholt/ncas-disdrometer-2_cao-sparsholt_"
        graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-disdrometer-2-sparsholt" + data_version + "/" + str(year) + '/' + strmon + '/'
        template_file_path_1 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-disdrometer-2_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1]
    path_out = [path_out_1]
    graph_path_out = [graph_path_out_1]
    data_product = "_precipitation"

    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to process files')
    n = 0	#Times

    for day_incr in range(3):
        n_wrong = -1

        nday_file = nday + day_incr - 1

        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now

        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "chds" + datestring_now[2:] + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print('Raw files = ',infiles)
            nfiles = len(infiles)
            #print("nfiles = ",nfiles)
            if nfiles == 0:
                print('No files on this day')
            #Check for any files whose name start with . and remove them from list. A rare problem but it happened on 03/03/21
            if nfiles > 0:
                for m in range(nfiles):
                    if infiles[m][0] == ".":
                        n_wrong = m
                        print("File ", n_wrong, " found with wrong format on ", datestring_now)
                        nfiles = nfiles - 1
                        print("Number of files reduced to ", nfiles)
                if n_wrong >= 0:
                    del infiles[n_wrong]
                    print("Modified infiles = ",infiles)


        else:
            nfiles = 0
            print('No directory found')

        #for nf in range(nfiles):    #Indent from here
            #f = open(path_in+infiles[nf], 'r')




        for nf in range(nfiles):    #Indent from here

            source_file_path = path_in + infiles[nf]
            disdro_handler = module_distrometer_format5.Handler(source_file_path)

            number_of_records = disdro_handler.return_number_of_records()
            extracted_spectral_data = disdro_handler.return_spectral_data(0)       #0 is the index. An index is needed
            #For some reason, when I read values into vals, the value of extracted_spectral_data is set to zero in the print statement below.
            extracted_raw_data = disdro_handler.return_raw_data(0)
            datetimes = disdro_handler.return_datetimes_for_records()
            #print('datetimes = ',datetimes)
            numtimes = np.rint((date2num(datetimes) - nday + epoch_offset) * 86400.0)	#Time in seconds since midnight 
            print('numtimes = ',numtimes)
            print('nfiles, n, number_of_records, len(numtimes), (n+number_of_records) = ',nfiles, n, number_of_records, len(numtimes), (n+number_of_records))
            timesecs[n : (n + number_of_records)] = numtimes 

            for nn in range(number_of_records):
                vals[n,:] = disdro_handler.return_spectral_data(nn)
                n = n + 1	#Counter for total number of points across all files for the days



            print('Number of records = ', number_of_records)
            #print('Spectral data = ', extracted_spectral_data)
            #print('Raw data = ', extracted_raw_data)
            print('First element of datetimes array = ', datetimes[0])
            print('First numeric time = ', date2num(datetimes[0]))
            print('Numeric times = ', date2num(datetimes))
            print(timesecs)
            #print(vals[0,:])
            #print(vals[1,:])

            print('Number of records read = ', n)
            print('Last time read (seconds since midnight that day) = ',timesecs[n-1])

            #if timesecs[n-1] > 86400.0: 
                #n = n - 1

            



    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        #timesecs = timesecs[(timesecs > 0) & (timesecs <= 86400.0)]
        vals = vals[(timesecs > 0) & (timesecs <= 86400.0),:]
        timesecs = timesecs[(timesecs > 0) & (timesecs <= 86400.0)]
        n = len(timesecs)
        print('n after time filtering = ',n)
        #vals = vals[0:n,:]
        #print('vals.shape = ',vals.shape)
        #timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        #print(timesecs[0:9])
        #print(vals[0:9,0])
        #print(vals[0:9,1])
        #print(vals[0:9,2])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        ch_count_vol_product = mean_vol*vals
        #print('ch_count_vol_product min, max = ', np.amin(ch_count_vol_product), np.amax(ch_count_vol_product))
        ch_accum_array = (np.sum(ch_count_vol_product,axis=1))/(disdro_A*1.0e+06)
        #print( ch_accum_array)
        #print( 'ch_accum_array min, max = ', np.amin(ch_accum_array), np.amax(ch_accum_array))
        #print( ch_accum_array.shape)

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs,vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 3
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            #if len(vals_ok) > 0:
                #valid_min_max[0,k] = np.amin(vals_ok)
                #valid_min_max[1,k] = np.amax(vals_ok)

        valid_min_max = corr_details[1]
        #print('valid min, max = ',valid_min_max)
        valid_min_counts = int(valid_min_max[0])
        valid_max_counts = int(valid_min_max[1])
        #print('qc_flag shape = ',qc_flag.shape)
        qc_flag_thickness = qc_flag[:,0]	#1d qc_flag for thickness_of_rainfall_amount
        qc_flag = np.broadcast_to(qc_flag,(n,nbins))
        print('Maximum qc_flag, qc_flag_thickness = ', np.amax(qc_flag), np.amax(qc_flag_thickness))


        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)	#Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 3, vals, missing_value)


        #print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values 


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #lengths_of_dimensions = time_details[8]
        #substitutions = time_details[9]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        #substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],sc_arr[0]), "end_time": datetime.datetime(year,month,day,hr_arr[n-1],mi_arr[n-1],sc_arr[n-1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "T_min": np.amin(vals[:,0]), "T_max": np.amax(vals[:,0]), "RH_min": np.amin(vals[:,1]), "RH_max": np.amax(vals[:,1]), "P_min": np.amin(vals[:,2]), "P_max": np.amax(vals[:,2]), "WS_min": np.amin(vals[:,3]), "WS_max": np.amax(vals[:,3]), "WD_min": np.amin(vals[:,4]), "WD_max": np.amax(vals[:,4]), "min_thick": np.amin(vals[:,5]), "max_thick": np.amax(vals[:,5]), "uname": user_id, "codename": __file__, "machinename": computer_id}


        for nn in range(nvar):	#Each data file

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "min_diam": np.amin(rd80_thresholds), "max_diam": np.amax(rd80_thresholds), "min_counts": valid_min_counts, "max_counts": valid_max_counts, "min_thick": np.amin(ch_accum_array), "max_thick": np.amax(ch_accum_array), "uname": user_id, "codename": __file__, "machinename": computer_id}

#rt_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_las
#t_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_dat
#etime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]),

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n, "diameter": nbins}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]
            data_object["variables"]["diameter"]["values"][:] = rd80_thresholds[:]

            data_object["variables"]["number_of_hydrometeors_per_size_channel"]["values"][:] = vals[:,:] 
            data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = ch_accum_array[:] 
            data_object["variables"]["qc_flag_hydrometeors"]["values"][:] = qc_flag[:,:]
            data_object["variables"]["qc_flag_thickness_of_rainfall_amount"]["values"][:] = qc_flag_thickness[:]
            variables_to_plot = ["number_of_hydrometeors_per_size_channel", "thickness_of_rainfall_amount"]
            #variables_to_plot = ["thickness_of_rainfall_amount"]

            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print('Data object handler exit code = ',exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)

            print('Changing permissions of file ',out_file)
            oscommand = "chgrp netcdf " + out_file 
            os.system(oscommand)
            #oscommand = "chmod g+w " + out_file
            oscommand = "chmod 775 " + out_file
            os.system(oscommand) 

# ------------------------------------------------------------------------
# Broadband radiometer from chpy* file netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_bbrad_f5(nday):

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    chids = ['pyrCM21_ch','pyr_CMP21_ch', 'pyrCG4_ch', 'pyrCG4_ch', 'pyrCP1_ch', 'pyrCP1_T_ch']	#total pyranometer, total pyrgeo, pyrgeo T (duplicate same chids as no separate corr file), diffuse pyran, direct pyrhel, direct pyrhel T
    nvar = np.size(chids)
    print('Number of variables to be read = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    A = -4.56582625e-7
    B =  8.97922514e-5
    C = -6.95640241e-3
    D =  2.74163515e-1
    E = -6.23224724
    F = 66.1824897

    A1=0.0010295
    B1=2.391e-4
    G1=1.568e-7
    Tk = 273.15
    SBconst = 5.67e-8
    scal = load_bbrad_calibrations(nday)    
    #scal = np.zeros((nvar-1))
    #scal[0] = 10.89e-6	#CM21
    #scal[1] = 8.7e-6	#CMP21
    #scal[2] = 13.96e-6	#CG4
    #scal[3] = 7.89e-6	#CHP1 
    missing_value = -1.0E+20
    epoch_offset = 719163
    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nvar))
    input_missing    = np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_moe_home2/range/broadband_radiometers/pyr/"
    path_out_1 = "/data/amof-netCDF/ncas-radiometer-1/ncas-radiometer-1_cao_"
    path_out_2 = "/data/amof-netCDF/ncas-radiometer-2/ncas-radiometer-2_cao_"
    path_out_3 = "/data/amof-netCDF/ncas-radiometer-3/ncas-radiometer-3_cao_"
    path_out_4 = "/data/amof-netCDF/ncas-radiometer-4/ncas-radiometer-4_cao_"

    template_file_path_1 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-radiometer-1_agilent_metadata-template.yaml"	#YAML template
    template_file_path_2 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-radiometer-2_agilent_metadata-template.yaml"	#YAML template
    template_file_path_3 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-radiometer-3_agilent_metadata-template.yaml"	#YAML template
    template_file_path_4 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-radiometer-4_agilent_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1, template_file_path_2, template_file_path_3, template_file_path_4]
    path_out = [path_out_1, path_out_2, path_out_3, path_out_4]
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-radiometer-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_2 = "/data/amof-netCDF/graphs/ncas-radiometer-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_3 = "/data/amof-netCDF/graphs/ncas-radiometer-3" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_4 = "/data/amof-netCDF/graphs/ncas-radiometer-4" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1, graph_path_out_2, graph_path_out_3, graph_path_out_4]
    data_product = "_radiation"

    # -------------------------------------------------------------------------
    # Loop to process raw data from previous, current and next day's data files
    # Data for the current day could be found in each of these format5 files.
    # contains data entirely from the previous day.
    # Hence, need to open files from all 3 days.
    # -------------------------------------------------------------------------

    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to read raw format5 files')
    n = 0       #Times

    for day_incr in range(3):
        n_wrong = -1

        nday_file = nday + day_incr - 1

        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now

        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "chpy" + datestring_now[2:] + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print('Raw files = ',infiles)
            nfiles = len(infiles)
            #print("nfiles = ",nfiles)
            if nfiles == 0:
                print('No files on this day')
            #Check for any files whose name start with . and remove them from list. A rare problem but it happened on 03/03/21
            if nfiles > 0:
                for m in range(nfiles):
                    if infiles[m][0] == ".":
                        n_wrong = m
                        print("File ", n_wrong, " found with wrong format on ", datestring_now)
                        nfiles = nfiles - 1
                        print("Number of files reduced to ", nfiles)
                if n_wrong >= 0:
                    del infiles[n_wrong]
                    print("Modified infiles = ",infiles)


        else:
            nfiles = 0
            print('No directory found')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            for z in range(8):      #8 header lines, may eventually want some of them
                line = f.readline()

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:
 
                line = f.readline()
                if not line: break
                #startnum = int(date2num(datetime.datetime(int(datestring[0:4]),int(line[0:2]),int(line[3:5]),0,0,0)))
                #print 'startnum, nday = ',startnum, nday
                #if startnum == nday:
                #May need to change to handle change of year. Or maybe easier to just lose 1 point once a year?!
                startnum = date2num(datetime.datetime(int(year_now),int(line[0:2]),int(line[3:5]),int(line[6:8]),int(line[9:11]),int(line[12:14])))
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum > nday and startnum <= (nday + 1):

                    timesecs[n] = 3600.0*float(line[6:8]) + 60.0*float(line[9:11]) + float(line[12:14])
                    if n == 0:
                        print('First point = ', nday, timesecs[0], startnum, datestring_now)
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0

                    fdata = line.split()	#Need to not specify ' ' or it treats sucessive delimiters as individuals
                    vals[n,0] = float(fdata[1])/scal[0]	#CM21
                    vals[n,1] = float(fdata[4])/scal[1]	#CMP21
                    vals[n,2] = float(fdata[2])/scal[2]	#CG4
                    vals[n,4] = float(fdata[5])/scal[3]	#CHP1
                    #Temperatures and corrections
                    R_CG4 = float(fdata[3])/1000.	#kohms. 10^5 ohms is sensible highest
                    if R_CG4 <= 100.0:
                        T_CG4 = (A*R_CG4**5 + B*R_CG4**4 + C*R_CG4**3 + D*R_CG4**2 + E*R_CG4 + F) + Tk
                        corr_CG4 = SBconst * T_CG4 ** 4               
                        vals[n,2] = vals[n,2] + corr_CG4
                        vals[n,3] = T_CG4
                    else:
                        vals[n,2] = missing_value	#This it should be missing as can't give it a value? Set a value for QC flag too?
                        vals[n,3] = missing_value
                    R_CHP1 = float(fdata[6])
                    if R_CHP1 <= 100000.0:
                        lnR_CHP1 = np.log(R_CHP1);
                        vals[n,5] = 1.0/(A1 + (B1 * lnR_CHP1) + (G1 * lnR_CHP1 ** 3)) 
                    else:
                        vals[n,5] = missing_value

                    n += 1

            f.close()
    
        print('No. values from raw format5 files = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        #print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        #valid_min_max = corr_details[1]
        #print('qc_flag shape = ',qc_flag.shape)

        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)	#Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 2, vals, missing_value)


        #print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values 


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "CM21_min": valid_min_max[0,0], "CM21_max": valid_min_max[1,0], "CMP21_min": valid_min_max[0,1], "CMP21_max": valid_min_max[1,1], "CG4_min": valid_min_max[0,2], "CG4_max": valid_min_max[1,2], "CG4_T_min": valid_min_max[0,3], "CG4_T_max": valid_min_max[1,3], "CHP1_min": valid_min_max[0,4], "CHP1_max": valid_min_max[1,4], "CHP1_T_min": valid_min_max[0,5], "CHP1_T_max": valid_min_max[1,5], "uname": user_id, "codename": __file__, "machinename": computer_id}


        for nn in range(nvar-2):	#Each data file

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print('Data_object_id = ',data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]


            if nn == 0:
                chids = ["pyrCM21_ch"]
                #Call corrections, but may want to re-order so only add that code once
                data_object["variables"]["downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,0]
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,0]
                variables_to_plot = ["downwelling_shortwave_flux_in_air"]
            if nn == 1:
                chids = ["pyr_CMP21_ch"]
                data_object["variables"]["diffuse_downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,1]	#Should diffuse be added to name?
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,1]
                variables_to_plot = ["diffuse_downwelling_shortwave_flux_in_air"]
            if nn == 2:
                chids = ["pyrCG4_ch, pyrCG4_ch"]
                data_object["variables"]["downwelling_longwave_flux_in_air"]["values"][:] = vals[:,2] 
                data_object["variables"]["body_temperature"]["values"][:] = vals[:,3] 
                data_object["variables"]["qc_flag_downwelling_longwave_flux_in_air"]["values"][:] = qc_flag[:,2]
                data_object["variables"]["qc_flag_body_temperature"]["values"][:] = qc_flag[:,3]
                variables_to_plot = ["downwelling_longwave_flux_in_air", "body_temperature"]
            if nn == 3:
                chids = ["pyrCP1_ch, pyrCP1_T_ch"]
                data_object["variables"]["direct_downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,4] 
                data_object["variables"]["body_temperature"]["values"][:] = vals[:,5] 
                data_object["variables"]["qc_flag_direct_downwelling_shortwave_flux_in_air"]["values"][:] = qc_flag[:,4]
                data_object["variables"]["qc_flag_body_temperature"]["values"][:] = qc_flag[:,5]
                variables_to_plot = ["direct_downwelling_shortwave_flux_in_air", "body_temperature"]
            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print('Data object handler exit code = ',exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)

            print('Changing permissions of file ',out_file)
            oscommand = "chgrp netcdf " + out_file 
            os.system(oscommand)
            #oscommand = "chmod g+w " + out_file
            oscommand = "chmod 775 " + out_file
            os.system(oscommand) 

# ------------------------------------------------------------------------
# Campbell PWS100 present weather sensor netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_pws100(nday):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon

    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals)
    chids = ['pws100','pws100','pws100_met','pws100_met','pws100_met','pws100_met','pws100_met','pws100_met']   #Need a file for each dataset read, but can be same one.
    nvar = np.size(chids)
    n = 0       #Number of data points

    # -------------------------
    # Define various parameters
    # -------------------------

    missing_value = -1.0E+20
    epoch_offset = 719163

    # ----------------------
    # Initialize data arrays
    # ----------------------
    nbins = 300 #Number of size bins
    vals        = missing_value * np.ones((10000,nvar))
    vals[:,1] = -127.0	#Missing value for synoptic code, which is a short integer
    drop_vals        = missing_value * np.ones((10000,300))
    timesecs    = np.zeros((10000))
    input_missing =  np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))
    bin_middle = 0.1 * np.arange(nbins) + 0.05

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/campbell_PWS/mirror_cl51sky_campbell_data/PWS100/"
    path_out = "/data/amof-netCDF/ncas-present-weather-1/ncas-present-weather-1_cao_"
    data_version = "_v1.0"      #Changed temporarily so as to show this was generated with v3 python
    graph_path_out = "/data/amof-netCDF/graphs/ncas-present-weather-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    data_product = "_present-weather"

    template_file_path = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-present-weather-1_metadata-template.yaml"      #YAML template

    #------------Setting up YAML data object----------------------
    #Need export PYTHONPATH=/home/jla/python/global_modules
    handler = module_data_object_python3.Handler()
    data_object_id = handler.load_template(template_file_path)
    print('Data_object_id = ',data_object_id)
    #handler.show_requirements_for_template(data_object_id)

    # ---------------------------------------------------------------
    # Loop to process raw data from current and next day's data files
    # Since first data point in file is at 00:00:00, first data point
    # contains data entirely from the previous day.
    # Hence, need to open file from next day too.
    # ---------------------------------------------------------------

    for day_incr in range(1):   #Used in other codes, should only need to read current day for this sensor

        nday_file = nday + day_incr


        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        #print(nday_file, datevals)


        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "pws100-" + datestring_now + '-'
            expr = re.compile(tempstr+'[0-9][0-9][0-9][0-9][0-9][0-9]'+'.txt')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print('Raw files = ',infiles)
            nfiles = len(infiles)
            if nfiles == 0:
                print('No files on this day')

        else:
            nfiles = 0
            print('No directory found')

        # ----------------------
        # Process data files
        # ----------------------
        print('Starting to process raw data files')

        for nf in range(nfiles):
            f = open(path_in+infiles[nf], 'r')

            #There are no header lines

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #print fdata
                fdata = line.split()
                #print(len(line))
                startnum = date2num(datetime.datetime(int(year_now),int(month_now),int(day_now),int(line[0:2]),int(line[3:5]),int(line[6:8])))
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum > nday and startnum <= (nday + 1):
                    timesecs[n] = 3600.0*float(line[0:2]) + 60.0*float(line[3:5]) + float(line[6:8])
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0
                    if len(line) > 3100:        # There's data in the line beyond the quantities we're writing to a file
                        #fdata=fdatatmp[1].split(';')
                        vals[n,0] = float(fdata[3])     #Visibility
                        vals[n,1] = float(fdata[4])     #WMO Synop code, change to int later
                        vals[n,2] = float(fdata[24]) + 273.15   #Temperature (K)
                        vals[n,3] = float(fdata[25])    #RH (%)
                        vals[n,4] = float(fdata[29])    #Rain rate (mm/hr)
                        vals[n,5] = float(fdata[30])    #Rain accumulation (mm)
                        vals[n,6] = float(fdata[331])   #Average velocity (m/s)
                        vals[n,7] = float(fdata[332])   #Average size (mm)
                        drop_vals[n,:] = np.array(fdata[31:331], dtype = np.float32)    #Drops in each bin
                    n += 1

            f.close()

        print('No. values from raw data files = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        drop_vals = drop_vals[0:n,:]
        #input_missing = input_missing[0:n,:]
        timesecs   = timesecs[0:n]

        #print(timesecs[0:9])
        #print(vals[0:9,0])
        #print(vals[0:9,1])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]    #Assumes evenly spaced, may want to improve this

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        #Make any qc_flag values which were missing at input == 2. These won't have been read from the correction file
        for k in range(nvar):
            sub_vals = vals[:,k]        #Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)
        #qc_flag = qc_flag + input_missing      #Add value of 1 from any values that were set as missing when reading data
        #Can't use valid_min_max from function as won't have the input_missing values
        #valid_min_max = corr_details[1]


        #print('qc_flag shape = ',qc_flag.shape)
        #print('valid_min_max = ',valid_min_max)

        #vals_qc = np.where(qc_flag <> 2, vals, 0)
        #vals_qc = np.where(qc_flag <> 3, vals_qc, missing_value)
        #vals = np.where(qc_flag <> 3, vals, missing_value)
        vals = np.where(qc_flag != 2, vals, missing_value)

        #qc_flag_drops = np.broadcast_to(qc_flag[:,5],(n,nbins))        #Broadcast qc_flag for accumulation to drop counts as it's likely to be the sa
        #qc_flag_drops = np.broadcast_to(qc_flag[:,5],(nbins,n))        #Broadcast qc_flag for accumulation to drop counts as it's likely to be the sa
        qc_flag_drops = np.zeros((n,nbins))
        #qc_flag_drops[:,0:nbins] = np.squeeze(qc_flag[:,5])

        #print(qc_flag[:,5])
        #print(qc_flag[:,5].shape)
        #print(qc_flag_drops.shape)

        for m in range(n):
            qc_flag_drops[m,:] = qc_flag[m,5]



        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)       #Will over-ride HOLDCAL=0 set in previous line
        #PWS100 is not prone to false positive counts, so if a HOLDCAL period is encountered in the correction file
        #it is treated in the same way as BADDATA and flagged as missing
        vals = np.where(qc_flag != 3, vals, missing_value)
        drop_vals = np.where(qc_flag_drops != 3, drop_vals, missing_value)
        valid_min_counts = int(np.amin(drop_vals))
        valid_max_counts = int(np.amax(drop_vals))
        #print(qc_flag[0:49,0])


    #--------------------------------------------
    # Sorting out day/time information
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)
        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #lengths_of_dimensions = time_details[8]
        #substitutions = time_details[9]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])
        #print(first_last_datetime)

        #print('vals shape = ',vals.shape)

        # -------------------------------
        # Open new netCDF file for output
        # -------------------------------
        cfarr_head = 'ncas-present-weather-1_cao_'
        out_file   = path_out+datestring+data_product+data_version+'.nc'
        print('Opening new NetCDF file ' + out_file)

        lengths_of_dimensions = {"time": n, "diameter": nbins}

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec": np.amax(sc_arr),"vis_min": valid_min_max[0,0], "vis_max": valid_min_max[1,0], "synop_min": valid_min_max[0,1], "synop_max": valid_min_max[1,1], "T_min": valid_min_max[0,2], "T_max": valid_min_max[1,2], "RH_min": valid_min_max[0,3], "RH_max": valid_min_max[1,3], "min_rr": valid_min_max[0,4], "max_rr": valid_min_max[1,4], "min_thick": valid_min_max[0,5], "max_thick": valid_min_max[1,5], "min_diam": np.amin(bin_middle), "max_diam": np.amax(bin_middle), "min_counts": valid_min_counts, "max_counts": valid_max_counts, "uname": user_id, "codename": __file__, "machinename": computer_id}

        data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
        data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
        data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]      #Need to include fraction of day
        data_object["variables"]["year"]["values"][:] = yr_arr[:]
        data_object["variables"]["month"]["values"][:] = mn_arr[:]
        data_object["variables"]["day"]["values"][:] = dy_arr[:]
        data_object["variables"]["hour"]["values"][:] = hr_arr[:]
        data_object["variables"]["minute"]["values"][:] = mi_arr[:]
        data_object["variables"]["second"]["values"][:] = sc_arr[:]
        data_object["variables"]["diameter"]["values"][:] = bin_middle[:]
        data_object["variables"]["visibility_in_air"]["values"][:] = vals[:,0]
        data_object["variables"]["wmo_synoptic_code"]["values"][:] = vals[:,1]
        data_object["variables"]["air_temperature"]["values"][:] = vals[:,2]
        data_object["variables"]["relative_humidity"]["values"][:] = vals[:,3]
        data_object["variables"]["rainfall_rate"]["values"][:] = vals[:,4]
        data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,5]
        data_object["variables"]["number_of_hydrometeors_per_size_channel"]["values"][:] = drop_vals[:,:]
        data_object["variables"]["qc_flag_visibility_in_air"]["values"][:] = qc_flag[:,0]
        data_object["variables"]["qc_flag_wmo_synoptic_code"]["values"][:] = qc_flag[:,1]
        data_object["variables"]["qc_flag_air_temperature"]["values"][:] = qc_flag[:,2]
        data_object["variables"]["qc_flag_relative_humidity"]["values"][:] = qc_flag[:,3]
        data_object["variables"]["qc_flag_rainfall_rate"]["values"][:] = qc_flag[:,4]
        data_object["variables"]["qc_flag_thickness_of_rainfall_amount"]["values"][:] = qc_flag[:,5]
        #number_of_hydrometeors_per_fallspeed_channel_per_size_channel, different velocity/size grid
        data_object["variables"]["qc_flag_hydrometeors"]["values"][:] = qc_flag_drops[:,:]
        exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
        print(exit_code)



        variables_to_plot = ["visibility_in_air", "wmo_synoptic_code", "air_temperature", "relative_humidity", "rainfall_rate", "thickness_of_rainfall_amount", "number_of_hydrometeors_per_size_channel"]
        generate_netcdf_graphs(out_file, graph_path_out, variables_to_plot, datestring)

        print('Changing permissions of file ',out_file)
        oscommand = "chgrp netcdf " + out_file 
        os.system(oscommand)
        oscommand = "chmod g+w " + out_file
        os.system(oscommand)


# -----------------------------------------------------------------------------------
# CNR4 net flux radiometer from datataker in flux compound netcdf generation function
# -----------------------------------------------------------------------------------
def generate_netcdf_cnr4_netflux(nday):

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    chids = ['cnr4_dsw','cnr4_dlw', 'cnr4_usw', 'cnr4_ulw', 'cnr4_T']   #downwelling shortwave, downwelling longwave, upwelling shortwave, upwelling longwave
    nvar = np.size(chids)
    print('Number of variables to be read = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    Tk = 273.15
    SBconst = 5.67e-8
    A = 0.003908
    B = -5.8019E-07  #A and B are coefficents to calculate T(K) from Pt100 resistance

    #scal = load_bbrad_calibrations(nday)
    scal = np.zeros((nvar-1))
    scal[0] = 17.41e-3  #Downwelling shortwave, values in raw file are in mV
    scal[1] = 8.99e-3   #Downwelling longwave
    scal[2] = 13.97e-3  #Upwelling shortwave
    scal[3] = 7.43e-3   #Upwelling longwave
    missing_value = -1.0E+20
    epoch_offset = 719163
    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((100000))
    vals        = missing_value*np.ones((100000,nvar))
    input_missing    = np.zeros((100000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range-tower/mirror_range-tower_data/net/"
    path_out_1 = "/data/amof-netCDF/ncas-radiometer-5/ncas-radiometer-5_cao_"

    template_file_path_1 = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-radiometer-5_metadata-template.yaml"   #YAML template

    tfp = [template_file_path_1]
    path_out = [path_out_1]
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-radiometer-5" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1]
    data_product = "_radiation"

    # -------------------------------------------------------------------------
    # Loop to process raw data from previous, current and next day's data files
    # Data for the current day could be found in each of these format5 files.
    # contains data entirely from the previous day.
    # Hence, need to open files from all 3 days.
    # -------------------------------------------------------------------------


    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to read flux compound Datataker raw files')
    n = 0       #Times

    for day_incr in range(1):   #May want to expand to read day before and day after to check there are no data in those files
        n_wrong = -1

        #nday_file = nday + day_incr - 1
        nday_file = nday

        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now

        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "net_" + datestring_now
            expr = re.compile(tempstr+'....'+'.dat')    #'....' signifies 4 characters
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print('Raw files = ',infiles)
            nfiles = len(infiles)
            #print("nfiles = ",nfiles)
            if nfiles == 0:
                print('No files on this day')
            #Check for any files whose name start with . and remove them from list. A rare problem but it happened on 03/03/21
            if nfiles > 0:
                for m in range(nfiles):
                    if infiles[m][0] == "." or os.stat(path_in+infiles[m]).st_size == 0:
                        n_wrong = m
                        print("File ", n_wrong, " found with wrong format or zero size on ", datestring_now)
                        nfiles = nfiles - 1
                        print("Number of files reduced to ", nfiles)
                if n_wrong >= 0:
                    del infiles[n_wrong]
                    print("Modified infiles = ",infiles)


        else:
            nfiles = 0
            print('No directory found')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            #No header lines in raw file

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #startnum = int(date2num(datetime.datetime(int(datestring[0:4]),int(line[0:2]),int(line[3:5]),0,0,0)))
                #print 'startnum, nday = ',startnum, nday
                #if startnum == nday:
                #May need to change to handle change of year. Or maybe easier to just lose 1 point once a year?!
                startnum = date2num(datetime.datetime(int(year_now),int(month_now),int(day_now),int(line[9:11]),int(line[12:14]),int(line[15:17])))
                #print('startnum, nday = ',startnum, nday)
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum > nday and startnum <= (nday + 1):
                    if len(line) > 175:
                        fdata = line.split(',')
                        date_time = fdata[0].split()
                        timesecs[n] = 3600.0 * float(date_time[1][0:2]) + 60.0 * float(date_time[1][3:5]) + float(date_time[1][6:len(date_time[1])])
                        tempstr = fdata[19].split(';')
                        res = float(tempstr[0])/100.0    #Pt-100 resistance, div. by 100 to use in temperature equation
                        sensor_T = 273.15 + (((-1.0 * A) + np.sqrt(A**2 - (4.0 * B * (1.0 - res))))/(2.0 * B))
                        vals[n,4] = sensor_T
                        for m in range(4):
                            vals[n,m] = float(fdata[15 + m])/scal[m]
                            if m%2.0 == 1.0:
                                vals[n,m] = vals[n,m] + (SBconst * sensor_T**4)

                        if n == 0:
                            print('First point = ', nday, timesecs[0], startnum, datestring_now)
                        if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                            timesecs[n] = timesecs[n] + 86400.0

                        n += 1

            f.close()

        print('No. values from raw net flux files = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        #print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]        #Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        #valid_min_max = corr_details[1]
        #print('qc_flag shape = ',qc_flag.shape)

        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)       #Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 2, vals, missing_value)


        #print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values


    #--------------------------------------------
    # Sorting out day/time information
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "valid_min_1": valid_min_max[0,0], "valid_max_1": valid_min_max[1,0], "valid_min_2": valid_min_max[0,1], "valid_max_2": valid_min_max[1,1], "valid_min_3": valid_min_max[0,2], "valid_max_3": valid_min_max[1,2], "valid_min_4": valid_min_max[0,3], "valid_max_4": valid_min_max[1,3], "valid_min_5": valid_min_max[0,4], "valid_max_5": valid_min_max[1,4], "uname": user_id, "codename": __file__, "machinename": computer_id}

        for nn in range(1):     #Each data file

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print('Data_object_id = ',data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]  #Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]


            if nn == 0:
                chids = ["pyrCM21_ch"]
                #Call corrections, but may want to re-order so only add that code once
                data_object["variables"]["downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,0]
                data_object["variables"]["qc_flag_downwelling_shortwave_flux_in_air"]["values"][:] = qc_flag[:,0]
                data_object["variables"]["downwelling_longwave_flux_in_air"]["values"][:] = vals[:,1]
                data_object["variables"]["qc_flag_downwelling_longwave_flux_in_air"]["values"][:] = qc_flag[:,1]
                data_object["variables"]["upwelling_shortwave_flux_in_air"]["values"][:] = vals[:,2]
                data_object["variables"]["qc_flag_upwelling_shortwave_flux_in_air"]["values"][:] = qc_flag[:,2]
                data_object["variables"]["upwelling_longwave_flux_in_air"]["values"][:] = vals[:,3]
                data_object["variables"]["qc_flag_upwelling_longwave_flux_in_air"]["values"][:] = qc_flag[:,3]
                data_object["variables"]["body_temperature"]["values"][:] = vals[:,4]
                data_object["variables"]["qc_flag_body_temperature"]["values"][:] = qc_flag[:,4]
                variables_to_plot = ["downwelling_shortwave_flux_in_air", "downwelling_longwave_flux_in_air", "upwelling_shortwave_flux_in_air", "upwelling_longwave_flux_in_air", "body_temperature"]
            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print('Data object handler exit code = ',exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)

            print('Changing permissions of file ',out_file)
            oscommand = "chgrp netcdf " + out_file
            os.system(oscommand)
            #oscommand = "chmod g+w " + out_file
            oscommand = "chmod 775 " + out_file
            os.system(oscommand)


# ------------------------------------------------------------------------
# Metek USA-1 sonic anemometer netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_sonic(nday):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon

    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals)
    #Variables are northward wind, eastward wind, upward air velocity and air temperature
    chids = ['sonic','sonic','sonic','sonic']   #Need a file for each dataset read, but can be same one.
    nvar = np.size(chids)
    n = 0       #Number of data points

    # -------------------------
    # Define various parameters
    # -------------------------

    missing_value = -1.0E+20
    epoch_offset = 719163

    # ----------------------
    # Initialize data arrays
    # ----------------------
    vals        = missing_value * np.ones((1800000,nvar))
    timesecs    = np.zeros((1800000))
    input_missing =  np.zeros((1800000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range-tower/mirror_range-tower_data/usa/"
    #path_in = "/home/jla/range-tower/"
    path_out = "/data/amof-netCDF/ncas-sonic-5/ncas-sonic-5_cao_"
    data_version = "_v1.0"      #Changed temporarily so as to show this was generated with v3 python
    graph_path_out = "/data/amof-netCDF/graphs/ncas-sonic-5" + data_version + "/" + str(year) + '/' + strmon + '/'
    data_product = "_mean-winds"
    #path_out = "/home/jla/netCDF/files/amof-pluvio_chilbolton/"
    #text_path_out = "/mnt/wilma_home/jla/halo/text_profile/"

    #template_file_path = "/home/jla/python/metsensors/ncas-sonic-5_metadata-template.yaml"      #YAML template
    template_file_path = "/home/chilbolton_software/python/ncas_python/metsensors/ncas-sonic-5_metadata-template.yaml"      #YAML template

    #------------Setting up YAML data object----------------------
    #Need export PYTHONPATH=/home/jla/python/global_modules
    handler = module_data_object_python3.Handler()
    data_object_id = handler.load_template(template_file_path)
    print('Data_object_id = ',data_object_id)
    #handler.show_requirements_for_template(data_object_id)

    # ---------------------------------------------------------------
    # Loop to process raw data from current and next day's data files
    # Since first data point in file is at 00:00:00, first data point
    # contains data entirely from the previous day.
    # Hence, need to open file from next day too.
    # ---------------------------------------------------------------

    for day_incr in range(1):   #Update to more days to get data from previous or subsequent day's file

        nday_file = nday + day_incr


        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        #print(nday_file, datevals)


        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "usa_" + datestring_now
            expr = re.compile(tempstr+'[0-9][0-9][0-9][0-9]'+'.dat')
            infiles = [elem for elem in files if expr.search(elem)]

            infiles.sort()
            print('Raw files = ',infiles)
            nfiles = len(infiles)
            if nfiles == 0:
                print('No files on this day')

        else:
            nfiles = 0
            print('No directory found')

        # ----------------------
        # Process data files
        # ----------------------
        print('Starting to process raw data files')

        for nf in range(nfiles):
            f = open(path_in+infiles[nf], 'r')

            #There are no header lines

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #print fdata
                fdata = line.split()
                #print(len(line))
                startnum = date2num(datetime.datetime(int(year_now),int(month_now),int(day_now),int(line[9:11]),int(line[12:14]),int(line[15:17])))
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum > nday and startnum <= (nday + 1):
                    timesecs[n] = 3600.0*float(line[9:11]) + 60.0*float(line[12:14]) + float(line[15:21])
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0
                    if len(line) > 66:        # There's data in the line beyond the quantities we're writing to a file
                        #fdata=fdatatmp[1].split(';')
                        #print(line[9:21], ' ',line[24:32])
                        if line[24:32] != '        ' and line[35:43] != '        ' and line[46:54] != '        ' and line[57:65] != '        ':
                            vals[n,0] = float(line[24:32])     #Northward wind
                            vals[n,1] = float(line[35:43])     #Eastward wind
                            vals[n,2] = float(line[46:54])         #Upward
                            vals[n,3] = float(line[57:65]) + 273.15   #Temperature (K)

                        else:
                            input_missing[n,:] = 1
                    n += 1

            f.close()

        print('No. values from raw data files = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        input_missing = input_missing[0:n,:]
        timesecs   = timesecs[0:n]

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]    #Assumes evenly spaced, may want to improve this

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        #Make any qc_flag values which were missing at input == 2. These won't have been read from the correction file
        for k in range(nvar):
            sub_vals = vals[:,k]        #Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)
        #qc_flag = qc_flag + input_missing      #Add value of 1 from any values that were set as missing when reading data
        #Can't use valid_min_max from function as won't have the input_missing values
        #valid_min_max = corr_details[1]


        #print('qc_flag shape = ',qc_flag.shape)
        #print('valid_min_max = ',valid_min_max)

        vals = np.where(qc_flag != 2, vals, missing_value)

        #print(qc_flag[0:49,0])


    #--------------------------------------------
    # Sorting out day/time information
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)
        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #lengths_of_dimensions = time_details[8]
        #substitutions = time_details[9]
        #print(yr_arr[0:9])
        #print(mn_arr[0:9])
        #print(dy_arr[0:9])
        #print(epoch_timesecs[0:9])
        #print(hr_arr[0:9])
        #print(mi_arr[0:9])
        #print(sc_arr[0:9])
        #print(dyyr_arr[0:9])
        #print(first_last_datetime)

        #print('vals shape = ',vals.shape)

        # -------------------------------
        # Open new netCDF file for output
        # -------------------------------
        cfarr_head = 'ncas-present-weather-1_cao_'
        out_file   = path_out+datestring+data_product+data_version+'.nc'
        print('Opening new NetCDF file ' + out_file)

        lengths_of_dimensions = {"time": n}

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec": np.amax(sc_arr),"north_min": valid_min_max[0,0], "north_max": valid_min_max[1,0], "east_min": valid_min_max[0,1], "east_max": valid_min_max[1,1], "up_min": valid_min_max[0,2], "up_max": valid_min_max[1,2], "T_min": valid_min_max[0,3], "T_max": valid_min_max[1,3], "uname": user_id, "codename": __file__, "machinename": computer_id}


        data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
        data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
        data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]      #Need to include fraction of day
        data_object["variables"]["year"]["values"][:] = yr_arr[:]
        data_object["variables"]["month"]["values"][:] = mn_arr[:]
        data_object["variables"]["day"]["values"][:] = dy_arr[:]
        data_object["variables"]["hour"]["values"][:] = hr_arr[:]
        data_object["variables"]["minute"]["values"][:] = mi_arr[:]
        data_object["variables"]["second"]["values"][:] = sc_arr[:]
        data_object["variables"]["northward_wind"]["values"][:] = vals[:,0]
        data_object["variables"]["eastward_wind"]["values"][:] = vals[:,1]
        data_object["variables"]["upward_air_velocity"]["values"][:] = vals[:,2]
        data_object["variables"]["air_temperature"]["values"][:] = vals[:,3]
        data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,0]
        exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
        print(exit_code)


        variables_to_plot = ["northward_wind", "eastward_wind", "upward_air_velocity", "air_temperature"]
        generate_netcdf_graphs(out_file, graph_path_out, variables_to_plot, datestring)

        print('Changing permissions of file ',out_file)
        oscommand = "chmod g+w " + out_file
        os.system(oscommand)




# ------------------------------------------------------------------------
# Define the general netcdf generation function, common to all sensors
# ------------------------------------------------------------------------
#def generate_netcdf_common(nday,year,month,day,datestring):
def generate_netcdf_common(nday_local):

    epoch_offset = 719163
    #print(matplotlib.dates.get_epoch())	#Not available in older versions so can't use
    str_today = datetime.datetime.strftime(datetime.datetime.now(), '%Y%m%d')     #Today's date as a string
    num_today = int(mpl.dates.date2num(datetime.datetime(int(str_today[0:4]),int(str_today[4:6]),int(str_today[6:8]),0,0,0)))
    #print("In common, str_today, num_today = ", str_today, num_today)

    if num_today < 100000:
        nday_local = nday_local - epoch_offset

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------
    ndat             = num2date(nday_local)
    year             = ndat.year
    month            = ndat.month
    day              = ndat.day
    dat              = datetime.datetime(year,month,day,0,0,0)
    nowdat           = datetime.datetime.now()
    nowstring        = nowdat.strftime("%Y-%m-%d %H:%M:%S")
    datestring       = dat.strftime("%Y%m%d")
    start_of_day_str = dat.strftime("%Y-%m-%d %H:%M")

    #print(nday_local,ndat,year,month,day,dat,nowdat,nowstring,datestring)
    return year,month,day,datestring,start_of_day_str,nowstring



#--------------------------------------------
# Sorting out day/time information 
#--------------------------------------------
def generate_netcdf_datetimeinfo(year,month,day,n,timesecs):

    first_last_datetime = np.zeros((6,2),dtype = np.int32)      #Values to use in datetime.datetime, derived from time.gmtime
    tt = (year,month,day,0,0,0,0,0,0)
    tt_upper = (year,month,(day+1),0,0,0,0,0,0)
    day_start_epoch = time.mktime(tt)
    yr_arr=np.zeros(n)+year
    mn_arr=np.zeros(n)+month
    dy_arr=np.zeros(n)+day
    dyyr_arr=np.zeros(n)
    epoch_timesecs=np.zeros(n)
    hr_arr=np.zeros(n,dtype=np.int32)
    mi_arr=np.zeros(n,dtype=np.int32)
    sc_arr=np.zeros(n,dtype=np.float32)
    epoch_timesecs = day_start_epoch + timesecs
    #In cases where parts of previous or next day are found in file, it would be best to derive the above from the data timestamp
    for mm in range(n):
        tt = time.gmtime(epoch_timesecs[mm])
        hr_arr[mm] = np.int32(tt[3])
        mi_arr[mm] = np.int32(tt[4])
        sc_arr[mm] = timesecs[mm]-3600*hr_arr[mm]-60*mi_arr[mm]	#Keeps digits after decimal place, unlike tt[5] 
        dyyr_arr[mm] = np.float32(tt[7]+timesecs[mm]/86400.0)
        #If the last time is exactly midnight, we need different in
        if timesecs[mm] == 86400.0:     #Point from midnight of the following day
            sc_arr[mm] = sc_arr[mm] - 86400.0
            hr_arr[mm] = 24
            dyyr_arr[mm] = np.float32(tt[7])    #Will be number for next day
        #These values are used for first, last times in the global attributes
        #It's safer to derive them from the tt array, as using the arrays can cause problems with the datetime.datetime function
        #if values aren't compatible, e.g. last hr_arr value being 24
        if mm == 0:     #First point
            for nn in range(6):
                first_last_datetime[nn,0] = np.int32(tt[nn])
        if mm == n-1:   #Last point
            for nn in range(6):
                first_last_datetime[nn,1] = np.int32(tt[nn])

    #print(np.amin(hr_arr),np.amax(hr_arr),np.amin(mi_arr),np.amax(mi_arr),np.amin(sc_arr),np.amax(sc_arr))
    #print(sc_arr)
    #print(first_last_datetime)

    return yr_arr,mn_arr,dy_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs,first_last_datetime#,lengths_of_dimensions,substitutions



#--------------------------------------------
# Reading corrections file 
#--------------------------------------------
def load_netcdf_corrections(nday, chids, n_values, t_interval, timesecs, vals):

    epoch_offset = 719163
    n_chids = len(chids)
    qualflag = np.ones((n_values,n_chids), dtype = int)
    valid_min_max = np.zeros((2,n_chids))

    for n_inst in range(n_chids):
        corr_file = '/data/netCDF/corrections/' + chids[n_inst] + '.corr'
        print('n_inst, corr_file = ', n_inst, corr_file)

        sub_vals = vals[:, n_inst]	#Data for each instrument in file
        f = open(corr_file, 'r')

        #Check for HMP155 purgetimes first, as want any points flagged as purge time to be superseded by baddata from a correction file if it exists
        #Get hmp155 purge times for oatnew_ch and rhnew_ch
        if np.char.find(chids[n_inst], "oatnew") >= 0 or np.char.find(chids[n_inst], "rhnew") >= 0:
            print('HMP155 purge times')
            purge_limit = load_hmp155_purgetime(nday)
            start_purge = int(round(purge_limit[0]/t_interval))
            end_purge = int(round(purge_limit[1]/t_interval))
            #print('start_purge, end_purge = ', start_purge, end_purge)
            qualflag[start_purge:(end_purge+1),n_inst] = 3



        while True:	#Reading correction file

            line = f.readline().strip()
            if not line: break
            if len(line) == 30:	#Standard length for correction file
            #if len(line) == 31:	#Standard length for correction file
                startnum = int(date2num(datetime.datetime(int(line[0:4]),int(line[4:6]),int(line[6:8]),0,0,0)))
                #print 'nday, startnum = ',nday, startnum
                if startnum < 100000:
                    startnum = startnum + epoch_offset
                if startnum == nday:
                    #print(len(line), line)
                    print(line.strip())
                    t_start = 3600.0*float(line[9:11]) + 60.0*float(line[11:13]) + float(line[13:15])
                    t_end   = 3600.0*float(line[16:18]) + 60.0*float(line[18:20]) + float(line[20:22])
               	    #print(t_start, t_end, t_interval)
               	    #start_index = int(round(t_start/t_interval))	#Time range in file is 0 to 86390 
                    #end_index = int(round(t_end/t_interval))	#Don't think we need +1 here to round up.
               	    start_index = np.argmin(np.absolute(timesecs-t_start))	#Time range in file is 0 to 86390 
                    end_index = np.argmin(np.absolute(timesecs-t_end))	#Time range in file is 0 to 86390 
                    #end_index = np.where(t_end == timesecs)		#Don't think we need +1 here to round up.
                    #print('Correction start, end index, no. values in file = ',start_index, end_index, n_values)
                    if end_index >= (n_values-2):       #May be rounding errors with max. index so reduce if necessary. Also avoid last point in day being missed.
                        end_index = n_values - 1
                    if start_index == 1:        #Avoid having 1 single point remaining at start of day
                        start_index = 0
                    #print(start_index, end_index, n_values)
                    instring = line[23:]
                    if np.char.equal(instring,'HOLDCAL'):
                        if np.char.find(chids[n_inst],"rg") >= 0 or np.char.find(chids[n_inst],"disdrom") >= 0:
                            #print("qualflag for a raingauge")	#For a raingauge, only flag values as HOLDCAL if they're > 0
                            for flag_ind in range((end_index-start_index+1)):
                                #Set qualflag to 3 if the data value is greater than 0 and the flag isn't already set to 2 (BADDATA) from a previous line
                                if sub_vals[flag_ind+start_index] > 0 and qualflag[(flag_ind+start_index),n_inst] <= 1:
                                    qualflag[(flag_ind+start_index),n_inst] = 3
                                    #print(flag_ind+start_index, sub_vals[flag_ind+start_index])

                    if np.char.equal(instring,'BADDATA'):
                        qualflag[start_index:(end_index+1),n_inst] = 2	#Need end_index+1 to make change all points from start_index to (and including) end_index
                    print('Correction t_start(s), t_end(s), start_index, end_index, correction string = ', t_start, t_end, start_index, end_index, instring)

        f.close()

        vals_ok = sub_vals[np.where(qualflag[:,n_inst] == 1)]
        if len(vals_ok) > 0:
            valid_min_max[0,n_inst] = np.amin(vals_ok)
            valid_min_max[1,n_inst] = np.amax(vals_ok)

    return qualflag, valid_min_max


#--------------------------------------------
# Reading hmp155 purgetimes file 
#--------------------------------------------
def load_hmp155_purgetime(nday):

    purge_limit = np.zeros((2))
    epoch_offset = 719163

    hmp_file = '/data/netCDF/corrections/hmp155_purgetime.txt'

    f = open(hmp_file, 'r')

    while True:

        line = f.readline()
        if not line: break
        fdata = line.split(' ')
        date_start = int(date2num(datetime.datetime(int(line[0:4]),int(line[4:6]),int(line[6:8]),0,0,0)))
        date_end = int(date2num(datetime.datetime(int(line[16:20]),int(line[20:22]),int(line[22:24]),0,0,0)))
        if date_start < 100000:
            date_start = date_start + epoch_offset
            date_end = date_end + epoch_offset
       


        if nday > date_start and nday <= date_end:
            print(line)
            purge_limit[0] = 3600.0*float(fdata[4])
            purge_limit[1] = 3600.0*float(fdata[5])
            print('purge_limit = ', purge_limit)

    f.close()

    return purge_limit 


#-----------------------------------------------
# Reading broadband radiometer calibration files 
#-----------------------------------------------
def load_bbrad_calibrations(nday):

    scal = np.zeros((4))
    epoch_offset = 719163

    cm21_file = '/data/netCDF/corrections/CM21_calibration.txt'
    cmp21_file = '/data/netCDF/corrections/CMP21_calibration.txt'
    cg4_file = '/data/netCDF/corrections/CG4_calibration.txt'
    chp1_file = '/data/netCDF/corrections/CHP1_calibration.txt'
    file_in = [cm21_file, cmp21_file, cg4_file, chp1_file]

    for n in range(4):

        f = open(file_in[n], 'r')

        while True:

            line = f.readline()
            if not line: break
            fdata = line.split(' ')
            date_start = int(date2num(datetime.datetime(int(line[0:4]),int(line[4:6]),int(line[6:8]),0,0,0)))
            date_end = int(date2num(datetime.datetime(int(line[9:13]),int(line[13:15]),int(line[15:17]),0,0,0)))
            if date_start < 100000:
                date_start = date_start + epoch_offset
                date_end = date_end + epoch_offset
            #print(date_start, date_end)
            if nday >= date_start and nday < date_end:
                print(line)
                scal[n] = float(fdata[2])

        f.close()
    print('Broadband radiometer calibration factors = ',scal)
    return scal 





#--------------------------------------------
# Generate a plot of the data 
#--------------------------------------------
def generate_netcdf_graphs(out_file, graph_path_out, variables_to_plot, datestring):

    qc_threshold = 1	#Maximum qc_flagged value to plot in a 2d plot
    n_plots = len(variables_to_plot)
    
    tt = (int(datestring[0:4]),int(datestring[4:6]),int(datestring[6:8]),0,0,0,0,0,0)
    day_start = time.mktime(tt)

    ncfile = nc4.Dataset(out_file, 'r',format='NETCDF4_CLASSIC')
    ncfile.set_auto_maskandscale(False)
    splitstr = out_file.split('/')	#Split data filename at /, want last string of this array (filename)
    out_name = splitstr[len(splitstr)-1]
    point_pos = out_name.rfind('.')
    uscore_pos = out_name.find('_')
    #print(out_name, uscore_pos)
    #plot_file_root = out_name[0:point_pos]
    plot_file_root = out_name[0:uscore_pos+1]
    temp_time = ncfile.variables['time']
    var_time = temp_time[:]
    #var_time = (var_time - var_time[0])/3600.0
    var_time = (var_time - day_start)/3600.0
    n_points = var_time.size
    #print(n_points,var_time[0:9])
    x_axis_title = 'Time (hours since midnight)' 

    for n_meas in range(n_plots):

        plot_2d = 0	#Default is a 1d plot, but e.g. disdrometer counts data require a 2d plot
        #Derive output plot filename and create directory for it if it doesn't exist
        out_plot_path = graph_path_out + variables_to_plot[n_meas] + '/'
        #If this directory doesn't exist, create it
        if not os.path.isdir(out_plot_path):
            #It may be that higher directories will also be created by the os.makedirs command
            #If this is the case we need to change the permissions of these too
            #Work our way up the tree changing permissions, for 4 levels (will be above yyyy directory level)
            test_plot_path = graph_path_out + variables_to_plot[n_meas]	#No '/' at end of path
            last_fslash = [0,0,0,0] 
            dir_exist = np.zeros(4)
            for k in range(4):
                last_fslash[k] = test_plot_path.rfind('/')
                print(last_fslash[k], test_plot_path)
                #test_plot_path = test_plot_path[0: (last_fslash[k] + 1)]
                test_plot_path = test_plot_path[0: last_fslash[k]]
                if os.path.isdir(test_plot_path):
                    dir_exist[k] = 1
                
            print('last forward slash positions = ', last_fslash)
            print('directory already exists = ', dir_exist)

            print("Creating plot directory ", out_plot_path)
            os.makedirs(out_plot_path)
            oscommand = "chgrp -R netcdf " + out_plot_path
            os.system(oscommand) 
            #oscommand = "chmod -R g+w " + out_plot_path
            oscommand = "chmod -R 775 " + out_plot_path
            os.system(oscommand)

            for k in range(4):
                if dir_exist[k] == 0:
                    print('Changing permissions of directory ', out_plot_path[0 : last_fslash[k] + 1])
                    oscommand = "chgrp -R netcdf " + out_plot_path[0 : last_fslash[k] + 1]
                    os.system(oscommand) 
                    #oscommand = "chmod -R g+w " + out_plot_path
                    oscommand = "chmod -R 775 " + out_plot_path[0 : last_fslash[k] + 1]
                    os.system(oscommand)
                


        #out_plot_file = out_plot_path + plot_file_root + '_' + variables_to_plot[n_meas] + '.png'
        out_plot_file = out_plot_path + datestring + '_' + plot_file_root + variables_to_plot[n_meas] + '.png'  #Shorter name, starting with date, better for quicklooks and thumbnails
        print('out_plot_file = ',out_plot_file)

        #Mark any flagged points with a separate marker. Needs qc_flag name and values.
        qc_flag_read = 'qc_flag'	#Most are called qc_flag, but not all
        if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' and out_file.find('rain-gauge-4') >= 0:
            qc_flag_read = 'qc_flag_thickness_of_rainfall_amount'
        if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' and out_file.find('present-weather-1') >= 0:
            qc_flag_read = 'qc_flag_thickness_of_rainfall_amount'
        if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' and out_file.find('disdrometer') >= 0:
            qc_flag_read = 'qc_flag_thickness_of_rainfall_amount'
        if variables_to_plot[n_meas] == 'number_of_hydrometeors_per_size_channel' and out_file.find('disdrometer') >= 0:
            plot_2d = 1	#Requires a 2d plot 
            qc_flag_read = 'qc_flag_hydrometeors'
        if variables_to_plot[n_meas] == 'number_of_hydrometeors_per_size_channel' and out_file.find('present-weather-1') >= 0:
            plot_2d = 1	#Requires a 2d plot 
            qc_flag_read = 'qc_flag_hydrometeors'
        if variables_to_plot[n_meas] == 'rainfall_rate':
            qc_flag_read = 'qc_flag_rainfall_rate'
        if variables_to_plot[n_meas] == 'air_temperature':
            qc_flag_read = 'qc_flag_air_temperature'
        if variables_to_plot[n_meas] == 'air_temperature' and out_file.find('sonic-5') >= 0:
            qc_flag_read = 'qc_flag'
        if variables_to_plot[n_meas] == 'relative_humidity':
            qc_flag_read = 'qc_flag_relative_humidity'
        if variables_to_plot[n_meas] == 'wind_speed':
            qc_flag_read = 'qc_flag_wind_speed'
        if variables_to_plot[n_meas] == 'wind_from_direction':
            qc_flag_read = 'qc_flag_wind_from_direction'
        if variables_to_plot[n_meas] == 'downwelling_longwave_flux_in_air':
            qc_flag_read = 'qc_flag_downwelling_longwave_flux_in_air'
        if variables_to_plot[n_meas] == 'direct_downwelling_shortwave_flux_in_air':
            qc_flag_read = 'qc_flag_direct_downwelling_shortwave_flux_in_air'
        if variables_to_plot[n_meas] == 'downwelling_shortwave_flux_in_air' and out_file.find('radiometer-5') >= 0:
            qc_flag_read = 'qc_flag_downwelling_shortwave_flux_in_air'
        if variables_to_plot[n_meas] == 'upwelling_shortwave_flux_in_air':
            qc_flag_read = 'qc_flag_upwelling_shortwave_flux_in_air'
        if variables_to_plot[n_meas] == 'upwelling_longwave_flux_in_air':
            qc_flag_read = 'qc_flag_upwelling_longwave_flux_in_air'
        if variables_to_plot[n_meas] == 'body_temperature':
            qc_flag_read = 'qc_flag_body_temperature'
        if variables_to_plot[n_meas] == 'visibility_in_air':
            qc_flag_read = 'qc_flag_visibility_in_air'
        if variables_to_plot[n_meas] == 'wmo_synoptic_code':
            qc_flag_read = 'qc_flag_wmo_synoptic_code'

        #Read variable data from file
        temp_var_y = ncfile.variables[variables_to_plot[n_meas]]
        var_y = temp_var_y[:]
        y_axis_title = ncfile.variables[variables_to_plot[n_meas]].long_name + ' (' + ncfile.variables[variables_to_plot[n_meas]].units + ')'
        if plot_2d == 1:
            temp_sizes = ncfile.variables["diameter"]
            y_ds_scale = temp_sizes[:]
            n_sizes = y_ds_scale.size
            y_axis_title = ncfile.variables["diameter"].long_name + ' (' + ncfile.variables["diameter"].units + ')'
            cbar_label = ncfile.variables[variables_to_plot[n_meas]].long_name + ' (' + ncfile.variables[variables_to_plot[n_meas]].units + ')'
        temp_qc_flag = ncfile.variables[qc_flag_read]
        qc_flag_in = temp_qc_flag[:]
        #print(x_axis_title, y_axis_title)
        fig_title = ncfile.source + ' ' + datestring
        
        if plot_2d == 0:
            flagged_time = var_time[np.where(qc_flag_in > 1)]
            flagged_var_y = var_y[np.where(qc_flag_in > 1)]
            missing_value = ncfile.variables[variables_to_plot[n_meas]]._FillValue
            #y_axis_title = ncfile.variables[variables_to_plot[n_meas]].long_name + ' (' + ncfile.variables[variables_to_plot[n_meas]].units + ')'
            #fig_title = ncfile.source + ' ' + datestring
            #print(var_y[0:9])
            #print(x_axis_title, y_axis_title)
            #Calculate % missing values
            suspect_vals = 100.0 * flagged_var_y.size/n_points
            suspect_text = "Suspect: %.3f %%" % (suspect_vals) 
            print('suspect_text = ', suspect_text)

        else:	#2d plot
            flagged_time = var_time[np.where(qc_flag_in[:,0] > 1)]
            flagged_var_y = var_y[np.where(qc_flag_in > 1)]
            missing_value = ncfile.variables[variables_to_plot[n_meas]]._FillValue
            suspect_vals = 100.0 * flagged_var_y.size/(n_points * n_sizes)
            suspect_text = "Suspect: %.3f %%" % (suspect_vals) 
            print('suspect_text = ', suspect_text)
            var_masked = np.ma.masked_where(qc_flag_in > qc_threshold, var_y[:,:])    

        plt.ion()
        plt.ioff()
        fig1 = plt.subplots(nrows=1, constrained_layout = True)
        if plot_2d == 0:
        #plt.figure()
            plt.plot(var_time, var_y, linewidth=1)	#see https://matplotlib.org/2.0.0/api/pyplot_api.html#matplotlib.pyplot.legend and example code for more than 1 plot on figure
            #plt.plot(var_time, var_y, linewidth=1)	#see https://matplotlib.org/2.0.0/api/pyplot_api.html#matplotlib.pyplot.legend and example code for more than 1 plot on figure
            plt.plot(flagged_time, flagged_var_y, 'ro', fillstyle = 'none', markersize = 4, label = suspect_text)	#ro- joins points together, not wanted
            plt.xlabel(x_axis_title)
            plt.ylabel(y_axis_title)
            plt.xlim(0, 24)
            if variables_to_plot[n_meas] == 'wind_from_direction':
                plt.ylim(0, 360)
            if variables_to_plot[n_meas] == 'visibility_in_air':
                plt.ylim(0, 20000)
            if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' or variables_to_plot[n_meas] == 'rainfall_rate':
                print('Rainfall plot')
                min_y = np.amin(var_y)
                max_y = np.amax(var_y)
                #print('min, max_y = ', min_y, max_y)
                if min_y < 0:	#baddata present
                    min_y = -0.001
                if max_y <= 0.001:
                    print('Using minimum y-axis range')
                    min_y = -0.0001	#So that can see line at zero
                    max_y = 0.001    #Arbitrary low value
                plt.ylim(min_y, 1.2*max_y)
                #plt.text(suspect_text, 12, (0.9 * max_y))
            if variables_to_plot[n_meas] == 'northward_wind' or variables_to_plot[n_meas] == 'eastward_wind' or variables_to_plot[n_meas] == 'upward_air_velocity':
                if np.amin(var_y) < -100.0:     #A massive wind velocity, anything less than this is likely to be a missing value
                    min_ok = ncfile.variables[variables_to_plot[n_meas]].valid_min
                    max_ok = ncfile.variables[variables_to_plot[n_meas]].valid_max
                    min_y = 5.0*np.floor(min_ok/5.0)
                    max_y = 5.0*np.ceil(max_ok/5.0)
                    plt.ylim(min_y, max_y)
 
            plt.xticks( 4.0 * arange(7) )
            plt.suptitle(fig_title)
            plt.minorticks_on()
            plt.grid()
            plt.legend(loc=9)
            #plt.text(suspect_text, 12, (0.9 * max_y))

        else:	#plot_2d == 1
            #levels = MaxNLocator(nbins=15).tick_values(0.0, np.ceil(var_y.max()))    #Otherwise min is -Inf when plotting logs
            print('np.amin, np.amax counts = ', np.amin(var_y), np.amax(var_y))
            #print('var_y.shape = ',var_y.shape)
            #Set points where qc_flag != 1 to missing value
            #Can't use var_y = np.where(qc_flag_in == 1 , var_y, missing_value_array) as the output is an array, not values
            var_y[np.where(qc_flag_in > 1)] = missing_value	
            #for m in range(n_points):	#Using long-winded, clunky way instead!
                #for k in range(n_sizes):
                   #if qc_flag_in[m, k] != 1:
                       #var_y[m, k] = missing_value
            print('After qc_flag applied, np.amin, np.amax counts = ', np.amin(var_y), np.amax(var_y))

            #print('Ch disdro levels = ',levels)
            #cmap = plt.get_cmap('jet')
            #cmap.set_under("white")
            colormap = 'pyart_HomeyerRainbow'
            #Set up range for colorbar
            #For disdrometer, I'm using steps of 10 counts for vmax, but might want to make this more general for other insts.
            min_cbar = 0	#Can set this to e.g. 0.01, then zero values are white, but then they can't be distinguished from missing values
            max_cbar = 6
            #print('max(var_y), max_cbar = ',np.amax(var_y), max_cbar)
            #if max_cbar <= 1.0:
                #max_cbar = 5.0
            #print('max(var_y), max_cbar = ',np.amax(var_y), max_cbar)

            fig, ax = plt.subplots(figsize=[12,7])

            #plt.pcolormesh(var_time, y_ds_scale, np.transpose(var_y), cmap=cmap, vmin = min_cbar, vmax = max_cbar)       #Remove norm, so that vmin,vmax are used
            #cax = ax.pcolormesh(var_time, y_ds_scale, np.transpose(var_y), cmap=cmap, vmin = min_cbar, vmax = max_cbar)       #Remove norm, so that vmin,vmax are used
            #Use of masked array is better with colormap, as it shows as white where masked, not minimum value color
            cax = ax.pcolormesh(var_time, y_ds_scale, np.transpose(var_masked), cmap=colormap, vmin = 0.0, vmax = np.amax(var_y), shading = 'auto')
            if len(flagged_time > 0):
                plt.plot(flagged_time, (np.zeros(len(flagged_time))+y_ds_scale[0]), 'ro', fillstyle = 'none', markersize = 4)	#ro- joins points together, not wanted
            #cbar = fig.colorbar(cax, ticks = [0, 1, 2, 3, 4, 5, 6])
            cbar = fig.colorbar(cax, orientation='vertical')
            #cbar = colorbar(ticks = [0,1,2,3,4,5])
            #cbar.ax.set_yticklabels(['$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$', '$10^{4}$', '$10^{5}$', '$10^{6}$'])
            cbar.set_label(cbar_label)
            ax.set_xlabel(x_axis_title)
            ax.set_ylabel(y_axis_title)
            ax.set_xlim(0, 24)
            ax.set_xticks( 4.0 * arange(7) )
            ax.set_title(fig_title)
            ax.minorticks_on()
            #ax.xaxis.set_minor_locator(MultipleLocator(4))
            

        plt.savefig(out_plot_file,format='png')

        print('Changing permissions of file ',out_plot_file)
        oscommand = "chgrp netcdf " + out_plot_file
        os.system(oscommand) 
        #oscommand = "chmod g+w " + out_plot_file
        oscommand = "chmod 775 " + out_plot_file
        os.system(oscommand) 

    ncfile.close()
