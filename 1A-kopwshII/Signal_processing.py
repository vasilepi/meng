# -*- coding: utf-8 -*-
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal,fftpack
from scipy.ndimage import label
from scipy.signal import butter, lfilter
import numpy as np
from math import ceil
import utilities as utils

#********************************CODE BEGINS HERE********************************************
def main():

    #__________________________INPUT DATA SECTION_________________________________________________________
    #*Everything you will need to change are located in this section
    
    filename = "Flucht" #The given .csv file with the data. Change the name with the name of the file you want to read

    #******* Input data for the lowpass filter*********
    order_lowpass = 4
    fs = 60.0    # sample rate of the raw data, Hz
    cutoff = 1  # desired cutoff frequency of the filter, Hz

    #******Input data for the Savitzky-Golay method****
    order_svgolay = 2
    no_of_poins_svgolay = 41

    #**********Reduction factor for downsampling process******
    reduction_factor = 1 #change for lower or higher sampling rate, for the downsampling process
    #______________________________________________________________________________________________________
    
    #*****It is not necessary to touch anything from the code below*********

    df_new = utils.reader(filename) #read the file      
    measurements_points = df_new.keys()[1:] #Extract the names of each measurement point
    
    #Extract the time range for the measurements
    try:
        time_range = df_new.Zeit
    except:
        time_range = df_new.iloc[:, 0]

    temp_dict = {}
    iter = 1
    
    for i in measurements_points:
        
        # Filter the data, by removing the "noise" above the cutoff frequency.
        filtered_signal = utils.butter_lowpass_filter(df_new[i].to_numpy(), cutoff, fs, order_lowpass)

        #Smooth the filtered signal, using the Savitzky-Golay method       
        smoothened_signal = utils.savitzky_golay_filter(filtered_signal,no_of_poins_svgolay,order_svgolay)

        #Downsample the smoothened signal
        x,y = utils.downsample_data(smoothened_signal,list(time_range),reduction_factor)

        #Add Time column in final results table
        if not 'Time' in temp_dict:
            temp_dict['Time'] = x

        temp_dict[i] = y
        
        #Plot the signals
        plt.figure(iter)
        plt.title(i)
        plt.plot(time_range, list(df_new[i]),label = 'Raw data')
        # plt.plot(time_range,filtered_signal,label = 'Filtered signal')
        # plt.plot(time_range,smoothened_signal, label = 'Smoothened signal')
        plt.plot(x,y, label = 'Downsampled Data')
        plt.grid()
        plt.legend()

        iter += 1
        
    plt.show()
    utils.writer(temp_dict,filename + '_filtered2.csv') #write the data to file
    

if __name__ == "__main__":
    main()