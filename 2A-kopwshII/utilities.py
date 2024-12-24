# -*- coding: utf-8 -*-
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal,fftpack
from scipy.ndimage.measurements import label
from scipy.signal import butter, lfilter
import numpy as np
from math import ceil

#***************************************FUNCTION SETUP*************************************************
def fft_filtering(sig,time_range):
    """This function IS NOT CURRENTLY IN USE 

    Args:
        sig (list): The input signal
        time_range (list): The time range
    """

    time_step = len(sig)/(time_range[1] - time_range[0])
    period = time_step

    # The FFT of the signal
    sig_fft = fftpack.fft(sig)

    # And the power (sig_fft is of complex dtype)
    power = np.abs(sig_fft)**2

    # The corresponding frequencies
    sample_freq = fftpack.fftfreq(len(sig), d=time_step)

    pos_mask = np.where(sample_freq > 0)
    freqs = sample_freq[pos_mask]
    peak_freq = freqs[power[pos_mask].argmax()]

    np.allclose(peak_freq, 1./period)

    high_freq_fft = sig_fft.copy()
    high_freq_fft[np.abs(sample_freq) > 1] = 0
    filtered_sig = fftpack.ifft(high_freq_fft)

    return(filtered_sig)


def butter_lowpass(cutoff, fs, order=5):
    """Calculates the parameters for the Butterworth lowpass filter

    Args:
        cutoff (float): The desired cutoff frequency for the filter application
        fs (float): Sampling frequency
        order (int, optional): [description]. Defaults to 5.

    Returns:
        float: the parameters a, b
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    """Filters the given data range. Uses the butter_lowpass function to calculate the specific parameters a,b

    Args:
        data (list): The data range
        cutoff (float): The desired cutoff frequency for the filter application
        fs (float): Sampling frequency
        order (int, optional): The order of the desired interpolation. Defaults to 5.

    Returns:
        list: The filtered data range
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def savitzky_golay_filter(data_range,no_of_points,order):
    """This function smoothens the given signal, using the Savitzky-Golay method

    Args:
        data_range (list): The data range
        no_of_points (int): The number of points for the interpolation procedure
        order (int): The order of the interpolation

    Returns:
        list: The smoothened signal
    """
    if not no_of_points % 2 or no_of_points <= 0 or not isinstance(no_of_points,int):
        print('Error: Number of  points must be positive, odd and integer')
    else:
        y = signal.savgol_filter(data_range,no_of_points,order)
        return y

def downsample_data(data_range,time_range,reduction_factor):
    """Performs a reduction of number of samples according to the specified reduction factor.
       For example, to get 100 times less data points, the reduction factor must be equal to 100

    Args:
        data_range (list): The data range
        time_range (list): The time range for the corresponding data
        reduction_factor (int): The downsizing factor

    Returns:
        list: The downsampled time and data ranges
    """
    No_data_points = len(data_range)
    reduced_length = ceil(No_data_points/ceil(reduction_factor))
    reduced_time_length = ceil(No_data_points/reduced_length)
    
    x_resampled = time_range[::reduced_time_length]
    y_resampled = signal.resample(data_range, reduced_length)

    return x_resampled,y_resampled

def reader(filename):
    """Reads the data from the specific .csv file

    Args:
        filename (string): The filename of the desired file

    Returns:
        dataframe: The data from the corresponding .csv file in a form of dataframe
    """   
    curr_path = os.getcwd() 
    os.chdir(curr_path) #change the working directory to current
    
    df = pd.read_csv(filename + '.csv', sep = ',')
    df_new = df.replace(',','.',regex = True)
    df_new = df_new.astype(float)

    return df_new

def writer(data,filename):
    """Writes the filtered data to a .csv file

    Args:
        data (dictionary): The data to be written in file
        filename (string): The name of the output file
    """
    df_final = pd.DataFrame(data)
    df_final.to_csv(filename)
