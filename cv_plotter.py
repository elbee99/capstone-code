# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 17:32:38 2021

@author: Larry
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from math import pi
from matplotlib.ticker import AutoMinorLocator
from scipy import integrate
import csv
from itertools import zip_longest
from shapely.geometry import Polygon

#basic plotting code
def plot_single_cv(file, saveas=None, title=None,smoothen=True,
                   scan_number=None,total_scans=None,
                   split_2_scans=False,normalize=False,
                   three_column=False,
                   log = False,
                   potential_shift = 0,
                   give_capacitance = False):
    """
    

    Parameters
    ----------
    file : path to the desired .txt file containing
    CV data type: raw string
    saveas : location of desired saved plot. Raw string again
        DESCRIPTION. The default is None.
    title : string of desired title
        DESCRIPTION. The default is None.
    smooth : if True will smooth the noise using np.convolve, else will use
    raw data
        DESCRIPTION. The default is True.
    scan_number: if multiple scans are done in the same dataset this will plot
    the scan number given here e.g. if scan_number = 3 it will plot the third
    scan
    
    split_2_scans: if set to True will show the first and second scan as
    separate colours with separate labels
    
    three_columns: set to True if the export is in the 3 column format XEI
    does for SECCM scans
    
    log: if True it will plot V vs log(i) instead of i vs V
    
    give_capacitance: if True will find the capacitance and print it
    Returns
    -------
    None.

    """
    if three_column == True:
        rows_to_skip = 3
    else:
        rows_to_skip = 4
    data = np.genfromtxt(file,skip_header=rows_to_skip)
    if three_column != True:
        data = np.hsplit(data,2)
    else:
        data = np.hsplit(data,4)

    if three_column != True:
        V = data[0].flatten()
        I = data[1].flatten()
    else:
        V = data[1].flatten()
        I = data[2].flatten()
        I = I*10**3
        
    I_smooth = smooth(I,5)
    if smoothen == True:
        I = I_smooth
        
    if scan_number!=None and total_scans!=None:
        start_index = int((1024/total_scans)*(scan_number-1))
        V = V[start_index:start_index+int((1024/total_scans))]
        I = I[start_index:start_index+int((1024/total_scans))]
        
    if abs(np.max(V)) < 10:
        V = V*1000
        
    if give_capacitance == True:
        capacitance = find_capacitance(I, V)
        capacitance_string = 'Capacitance = {:.4f} pF'.format(capacitance)
        #print(capacitance_string)
    
    if normalize == True:
        I = normalize_func(True)
        
    
    
    V = V + potential_shift
    y_label = 'Current (pA)'
    if log == True:
        I = np.log10(-I)
        # I,V = V,I #to plot potential vs log(current) instead of vice versa
        y_label = 'log(current)'
    
    if split_2_scans!=True:
        fig, ax = plt.subplots(figsize=(8,4.5),dpi=100)
        ax.plot(V,I)
        ax.set_xlabel('Potential vs PdH (mV)')
        ax.set_ylabel(y_label)
        #ax.set_ylim(-3,1)
        #ax.set_xlim(-200,400)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_ticks_position('both') 
        ax.grid(True)
    
    if split_2_scans==True:
        V_list = [V[:512],V[512:]]
        I_list = [I[:512],I[512:]]
        fig, ax = plt.subplots(figsize=(8,4.5),dpi=100)
        for i in range(0,2):
            
            
            ax.plot(V_list[i],I_list[i],label='Cycle '+str(i+1))
            ax.set_xlabel('Potential vs PdH (mV)')
            ax.set_ylabel(y_label)
            #ax.set_ylim(-5.5,0)
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.xaxis.set_ticks_position('both') 
            ax.grid(True)
            ax.legend()
            
            
    if type(title)==str:
        ax.set_title(title)
    if give_capacitance == True:
        ax.text(-800,-0.5,capacitance_string,
                bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   ))
    if type(saveas)==str:
        fig.savefig(saveas)
    
def gen_data_from_folder(folder,three_column=True,double_scan=True):
    data_files = []
    file_names = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.txt'):
                if three_column != True:
                    file_names.append(file)
                    data_files.append(os.path.join(root,file))
                else:
                    if file[-8:-4] == 'info':
                        pass
                    else:
                        file_names.append(file)
                        data_files.append(os.path.join(root,file))

    for i in range(len(data_files)):
        iteration+=1
        file_name = file_names[i]
        scan_label = "Scan "+str(file_name) #labels each plot as another scan
        if three_column == True:
            rows_to_skip = 3
        else:
            rows_to_skip = 4
        data = np.genfromtxt(file,skip_header=rows_to_skip)
        
        #print(data[0])
        if three_column != True:
            data = np.hsplit(data,2)
        else:
            data = np.hsplit(data,4)

        if three_column != True:
            V = data[0].flatten()
            I = data[1].flatten()
        else:
            V = data[1].flatten()
            I = data[2].flatten()
        
        #if voltage is in volts rather than millivolts this converts to millivolts
        if abs(np.max(V)) < 10:
            V = V*1000
            
        if double_scan == True:
            V = V[512:]
            I = I[512:]
        if three_column == True:
            I = I*10**3
        return (V,I)
    
    
def plot_all_cv_in_folder(folder,
                       saveas=None,
                       title=None,
                       double_scan=True,
                       scan_labels = [],
                       current_density=False,
                       diameter=7e-5,
                       linewidth=1.0,
                       legend=True,
                       normalize=False,
                       smoothen = False,
                       three_column = True,
                       potential_shift = 0,
                       capacitance_normalize = False,
                       first_half_only = False,
                       find_peak_current = False,
                       custom_normalize = [],
                       colour = None,
                       no_show = False):
                        
    """

    Parameters
    ----------
    folder : The path to the folder containing the .txt files from the 
    CVs you need to plot.
    
    saveas: string of location to save image of plot
    
    title: string of desired title
    
    double_scan: some CVs are scanned twice and the data must be halved
    Usually the second scan is preferred so I have that here. Set to True if
    a double scan is used
    
    scan_labels: self defined labels if required, if left empty will simply
    number the scans Scan 1, Scan 2 etc
    
    current_density: converts current to current density given the pipette
    diameter IN CENTIMETRES, by default a small pipette of diameter 700 nm,
    
    three_column: set to True if data is exported in three column format
    
    capacitance_normalize: INT. Finds the capacitance at a given integer potential
    and divides the current by this in an attempt to normalize all plots by this
    
    first_scan_only: boolean, if set to True will only plot the first half of
    the cycle (forward scan only)
    
    find_peak_current: boolean, if set to True will record the peak current of
    each cycle, find the median and std dev and add these to a .csv file in 
    the folder defined by the folder argument

    custom_normalize: including a list of floats here will result in normalization
    of each of the plots by the value in the list. 
    Returns
    -------
    None.

    """
    data_files = []
    file_names = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.txt'):
                if three_column != True:
                    file_names.append(file)
                    data_files.append(os.path.join(root,file))
                else:
                    if file[-8:-4] == 'info':
                        pass
                    else:
                        file_names.append(file)
                        data_files.append(os.path.join(root,file))
    
    fig, ax = plt.subplots(figsize=(8,4.5),dpi=100)
    iteration = 0
    peak_currents = []
    for i in range(len(data_files)):
        iteration+=1
        file_name = file_names[i]
        scan_label = "Scan "+str(file_name) #labels each plot as another scan
        #changes scan labels to user defined ones 
        if scan_labels!=[]:
            scan_label = scan_labels[iteration-1]
        file = data_files[i]
        
        if three_column == True:
            rows_to_skip = 3
        else:
            rows_to_skip = 4
        data = np.genfromtxt(file,skip_header=rows_to_skip)
        
        #print(data[0])
        if three_column != True:
            data = np.hsplit(data,2)
        else:
            data = np.hsplit(data,4)

        if three_column != True:
            V = data[0].flatten()
            I = data[1].flatten()
        else:
            V = data[1].flatten()
            I = data[2].flatten()
        
        #if voltage is in volts rather than millivolts this converts to millivolts
        if abs(np.max(V)) < 10:
            V = V*1000
            
        if double_scan == True:
            V = V[512:]
            I = I[512:]
        if three_column == True:
            I = I*10**3
            
        #used if we want to report in current density rather than current
        if current_density==True:
            area = (pi/4)*(diameter**2) #droplet area in cm^2
            I = I*10e-9 #converts current from picoamps to milliamps
            I = I/area #converts current to current density
        
        if smoothen == True:
            I = smooth(I,10)

            
        if capacitance_normalize != False:
            capacitance = find_capacitance(I,V)
            #print(capacitance)
            scan_label += ", Capacitance = {:.3f} pF".format(capacitance)
            if capacitance != 0:
                
                #I = I / capacitance
                shift = np.mean(I[5:25])
                I = I - shift
            
        if normalize == True:
            shift = np.mean(I[5:50])
            I = I - shift
            I = normalize_func(I)
        
        if first_half_only == True:
            V = V[:len(V)//2]
            I = I[:len(I)//2]

        if find_peak_current == True:
            peak_current = I[len(I)//2]
            peak_currents.append(peak_current)
        
        if len(custom_normalize) == len(data_files):
            I = I/custom_normalize[i]

        V = V+potential_shift

        if type(colour)==str:
            ax.plot(V,I,label=scan_label,linewidth=linewidth,color=colour)
        else:
            ax.plot(V,I,label=scan_label,linewidth=linewidth)
        
    
    ax.set_xlabel('Potential vs PdH (mV)')
    ax.set_ylabel('Current (pA)')
        
    #finds the median and std dev of the peak currents reached
    #writes a file with the current values, median and std dev
    if find_peak_current == True:
        peak_currents = peak_currents
        median_peak_current = np.median(np.array(sorted(peak_currents)))
        stddev_peak_current = np.std(np.array(sorted(peak_currents)))
        header = ['Currents (pA)','Median','Standard Deviation']
        columns = [peak_currents,[median_peak_current],[stddev_peak_current]]
        with open(folder+r'\currents.csv','w+',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            
            for i in zip_longest(*columns):
                writer.writerow(i)
        
        
    #changes y label if we are reporting current density
    if current_density == True:
        ax.set_ylabel('Current density (mA/cm$^2$)')
    #ax.set_ylim(-0.2,0.2)
    #ax.set_xlim(-400,250)
    if legend ==True:
        ax.legend(loc='lower right')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.set_xlim(-1000,0)
    #ax.set_ylim(-50,0)
    ax.grid(True)
    if title!=None:
        ax.set_title(title)
    
    if saveas!=None:
        fig.savefig(saveas)
    if no_show == False:
        plt.show()
            


def smooth(y,box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y,box,mode='same')
    return y_smooth

def normalize_func(array):
    """

    Parameters
    ----------
    array : a NumPy array one wishes to normalize (make such that the max
    value = 1)

    Returns
    -------
    Returns the given array normalized such that the max value = 1

    """
    return array/max(abs(array))

def find_capacitance(I,V,capacitance_pot_range = [-0.4,-0.03]):
    """

    Parameters
    ----------
    I : ARRAY. The current array for which we want to find the capacitance
    V : ARRAY. The voltage array for which we want to find the capacitance
    capacitance_pot : list of 2 ints. The potential range over which we will
    determine the capacitance

    Returns
    -------
    capacitance, INT

    """
    # print("len(v)/2)=",len(V)/2)
    first_half_V = V[0:int(len(V)/2)]
    closest_pot_min = V[0]
    for i in first_half_V:
        #print('Checking',i,'mv')
        if abs(i-capacitance_pot_range[0]) <= abs(closest_pot_min-capacitance_pot_range[0]):
            closest_pot_min = i
            
    closest_pot_max = V[0]
    for i in first_half_V: 
        #print('Checking',i,'mv')
        if abs(i-capacitance_pot_range[1]) <= abs(closest_pot_max-capacitance_pot_range[1]):
            closest_pot_max = i
    # print('Closest potential found=',closest_pot_min,closest_pot_max)
    #finds the indices of both values of potential which are closest to the desired
    index_min = int(np.where(V == closest_pot_min)[0])
    index_max = int(np.where(V == closest_pot_max)[0])
    # print(index_min,V[index_min],index_max,V[index_max])
    
    #convert potential to time based on scan rate of 500 mV/s
    V = V/0.5
    
    #values of V and I for the forward and backward slices
    V_forward = V[index_max:index_min+1]
    # print("index_max,index_min=",index_max,index_min)
    I_forward = I[index_max:index_min+1]
    V_backward = V[-index_min:-index_max+1]
    I_backward = I[-index_min:-index_max+1]
    
    #define polygon points on forward curve
    polygon_points = []
    for i in range(len(V_forward)):
        polygon_points.append((V_forward[i],I_forward[i]))
        
    #define polygon points on reverse curve
    for i in range(len(V_backward)):
        polygon_points.append((V_backward[i],I_backward[i]))
    
    #define first point again to close the loop
    polygon_points.append((V_forward[0],I_forward[0]))
        
    #create the polygon and find its area
    polygon = Polygon(polygon_points)
    area = polygon.area
    # print(area)
    return area
    



if __name__=='__main__':
    pass

    folder_path = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Cu NWs (0.1M acid)/Export line scan 3 edit"
    plot_all_cv_in_folder(folder_path,
            title = "HER SECCM on Cu NWs and anC substrate in 0.1 M H$_2$SO$_4$"
            )
