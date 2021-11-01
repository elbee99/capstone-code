import numpy as np
import matplotlib.pyplot as plt
import os
from math import pi
from matplotlib.ticker import AutoMinorLocator
import csv
from itertools import zip_longest
from shapely.geometry import Polygon
from tkinter import *
from tkinter import filedialog
 
    
def plot_all_cv_in_folder(folder,
                       saveas=None,
                       title=None,
                       double_scan=False,
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
                       find_peak_current = False):
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
        scan_label = "Scan "+str(iteration) #labels each plot as another scan
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
        
        V = V+potential_shift
        ax.plot(V,I,label=scan_label,linewidth=linewidth)
        
    
    ax.set_xlabel('Potential vs Ag/AgCl (mV)')
    ax.set_ylabel('Current (pA)')
    
    if normalize == True:
        if type(title) == str:
            title += ", current normalized"
        if type(saveas)== str:
            saveas += "_currentnorm"
            
    if capacitance_normalize == True:
        if type(title) == str:
            title += ", capacitance normlized"
        if type(saveas)== str:
            saveas += "_capacitancenorm"
            
            
    
    #finds the median and std dev of the peak currents reached
    #writes a file with the current values, median and std dev
    if find_peak_current == True:
        #peak_currents = sorted(peak_currents)
        median_peak_current = np.median(np.array(peak_currents))
        stddev_peak_current = np.std(np.array(peak_currents))
        header = ['Currents (pA)','Median','Standard Deviation']
        columns = [peak_currents,[median_peak_current],[stddev_peak_current]]
        with open(folder+r'/currents.csv','w+',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            
            for i in zip_longest(*columns):
                writer.writerow(i)
        
        
    #changes y label if we are reporting current density
    if current_density == True:
        ax.set_ylabel('Current density (mA/cm$^2$)')
    # ax.set_ylim(-0.2,0.2)
    # ax.set_xlim(-400,250)
    if legend ==True:
        ax.legend(loc='lower right')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_ticks_position('both') 

    ax.grid(True)
    if title!=None:
        ax.set_title(title)
    
    if saveas!=None:
        fig.savefig(saveas)
    
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

def find_capacitance(I,V,capacitance_pot_range = [-200,300],
                     scan_rate = 500):
    """

    Parameters
    ----------
    I : ARRAY. The current array for which we want to find the capacitance
    V : ARRAY. The voltage array for which we want to find the capacitance
    capacitance_pot_range: range of potential capacitance wants to be found
    in
    scan_rate: scan rate of CV, used to convert potential array to time array

    Returns
    -------
    capacitance, INT

    """
    
    closest_pot_min = V[0]
    for i in V:
        #print('Checking',i,'mv')
        if abs(i-capacitance_pot_range[0]) <= abs(closest_pot_min-capacitance_pot_range[0]):
            closest_pot_min = i
            
    closest_pot_max = V[0]
    for i in V:
        #print('Checking',i,'mv')
        if abs(i-capacitance_pot_range[1]) <= abs(closest_pot_max-capacitance_pot_range[1]):
            closest_pot_max = i
    #print('Closest potential found=',closest_pot)
    #finds the indices of both values of potential which are closest to the desired
    index_min = int(np.where(V == closest_pot_min)[0])
    index_max = int(np.where(V == closest_pot_max)[0])
    #print(index_min,V[index_min],index_max,V[index_max])
    
    #convert potential to time based on scan rate of 500 mV/s
    V = V/scan_rate
    potential_window = 0.001*(capacitance_pot_range[1] - capacitance_pot_range[0])
    
    #values of V and I for the forward and backward slices
    V_forward = V[index_max:index_min+1]
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
    integral = polygon.area
    capacitance = integral / (2*potential_window)
    
    
    return capacitance

def plot_dummy():
    x = np.linspace(0,1,100)
    y = x**2
    plt.figure(figsize = (1,1))
    plt.plot(x,y)
    plt.show()




def main():

    #gets default parameters from parameters.txt
    parameters = open('parameters.txt','r')
    parameter_lines = parameters.readlines()
    parameter_lines = [line.rstrip().split() for line in parameter_lines]
    double_scan = bool(int(parameter_lines[0][1]))
    linewidth = float(parameter_lines[1][1])
    legend = bool(int(parameter_lines[2][1]))
    normalize = bool(int(parameter_lines[3][1]))
    smoothen = bool(int(parameter_lines[4][1]))
    three_column = bool(int(parameter_lines[5][1]))
    potential_shift = float(parameter_lines[6][1])
    capacitance_normalize = bool(int(parameter_lines[6][1]))
    first_half_only = bool(int(parameter_lines[7][1]))
    find_peak_current = bool(int(parameter_lines[8][1]))


    root = Tk()
    #sets title of window
    root.title("CV Plotter")
    #allows use of file browser dialogue box to select desired directory
    root.dir_name = filedialog.askdirectory(initialdir=os.path.expanduser('~'))

    #defining all the widgets and their default values 
    
    #shows the selected folder path
    folder_label = Label(root,text="Selected folder: "+root.dir_name).grid(row=0,column=0,columnspan=3,sticky='w')

    #below is to get title 
    title_label = Label(root, text="Insert desired title here (LaTex syntax works for symbols): ")
    title_label.grid(row=1,column=0,sticky='w')
    title_input = Entry(root)
    title_input.grid(row=1,column=1,columnspan=2,sticky='w')
    
    potential_shift_label = Label(root, text="Insert potential shift (if desired):").grid(row=2,column=0,sticky='w')
    potential_shift_input = Entry(root)
    potential_shift_input.grid(row=2,column=1,columnspan=2,sticky='w')
    potential_shift_input.insert(0,str(potential_shift))

    legend_var = BooleanVar(root,legend) 
    legend_box = Checkbutton(root, text="Show legend on plot",variable = legend_var).grid(row=3,column=0,sticky='w')
    
    scan_labels_label = Label(root, text="Insert custom legend labels here if desired, seperated by commas (e.g. 'Label 1,Label 2, etc.'):")
    scan_labels_label.grid(row=3,column=1,sticky='w')
    scan_labels_input = Entry(root)
    scan_labels_input.grid(row=3,column=2,sticky='w')

    smoothen_var = BooleanVar(root,smoothen)
    smoothen_box = Checkbutton(root, text="Smoothen each CV to reduce noise", variable = smoothen_var).grid(row=4,column=0,sticky='w')

    linewidth_label = Label(root, text="Insert desired linewidth here").grid(row=4,column=1,sticky='w')
    linewidth_input = Entry(root)
    linewidth_input.grid(row=4,column=2,sticky='w')
    linewidth_input.insert(0,str(linewidth))


    
    normalize_var = BooleanVar(root,normalize)
    normalize_box = Checkbutton(root, text="Normalize each CV with respect to maximum current reached",variable=normalize_var).grid(row=5,column=0,sticky='w')

    capacitance_normalize_var = BooleanVar(root,capacitance_normalize)
    capacitance_normalize_box = Checkbutton(root, text="Normalize the CVs with respect to their capacitance",variable = capacitance_normalize_var).grid(row=5,column=2,sticky='w')


    double_scan_var = BooleanVar(root,double_scan)
    double_scan_box = Checkbutton(root, text = "Plot second of two scans taken",variable = double_scan_var).grid(row=6,column=0,sticky='w')
    
    first_half_only_var = BooleanVar(root,first_half_only)
    first_half_only_box = Checkbutton(root, text="Only show to forward scan",variable = first_half_only_var).grid(row=6,column=2,sticky='w')

    find_peak_current_var = BooleanVar(root,find_peak_current)
    find_peak_current_box = Checkbutton(root, text="Tabulate the largest current reached by each scan in a .csv file",variable = find_peak_current_var).grid(row=7,column=0,sticky='w')

    three_column_var = BooleanVar(root,three_column)
    three_column_box = Checkbutton(root, text="Data is in normal three column format", variable = three_column_var).grid(row=7,column=2,sticky='w')


    def plot_button_func():
        
        #getting scan labels if included
        scan_labels = scan_labels_input.get()
        if scan_labels == '':
            scan_labels = [] 
        else:
            scan_labels = scan_labels.split(',')

        #reads folder from the input given from the dialog box
        folder = root.dir_name
        plot_all_cv_in_folder(folder,
                   title=title_input.get(),
                   double_scan=double_scan_var.get(),
                   scan_labels = scan_labels,
                   current_density=False,
                   diameter=7e-5,
                   linewidth =float(linewidth_input.get()),
                   legend=legend_var.get(),
                   normalize=normalize_var.get(),
                   smoothen = smoothen_var.get(),
                   three_column = three_column_var.get(),
                   potential_shift = float(potential_shift_input.get()),
                   capacitance_normalize = capacitance_normalize_var.get(),
                   first_half_only = first_half_only_var.get(),
                   find_peak_current = find_peak_current_var.get())
    
    #this button is for testing only

    plot_button = Button(root, text = "Plot",command = plot_button_func) 
    plot_button.grid(row=8,column=1,sticky='w')
    root.mainloop()
    

if __name__=='__main__':
#    plot_all_cv_in_folder(r"/home/larry/Documents/Summer Project/05-08-21 1' a-C HER SECCM/Pipette 1/Sputtered/CVs/Export grid scan 1")
    main()
