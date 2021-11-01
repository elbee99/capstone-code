import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from matplotlib.ticker import AutoMinorLocator
from shapely.geometry import Polygon
from matplotlib.lines import Line2D
print('weird')


def gen_text_from_folder(folder):
    """
    Takes a folder as argument and returns a tuple of (file_names
    , text files) themselves as lists
    """
    data_files = []
    file_names = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.txt'):
                if file[-8:-4] == 'info':
                    pass
                else:
                    file_names.append(file)
                    data_files.append(os.path.join(root,file))
    return file_names,data_files

def gen_array_from_text(data_file, file_name):
    """
    Takes a data_file and file_name as arguments as strings
    Returns a (V,I) tuple of potential and current arrays
    """
    data = np.genfromtxt(data_file,skip_header=3)
    data = np.hsplit(data,4)
    V = data[1].flatten()
    I = data[2].flatten()
    return V,I


def convert_cv_to_time(V,I,scan_rate):
    """

    Parameters
    ----------
    V: potential array
    I: current array
    scan_rate: float of scan rate used
    Returns
    -------
    Returns I,T arrays of current and time respectively

    """
    potential_window=V[20]-V[5]
    total_time = abs(potential_window/scan_rate)
    time_between_points = total_time/len(V[5:20])
    T = []
    time = 0
    for i in range(len(I)):
        T.append(time)
        time+=time_between_points

    T = np.array(T)
    return I,T

def array_modifier(V,I,
                double_scan = True,
                smoothen = False,
                normalize = False,
                capacitance_normalize=False,
                first_half_only = False):

        """
        Modifies the arrays in desired ways, returning the (V,I) tuple
        """

        if double_scan == True:
            V = V[512:]
            I = I[512:]

        if smoothen == True:
            I = smooth(I,10)

        if capacitance_normalize != False:
            capacitance = find_capacitance(I,V)
            # print(capacitance)
            #scan_label += ", Capacitance = {:.3f} pF".format(capacitance)
            if capacitance != 0:
                
                I = I / capacitance
                shift = np.mean(I[5:25])
                I = I - shift
            
        if normalize == True:
            shift = np.mean(I[5:50])
            I = I - shift
            I = normalize_func(I)
        
        if first_half_only == True:
            V = V[:len(V)//2]
            I = I[:len(I)//2]

        return V,I

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

def plotter(folders=[],colors=[],
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
           capacitance_normalize = False,
           first_half_only = False,
        custom_legend = False,
        xlabel = None,
        ylabel = None):


    fig, ax = plt.subplots(figsize=(8,4.5),dpi=100)
    print('plot')
    for i in range(len(folders)):
        folder = folders[i]
        file_names,data_files = gen_text_from_folder(folder)
        V_list = []
        I_list = []
        # print(data_files)
        for j in range(len(data_files)):
            data_file = data_files[j]
            file_name = file_names[j]
            V, I = gen_array_from_text(data_file,file_name)
            V,I = array_modifier(V,I,
                    double_scan = double_scan,
                    smoothen=smoothen,
                    normalize = normalize,
                    capacitance_normalize=capacitance_normalize,
                    first_half_only=first_half_only)
            ax.plot(V,I,color=colors[i],linewidth=linewidth)

                    
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title!=None:
        if normalize == True:
            title = title + ", normalized by current"
        if capacitance_normalize == True:
            title = title + ", capacitance normalized"
        ax.set_title(title)
    
    if custom_legend == True:
        custom_lines = [Line2D([0], [0], color=colors[0], lw=4),
                Line2D([0], [0], color=colors[1], lw=4)]
        ax.legend(custom_lines, ['Substrate only','Cu NWs on substrate'])
    if saveas!=None:
        fig.savefig(saveas)
        print('saved')
    plt.show()

if __name__ == "__main__":

    #SECCM on Cu NWs in 0.1M H2SO4
    substrate_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Cu NWs (0.1M acid)/linescan_3_substrate" 
    dropcast_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Cu NWs (0.1M acid)/linescan_3_dropcast"
    saveas = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Cu NWs (0.1M acid)/linescan3_comparison"
    plotter(folders=[substrate_folder,dropcast_folder],colors=['r','b'],
            custom_legend=True,linewidth=0.7,
            saveas=saveas,
            title = "HER SECCM on anC vs anC/Ni NPs w/ 0.1 M H$_2$SO$_4$",
            xlabel="Potential vs PdH (V)",ylabel="Current (nA)",
            # first_half_only=True,
            # capacitance_normalize=True,
            )


    ##SECCM on Ni NPs in 0.1M H2SO4
    #substrate_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Ni NPs (0.1M acid)/Pipette 2/linescan_2_substrate" 
    #dropcast_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Ni NPs (0.1M acid)/Pipette 2/linescan_2_dropcast"
    #saveas = r"/home/larry/Documents/Capstone Project/Data/21-10-07 SECCM on Ni NPs (0.1M acid)/Pipette 2/linescan2_comparison"
    #plotter(folders=[substrate_folder,dropcast_folder],colors=['r','b'],
    #        custom_legend=True,linewidth=0.7,
    #        saveas=saveas,
    #        title = "HER SECCM on anC vs anC/Ni NPs w/ 0.1 M H$_2$SO$_4$",
    #        xlabel="Potential vs PdH (V)",ylabel="Current (nA)",
    #        # first_half_only=True,
    #        # capacitance_normalize=True,
    #        )

    ##SECCM on Ni NPs at pH=12
    #substrate_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-06 SECCM on Ni NPs (base)/Pipette 1/Substrate and dropcast/scan_1_substrate" 
    #dropcast_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-06 SECCM on Ni NPs (base)/Pipette 1/Substrate and dropcast/scan_1_dropcast"
    #saveas = r"/home/larry/Documents/Capstone Project/Data/21-10-06 SECCM on Ni NPs (base)/Pipette 1/Substrate and dropcast/scan_1_comparison_capacitancenorm"
    #plotter(folders=[substrate_folder,dropcast_folder],colors=['r','b'],
    #        custom_legend=True,linewidth=0.7,
    #        saveas=saveas,
    #        title = "HER SECCM on anC vs anC/Ni NPs at pH=12",
    #        xlabel="Potential vs Pt (V)",ylabel="Current (nA)",
    #        # first_half_only=True,
    #        capacitance_normalize=True,
    #        )



    ##SECCM on Ni NPs at pH=7
    #substrate_folder = r"/home/larry/Documents/Capstone Project/Data/21-09-28 SECCM on Ni NPs/Medium pipette 2/scan_1_substrate" 
    #dropcast_folder = r"/home/larry/Documents/Capstone Project/Data/21-09-28 SECCM on Ni NPs/Medium pipette 2/scan_1_dropcast"
    #saveas = r"/home/larry/Documents/Capstone Project/Data/21-09-28 SECCM on Ni NPs/Medium pipette 2/scan_1_comparison_capacitancenorm"
    #plotter(folders=[substrate_folder,dropcast_folder],colors=['r','b'],
    #        custom_legend=True,linewidth=0.7,
    #        saveas=saveas,
    #        title = "HER SECCM on anC vs anC/Ni NPs at pH=7",
    #        xlabel="Potential vs Ag/AgCl (V)",ylabel="Current, capacitance normalized",
    #        # first_half_only=True,
    #        capacitance_normalize=True,
    #        )

     # # SECCM on Cu NWs at pH=7
     # substrate_folder = r"/home/larry/Documents/Capstone Project/Data/21-09-29 SECCM on Cu NWs/Pipette 3/Substrate and dropcast/scan_1_substrate"
     # dropcast_folder = r"/home/larry/Documents/Capstone Project/Data/21-09-29 SECCM on Cu NWs/Pipette 3/Substrate and dropcast/scan_1_dropcast"
     # saveas = r"/home/larry/Documents/Capstone Project/Data/21-09-29 SECCM on Cu NWs/Pipette 3/Substrate and dropcast/scan_1_comparison_capacitancenorm"
     # plotter(folders=[substrate_folder,dropcast_folder],colors=['r','b'],
     #         custom_legend=True,linewidth=0.7,
     #         saveas=saveas,
     #         title = "HER SECCM on anC vs anC/Cu NWs at pH=7",
     #         xlabel="Potential vs Ag/AgCl (V)",ylabel="Current, capacitance normalized",
     #         # first_half_only=True,
     #         capacitance_normalize=True,
     #         )
    
    # # SECCM on Ni NPs at pH=2
    # substrate_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-04 SECCM on Ni NPs (acid)/Pipette 4/Substrate and drop-cast/scan_1_substrate"
    # dropcast_folder = r"/home/larry/Documents/Capstone Project/Data/21-10-04 SECCM on Ni NPs (acid)/Pipette 4/Substrate and drop-cast/scan_1_dropcast"
    # saveas = r"/home/larry/Documents/Capstone Project/Data/21-10-04 SECCM on Ni NPs (acid)/Pipette 4/Substrate and drop-cast/scan_1_comparison_capacitancenorm_linear"
    # # print(dropcast_folder)
    # plotter(folders=[substrate_folder,dropcast_folder],colors=['r','b'],
    # custom_legend=True,
    # saveas = saveas,
    # title = "HER SECCM on anC vs anC/Ni NPs at pH=2",
    # first_half_only=True,
    # xlabel="Potential vs PdH (V)",ylabel="Current normalized by capacitance",
    # capacitance_normalize=True,
    # )
