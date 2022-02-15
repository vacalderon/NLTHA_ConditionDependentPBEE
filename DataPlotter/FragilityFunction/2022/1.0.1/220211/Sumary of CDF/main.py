# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 10:10:49 2022

@author: VACALDER
"""

#PROGRAM TO PLOT CDF PLOTS IN CONDITION DEPENDENT PERFORMANCE BASED DESIGN
#V1.0.1
#Ph.D. Victor A. Calderon
#2022(c)


import pandas as pd
import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt

corrosion_levels = [0,5,10,15,20]
axial_load_ratios = [0.05,0.10,0.15,0.20]
limit_states = ['serviciability','damage_control','ultimate']
data_file ='cdf_summary.csv'
colors = ['#d94545','#519872','#73a8d4','#f7b76d','#545066']
main_DF = pd.read_csv(data_file)

number_of_plots_in_X = len(limit_states)
number_of_plots_in_Y = len(axial_load_ratios)
fig,axs = plt.subplots(number_of_plots_in_X, number_of_plots_in_Y)

counter_x = -1
counter_y = 0

#Plot Controls

linestyle_str = ['solid', 'dotted','dashed','dashdot','dashed']  # Same as '-.'

linestyle_tuple = [
     ('solid',                 (0, ())), 
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]


for i in limit_states:
    counter_x+=1 
    counter_y = 0
    for j in axial_load_ratios:
        z=-1
        for k in corrosion_levels:
            
            z+=1                
            lyne_type=linestyle_str[z]
            plot_color = colors[z]
            corrosion_level_filter=(main_DF.CL==k)
            ALR_filter=(main_DF.ALR==j)
            limit_state_filter=(main_DF.limit_state==i)
            CDF_info=main_DF[corrosion_level_filter & ALR_filter & limit_state_filter]
            mean_value = float(CDF_info['mean'])
            standard_deviation = float(CDF_info['std_dev'])
            
            data = np.random.randint(2,100,size=1000)
            
            #sort data
            x = np.sort(data)
            
            #calculate CDF values
            y = sps.norm.cdf(np.log(x),loc=np.log(mean_value),scale=standard_deviation)
            
            
            #plot CDF
            axs[counter_x,counter_y].plot(x, y,color=plot_color,linestyle=lyne_type)
            print('x=',counter_x)
            print('y=',counter_y)
        counter_y+=1