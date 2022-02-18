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
fig,axs = plt.subplots(number_of_plots_in_X, number_of_plots_in_Y,figsize=(10.0, 6.5),sharex=True,sharey=True, gridspec_kw=dict(wspace=0.30, hspace=0.30))
plt.rcParams.update({'font.family':'serif'})
plt.rcParams.update({'font.serif':'Times New Roman'})
counter_x = -1
counter_y = 0

#Plot Controls

linestyle_str = ['solid', 'dotted','dashed','dashdot',(0, (3, 5, 1, 5, 1, 5))]  # Same as '-.'

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
            
            if counter_x== 0:
                axs[counter_x,counter_y].set_title("ALR="+str(j*100)+' %', fontsize=11)

            if counter_x== 2:

                axs[counter_x,counter_y].set_xlabel(r'Sd($t_{eff},\xi$) (in)')
                
            if counter_y== 0 and counter_x== 0:
                axs[counter_x,counter_y].set_ylabel(r'Serviciability P[$\varepsilon>\varepsilon_{S}$]')
            if counter_y== 0 and counter_x== 1:
                axs[counter_x,counter_y].set_ylabel(r'Damage Control P[$\varepsilon>\varepsilon_{DC}$]')
            if counter_y== 0 and counter_x== 2:
                axs[counter_x,counter_y].set_ylabel(r'Ultimate P[$\varepsilon>\varepsilon_{u}$]')

         
        counter_y+=1



fig.legend(['CL = 0%','CL = 5%', 'CL = 10%','CL = 15%','CL = 20%'],loc="lower center",ncol=5,
    bbox_to_anchor=(0.5,0.925),
    bbox_transform=fig.transFigure,edgecolor='black' )
plt.savefig("CDF_summary.pdf",dpi=600,bbox_inches='tight', pad_inches=0)