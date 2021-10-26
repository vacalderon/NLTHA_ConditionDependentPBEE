# -*- coding: utf= -8 -*-
"""
Created on Tue Oct 26 09:17:25 2021

@author: VACALDER
"""
# Program to select reduced number of mainshocks


import pandas as pd


data_file ='MainDF.csv'
column_file = 'column_database.csv'
groundmotions_file = 'mainshock_file_database_master.csv'
main_DF = pd.read_csv(data_file)
column_DF = pd.read_csv(column_file)
groundmotion_df = pd.read_csv(groundmotions_file)
iShapeFactor = [4, 6, 8]
iCL = [0, 5, 10, 15, 20]
steel_strain_limits = [0.0,0.015,0.05,0.15]
color_codes = ['#D94545','#73A8D4','#F7B76D','#519872','#57467B','#545066']

for j in range(0,3):
    corrosion_limit = (main_DF.CorrosionLvl_Long==5)
    lower_limit = steel_strain_limits[j]
    upper_limit = steel_strain_limits[j+1]
    limits_filter=(main_DF.Steel_Strain.between(lower_limit,upper_limit))
    selection_of_EQ_selection = main_DF[limits_filter & corrosion_limit]
    EQselection = selection_of_EQ_selection.sample(11)
    EQselection.to_csv('EQ_selection_LS_'+str(j)+'.csv')
