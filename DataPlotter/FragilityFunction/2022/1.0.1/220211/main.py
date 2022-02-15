# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 13:35:07 2022

@author: VACALDER
"""
import pandas as pd
import numpy as np

data_file ='PosprocData_ALR0.20.csv'
iShapeFactor = [4, 6, 8]
iCL = [0, 5, 10, 15, 20]
iALR = [0.20]
main_DF = pd.read_csv(data_file)
# CODE TO GENERATE HYSTOGRAM
bins=np.arange(0, 250, 2)
length_of_bins=len(bins)

for i in iCL:
    counter =0
    SDs = []
    analyses = []
    LS_steel_collapses = []
    LS_concrete_collapses = []
    DC_steel_collapses = []
    DC_concrete_collapses = []
    Ultimate_collapses = []
    strains = []
    for j in range(length_of_bins-1):
        lower_limit = bins[counter]
        upper_limit = bins[counter+1]
        corrosion_level_filter=(main_DF.CorrosionLvl_Long==i)
        limits_filter=(main_DF.SD_Teff_xi.between(int(lower_limit),int(upper_limit)))
        corrosion_collapses=main_DF[corrosion_level_filter & limits_filter]
        mean_sd = (lower_limit+upper_limit)/2
        number_of_analyses = len(corrosion_collapses)
        Serviciability_Steel_collapses_count = corrosion_collapses['ServciabilitySteel'].sum()
        Serviciability_Concrete_collapses_count = corrosion_collapses['ServiciabilityConcrete'].sum()
        DamageControl_Steel_collapses_count = corrosion_collapses['DamageControlSteel'].sum()
        DamageControl_Concrete_collapses_count = corrosion_collapses['DamageControlConcrete'].sum()
        Ultimate_collapses_count = corrosion_collapses['UltimateSteel'].sum()
        SDs.append(mean_sd)
        analyses.append(number_of_analyses)
        LS_steel_collapses.append(Serviciability_Steel_collapses_count)
        LS_concrete_collapses.append(Serviciability_Concrete_collapses_count)
        DC_steel_collapses.append(DamageControl_Steel_collapses_count)
        DC_concrete_collapses.append(DamageControl_Concrete_collapses_count)
        Ultimate_collapses.append(Ultimate_collapses_count)
        counter = counter + 1
        
    ColapsesDict = {'SD_teff_xi': SDs,
                    'Number of Analyses': analyses,
                    'Collapses_LS_steel': LS_steel_collapses,
                    'Collapses_LS_concrete': LS_concrete_collapses,
                    'Collapses_DC_steel': DC_steel_collapses,
                    'Collapses_DC_concrete': DC_concrete_collapses,
                    'Collapses_ultiamte': Ultimate_collapses}
    DataFrame_CollapseResults = pd.DataFrame(ColapsesDict)
    DataFrame_CollapseResults.to_csv('collapseinfo_CL'+str(i)+'.csv')
    