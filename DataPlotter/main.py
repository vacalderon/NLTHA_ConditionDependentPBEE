# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 13:35:27 2021

@author: VACALDER
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from LibUnitsMUS import *
import seaborn as sns

data_file ='PosprocData.csv'
column_file = 'column_database.csv'
groundmotions_file = 'mainshock_file_database_master.csv'
main_DF = pd.read_csv(data_file)
column_DF = pd.read_csv(column_file)
groundmotion_df = pd.read_csv(groundmotions_file)
iShapeFactor = [4, 6, 8]
iCL = [0, 5, 10, 15, 20]
iALR = [0.05, 0.1, 0.15, 0.2]
compressive_strength_concrete = 5 * ksi
SD_effs = []
strain_LS = []
strain_DC = []
color_codes = ['#D94545','#73A8D4','#F7B76D','#519872','#57467B','#545066']

strain_ratios = []
Sd_eff_eq_xi = []
CorrosionLevels = []
final_vol_ratios = []
final_trans_ratios = []
alrs = []
shape_factors =[]

for i, row in main_DF.iterrows():
    D = float(row['DCol_in'])
    rhol = float(row['rhol'])
    CL = float(row['CorrosionLvl_Long'])
    delta_u = float(row['MaxDisplacement at MaxForce'])
    diameter_filter = (column_DF.column_diameter==D)
    volumetric_ratio_filter = (column_DF.rho_l==rhol)
    column_data = column_DF[diameter_filter & volumetric_ratio_filter]
    db = float(column_data['long_bar_diameter'].to_string(index=False))
    dbc = db*np.sqrt(1-CL*0.01)
    Fy = float(row['Fy_ksi'])
    L = float(row['height_of_column'])
    Lsp = 0.15*Fy*dbc
    epsilon_y = Fy/29000
    phi_y = 2.25*epsilon_y/D
    delta_y = phi_y*(L+Lsp)**2/3
    mu = abs(delta_u)/delta_y
    SD_el = float(row['SD_at_Teff'])
    e_steel = float(row['Steel_Strain'])
    e_steel_LS = 0.015
    e_steel_DC = float(row['LongitudinlSteelBuckling'])
    
    if e_steel<e_steel_LS:
        LS = 0
    elif e_steel>e_steel_LS:
        LS = 1
    strain_LS.append(LS)
    
    if e_steel<e_steel_DC:
        DC = 0
    elif e_steel>e_steel_DC:
        DC = 1    
    strain_DC.append(DC)
    
    if mu >1:
        xi_eq = 0.05 + 0.565*(mu-1)/(mu*np.pi)
        DF = np.sqrt((0.07)/(0.05+xi_eq))
        SD_Teff_xi_eq = DF * SD_el
        
    elif mu<=1:
        SD_Teff_xi_eq = SD_el
    
    SD_effs.append(SD_Teff_xi_eq)

main_DF.insert(27,"SD_Teff_xi",SD_effs,True)
main_DF.insert(28,"LS_collapse",strain_LS,True)
main_DF.insert(29,"DC_collapse",strain_DC,True)
#main_DF.to_csv('MainDF.csv')

rhol = [0.01,0.02,0.03,0.04]    
figure_number = 0
for i in iCL:
    figure_number = figure_number+1
    color = -1
    for j in rhol:
        color = color +1
        corrosion_level_filter=(main_DF.CorrosionLvl_Long==i)
        long_steel_rho_filter=(main_DF.rhol==j)
        corrosion=main_DF[corrosion_level_filter & long_steel_rho_filter]
        plt.figure(figure_number)
        plt.figure(figure_number).set_size_inches(5,5)
        plt.scatter(corrosion['SD_Teff_xi'],corrosion['Steel_Strain'],c=color_codes[color],label=r'$\rho_{l}=$'+str(rhol[color]))
        plt.title('Results for CL '+str(i), fontsize=12)
        plt.xlabel('Spectral Diplacement at Effective Period\n and Equivalent Damping'+r' $Sd(T_{eff},\xi_{eq}$)(in)', fontsize=10)
        plt.ylabel('Maximum steel strain'+r'$\varepsilon$ (in/in)', fontsize=10)
        plt.tick_params(direction='out', axis='both', labelsize=8)
        plt.grid()
        plt.legend(fontsize=6)
        plt.show()
        plt.figure(10+figure_number)
        count, bins, ignored = plt.hist(corrosion['SD_Teff_xi'], 30)
        plt.show()
        
        
for k, row in groundmotion_df.iterrows():
    groundmotion_name=row['horizontal_1_filename']
    
    for h, crow in column_DF.iterrows():
        D = crow['column_diameter']
        vol_ratio =crow['rho_l']
        trans_ratio = crow['rho_v']
        
        for shapefactor in iShapeFactor:
            for axial_load_ratio in iALR:
                
                Diameter_Column = float(D)*inch
                H = shapefactor*Diameter_Column
                Ag = 0.25 * np.pi * Diameter_Column ** 2
                AxialLoad = compressive_strength_concrete * Ag * axial_load_ratio
                Height_filter = (main_DF.height_of_column == H)
                Load_filter = (main_DF.PCol_kip == AxialLoad)
                GM_filter = (main_DF.earthquake == groundmotion_name)
                D_filter = (main_DF.DCol_in==D)
                volr_filter = (main_DF.rhol == vol_ratio)
                trans_filter = (main_DF.rhov == trans_ratio)
                subset_DF = main_DF[GM_filter & D_filter & volr_filter & trans_filter & Height_filter & Load_filter]
                
                for i in iCL:
                    if i == 0:
                        corrosion_filter = (subset_DF.CorrosionLvl_Long==i)
                        benchmark_data = subset_DF[corrosion_filter]
                        if benchmark_data.empty == True:
                            benchmark_strain =0
                            continue
                        else:
                            benchmark_strain = float(benchmark_data['Steel_Strain'].to_string(index=False))
                            corrosion_strain = float(benchmark_data['Steel_Strain'].to_string(index=False))
                            strain_ratio = (corrosion_strain/benchmark_strain)-1
                            strain_ratios.append(strain_ratio)
                            Sd_eff_eq_xi.append(float(benchmark_data['SD_Teff_xi'].to_string(index=False)))
                            alrs.append(axial_load_ratio)
                            CorrosionLevels.append(i)
                            final_vol_ratios.append(vol_ratio)
                            final_trans_ratios.append(trans_ratio)
                            shape_factors.append(shapefactor)
                                                       
                        
                    elif i !=0:
                        corrosion_filter = (subset_DF.CorrosionLvl_Long==i)
                        corrosion_data = subset_DF[corrosion_filter]
                        if corrosion_data.empty == True or benchmark_strain == 0:
                            continue
                        else:
                            corrosion_filter = (subset_DF.CorrosionLvl_Long==i)
                            corrosion_data = subset_DF[corrosion_filter]
                            corrosion_strain = float(corrosion_data['Steel_Strain'].to_string(index=False))
                            strain_ratio = (corrosion_strain/benchmark_strain)-1
                            strain_ratios.append(strain_ratio)
                            Sd_eff_eq_xi.append(float(benchmark_data['SD_Teff_xi'].to_string(index=False)))
                            alrs.append(axial_load_ratio)
                            CorrosionLevels.append(i)
                            final_vol_ratios.append(vol_ratio)
                            final_trans_ratios.append(trans_ratio)
                            shape_factors.append(shapefactor)
            
dataDict = {'steel_strain_ratio': strain_ratios,
            'spectral_displacement_teff_xi': Sd_eff_eq_xi,
            'CorrosionLvl_Long': CorrosionLevels,
            'volumetric_ratio': final_vol_ratios,
            'transverse_ratio': final_trans_ratios,
            'axial_load_ratio': alrs,
            'shape_factor': shape_factors,}
    
DataFrame_RatiosResults = pd.DataFrame(dataDict)        

color = -1


for i in iCL:
    figure_number = figure_number+1
    color = color +1
    corrosion_level_filter=(DataFrame_RatiosResults.CorrosionLvl_Long==i)
    corrosion=DataFrame_RatiosResults[corrosion_level_filter]
    plt.figure(figure_number)
    plt.figure(figure_number).set_size_inches(5,5)
    plt.scatter(corrosion['spectral_displacement_teff_xi'],corrosion['steel_strain_ratio'],c=color_codes[color],label=r'CL='+str(i))
    plt.title('Strain increase due to corrosion', fontsize=12)
    plt.xlabel('Spectral Diplacement at Effective Period\n and Equivalent Damping'+r' $Sd(T_{eff},\xi_{eq}$)(in)', fontsize=10)
    plt.ylabel('Steel strain ratio'+r'$\varepsilon_{CL}$/\varepsilon_{0} (in/in)', fontsize=10)
    plt.tick_params(direction='out', axis='both', labelsize=8)
    plt.grid()
    plt.legend(fontsize=6)
    plt.show()

plt.figure(figure_number+1)
sns.boxplot(x=DataFrame_RatiosResults["CorrosionLvl_Long"], y=DataFrame_RatiosResults["steel_strain_ratio"])
plt.ylim(-0.5,2)
plt.show()

# CODE TO GENERATE HYSTOGRAM
bins=np.arange(0, 250, 2)
length_of_bins=len(bins)

for i in iCL:
    counter =0
    SDs = []
    analyses = []
    LS_collapses = []
    DC_collapses = []
    strains = []
    for j in range(length_of_bins-1):
        lower_limit = bins[counter]
        upper_limit = bins[counter+1]
        corrosion_level_filter=(main_DF.CorrosionLvl_Long==i)
        limits_filter=(main_DF.SD_Teff_xi.between(int(lower_limit),int(upper_limit)))
        corrosion_collapses=main_DF[corrosion_level_filter & limits_filter]
        mean_sd = (lower_limit+upper_limit)/2
        number_of_analyses = len(corrosion_collapses)
        LS_collapses_count = corrosion_collapses['LS_collapse'].sum()
        DC_collapses_count = corrosion_collapses['DC_collapse'].sum()
        SDs.append(mean_sd)
        analyses.append(number_of_analyses)
        LS_collapses.append(LS_collapses_count)
        DC_collapses.append(DC_collapses_count)
        counter = counter + 1
        
    ColapsesDict = {'SD_teff_xi': SDs,
                    'Number of Analyses': analyses,
                    'Collapses_LS': LS_collapses,
                    'Collapses_DC': DC_collapses}
    DataFrame_CollapseResults = pd.DataFrame(ColapsesDict)
    DataFrame_CollapseResults.to_csv('collapseinfo_CL'+str(i)+'.csv')
    
    
        
        
        
        
        
    
    