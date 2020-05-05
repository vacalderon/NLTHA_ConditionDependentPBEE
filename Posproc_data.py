# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 13:23:30 2019

@author: VACALDER
"""

# PROGRAM TO ANALYZE DATA FROM BATCH RUN of NLTHA FOR TDPBEE
#   Victor A Calderon
#   PhD Student/ Research Assistant
#   NC STATE UNIVERSITY 
#   2019 (c)
#
#
# ----------------------------------------------------------------------------
#|                             IMPORTS
# ----------------------------------------------------------------------------

import time
start_time = time.time()
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#-----------------------------------------------------------------------------

# Opening folder to acces data


color=['#D94545','#73A8D4','#F7B76D']
labels=['CL=1.7%','CL=9.8%','CL=13.1%']
legends=[]
MS_path=r'C:\ConditionDependent_PBD\EarthquakeSelection\MainShock_Test'
MSListing = os.listdir(MS_path)
icover=[4]#[4.0,5.0,7.5]
iTcorr= [1.1307]#[1.1307,1.7667,3.975]
iTime=[5,50.,75.] #[5.,10.,15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.]
iwcr= [0.4]#[0.40, 0.45, 0.50, 0.55, 0.60]
rootdir=r'C:\ConditionDependent_PBD\NLTHA_OpenseesPy\data'
clr=-1

Es=29000
earthquake=[]
PGA_MS=[]
covers=[]
times=[]
WaterCement_Ratios=[]
CorrosionLvls_Long=[]
CorrosionLvls_Trans=[]
Steel_Strains=[]
CConc_Strains=[]
UConc_Strains=[]
YieldStresses=[]
YielStressesTrans=[]
AreaOfSteels=[]
spacings=[]
CoreDiameters=[]
AxialLoads=[]
Diameters=[]
AreaRebars=[]
CompStrength=[]
LS_ConcCover=[]
LS_SteelBB=[]
LS_ConfYield=[]
FirstPeriods=[]

#GM2=[r"RSN1231_CHICHI_CHY080-E.AT2"]

for GM in MSListing:
#cover=4.0
    for cover in icover:
        for wcr in iwcr:
#wcr=0.4
            i=-1
            for Time in iTime:
                i=i+1
                clr=clr+1
                datadir=rootdir+"\\"+GM+"\\"+str(cover)+"\\"+str(wcr)+"\\"+str(Time)
                
                
                #Read Conditions
                groundmotion=GM
                pgafile=conditions=open(datadir+"\\PGA.out")
                linespgafile=conditions.readline()
                pga=float(linespgafile.split()[0])
                conditions=open(datadir+"\\Conditions.out")
                linesconditions=conditions.readline()
                
                cov=float(linesconditions.split()[0])
                t=float(linesconditions.split()[1])
                wc=float(linesconditions.split()[2])
                CLl=float(linesconditions.split()[3])
                CLt=float(linesconditions.split()[4])
                
                #Read Period of the Structure
                Period_01=open(datadir+"\\Period.out")
                lines_Period_01=Period_01.readline()
                T1=float(lines_Period_01.split()[0])
                
                #Read Material Properties for run
                
                material_prop=open(datadir+"\\mat.out")
                lines_material_prop=material_prop.readline()
                
                
                YieldStress_Long=float(lines_material_prop.split()[0])
                YieldStress_Trans=float(lines_material_prop.split()[1])
                AreaOfSteel=float(lines_material_prop.split()[2])
                spacing_of_steel=float(lines_material_prop.split()[3])
                CoreDiameter=float(lines_material_prop.split()[4])
                AxialLoad=float(lines_material_prop.split()[5])
                Diameter=float(lines_material_prop.split()[6])
                AreaRebar=float(lines_material_prop.split()[7])
                CompStrengths=float(lines_material_prop.split()[8])
                
                ros=(4*AreaOfSteel)/(CoreDiameter*spacing_of_steel)
                Ag=0.25*math.pi*Diameter**2
                
                e_ccc=0.004
                e_bb=0.03+700*ros*YieldStress_Trans/Es-0.1*AxialLoad/(CompStrengths*Ag)
                e_csy=0.009-0.3*AreaRebar/Ag+3.9*YieldStress_Trans/Es
                
                #Force Displacement Plot
                
                
                d=open(datadir+"\\"+"DFree.out")
                F=open(datadir+"\\"+"RBase.out")
                
            
                linesd = d.readlines()
                linesf = F.readlines()
                x = [line.split()[1] for line in linesd]
                y = [line.split()[1] for line in linesf]
                
                X=[float(i) for i in x]
                Y=[-float(i) for i in y]
                
                plt.figure(1)    
                plt.plot(X,Y,color[clr],label=labels[clr],linewidth=3.0)
#                plt.title('Example for ChiChi EQ w/c=0.4', fontsize=32)
                plt.xlabel('Diplacement (in)', fontsize=60)
                plt.ylabel('BaseShear (kip)', fontsize=60)
                plt.tick_params(direction='out',axis='both',labelsize=48)
                plt.grid()
                plt.legend(fontsize=30)
                plt.show()
                
                
                
                # Steel Stress Strain Analysis
                SteelStressStrain=open(datadir+"\\"+"StressStrain.out")
                linesSteelStressStrain=SteelStressStrain.readlines()
                StlStress=[line.split()[1] for line in linesSteelStressStrain]
                StlStrain=[line.split()[-1] for line in linesSteelStressStrain]
                sigmaStl=[float(i) for i in StlStress]
                epsilonStl=[float(i) for i in StlStrain]
                
                plt.figure(2)
                plt.plot(epsilonStl,sigmaStl,color[clr],label=labels[clr],linewidth=3.0)
                plt.xlabel('strain (in/in)', fontsize=60)
                plt.ylabel('Stress (ksi)', fontsize=60)
                plt.tick_params(direction='out',axis='both',labelsize=48,width=1.5)
                plt.grid()
                plt.legend(fontsize=30)
                plt.show()
                
                # Confined Concrete Stress Strain Analysis
                CConcStressStrain=open(datadir+"\\"+"StressStrain2.out")
                linesCConcStressStrain=CConcStressStrain.readlines()
                CConcStress=[line.split()[1] for line in linesCConcStressStrain]
                CConcStrain=[line.split()[2] for line in linesCConcStressStrain]
                sigmaCConc=[float(i) for i in CConcStress]
                epsilonCConc=[float(i) for i in CConcStrain]
                
#                plt.figure(3)
#                plt.plot(epsilonCConc,sigmaCConc)
#                plt.title('Example Conf Concrete Response for ChiChi EQ w/c=0.4', fontsize=32)
#                plt.xlabel('strain (in/in)', fontsize=24)
#                plt.xlim(-0.015,0)
#                plt.ylabel('Stress (ksi)', fontsize=24)
#                plt.tick_params(direction='out',axis='both',labelsize=20)
#                plt.grid()
#                plt.show()
            
            
                # UncConfined Concrete Stress Strain Analysis
                UnConcStressStrain=open(datadir+"\\"+"StressStrain3.out")
                linesUnConcStressStrain=UnConcStressStrain.readlines()
                UnConcStress=[line.split()[1] for line in linesUnConcStressStrain]
                UnConcStrain=[line.split()[2] for line in linesUnConcStressStrain]
                sigmaUnConc=[float(i) for i in UnConcStress]
                epsilonUnConc=[float(i) for i in UnConcStrain]
                
#                plt.figure(4)
#                plt.plot(epsilonUnConc,sigmaUnConc)
#                plt.title('Example Unconf. Concret Strain Response for ChiChi EQ w/c=0.4', fontsize=32)
#                plt.xlabel('strain (in/in)', fontsize=24)
#                plt.ylabel('Stress (ksi)', fontsize=24)
#                plt.tick_params(direction='out',axis='both',labelsize=20)
#                plt.grid()
#                plt.show()
                earthquake.append(groundmotion)
                PGA_MS.append(pga)
                covers.append(cov)
                times.append(Time)
                WaterCement_Ratios.append(wcr)
                CorrosionLvls_Long.append(CLl)
                CorrosionLvls_Trans.append(CLt)
                Steel_Strains.append(max(max(epsilonStl),abs(min(epsilonStl))))
                CConc_Strains.append(-min(epsilonCConc))
                UConc_Strains.append(-min(epsilonUnConc))
                YieldStresses.append(YieldStress_Long)
                YielStressesTrans.append(YieldStress_Trans)
                AreaOfSteels.append(AreaOfSteel)
                spacings.append(spacing_of_steel)
                CoreDiameters.append(CoreDiameter)
                AxialLoads.append(AxialLoad)
                Diameters.append(Diameter)
                AreaRebars.append(AreaRebar)
                CompStrength.append(-CompStrengths)    
                LS_ConcCover.append(e_ccc)
                LS_ConfYield.append(e_csy)
                LS_SteelBB.append(e_bb)
                FirstPeriods.append(T1)

dataDict={'earthquake':earthquake,
          'pga_(g)':PGA_MS,
          'cover_cm':covers,
          'water_cement_ratio':WaterCement_Ratios,
          'time_yrs':times,
          'CorrosionLvl_Long':CorrosionLvls_Long,
          'CorrosionLvl_Transv':CorrosionLvls_Trans,
          'First_Period_s':FirstPeriods,
          'Steel_Strain':Steel_Strains,
          'Conf_Conc_Strain':CConc_Strains,
          'Unc_Conc_srain':UConc_Strains,
          'Fy_ksi':YieldStresses, 
          'fyt_ksi':YielStressesTrans, 
          'Ast_in2':AreaOfSteels, 
          'st_in':spacings, 
          'Dprime_in':CoreDiameters, 
          'PCol_kip':AxialLoads, 
          'DCol_in':Diameters, 
          'barAreaSec_in2':AreaRebars, 
          'fc_ksi':CompStrength,
          'LimitState_ConcreteCoverCrushing':LS_ConcCover,
          'ConfinementSteelYielding':LS_ConfYield, 
          'LongitudinlSteelBuckling':LS_SteelBB}

DataFrame_Out=pd.DataFrame(dataDict)
#DataFrame_Out.plot.line(x='time_yrs',y='Steel_Strain')#,s=20,c='time_yrs',colormap='viridis')
#plt.title('Confined Concrete maax Strain vs Time',fontsize=32)
#plt.xlabel('Time (yrs)',fontsize=24)
#plt.ylabel('Cnfined Concrete Strain (in/in)',fontsize=24)
#plt.tick_params(direction='out',axis='both',labelsize=20)
#plt.show
DataFrame_Out.to_csv(r'C:\ConditionDependent_PBD\NLTHA_OpenseesPy\Postproc\PosprocData.csv')
print("ALL ANALYSIS COMPLETE")
print("--- %s minutes ---" % ((time.time() - start_time)/60))