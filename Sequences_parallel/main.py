# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:32:05 2019

@author: VACALDER
"""

#------------------------------------------------------------------------------
#|      PROGRAM TO CHECK TIME DEPENDENT PROPERTIES EFFECTS ON STRUCTURES      |
#|      
#|
#|          Victor A Calderon
#|          PhD Student/ Research Assistant
#|          NC STATE UNIVERSITY 
#|          2019 (c)
#|
#|
#------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
#|                             IMPORTS
# ----------------------------------------------------------------------------

#import the os module
import time
start_time = time.time()
import os
import math
from LibUnitsMUS import *
import Build_RC_Column
import Postprocessor_of_data
import pandas as pd
from openseespy.opensees import *

# ----------------------------------------------------------------------------
#| VARIABLES THAT CHANGE WITH TIME
# ----------------------------------------------------------------------------
#
#
# *cover = Cover of concrete in cm
# *Tcorr = Time to corrosion in yrs
# *Time  = Different times that are being analyzed
# *wcr   = Water to cement ratio
# *dbi   = Initial longitudinal bar diameter
# *dti   = Initial transverse steel diameter


compressive_strength_concrete = 5*ksi
yield_strength_long_steel = 60*ksi
yield_strength_trans_steel = 60*ksi
iShapeFactor = [4]
icover = [4] #[4.,5.,7.5] #
iTcorr =  [1.1307] #[1.1307,1.7667,3.975] #
iTime = [5.] #[5.,10.,15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.] #
iwcr = [0.40] #[0.40, 0.45, 0.50, 0.55, 0.60] #
pid=getPID()
np=getNP()
GM_Path = r'.\sequences_records'
GMListing = os.listdir(GM_Path)
rootdir = r'.\\'
iALR = [0.10] #[0.05,0.1,0.15,0.2]
SeqDB = pd.read_csv('focus_ms_as_file_database.csv')
GeomDB = pd.read_csv('column_database_04.csv')
counter = 0
# ----------------------------------------------------------------------------
#|                             BATCH RUN
# ----------------------------------------------------------------------------

for column,Crow in GeomDB.iterrows():
    D=float(Crow['column_diameter'])
    dbi=float(Crow['long_bar_diameter'])
    nbi=float(Crow['number_of_bars_longitudinal'])
    dti=float(Crow['trans_bar_diameter'])
    sti=float(Crow['spacing_trans_steel'])
    rhol=float(Crow['rho_l'])
    rhov=float(Crow['rho_v'])
    for shapefactor in iShapeFactor:
        Height_of_Column=shapefactor*D       
        for ALR in iALR:
            
            Ag=0.25*math.pi*D**2
            AxialLoad=compressive_strength_concrete*Ag*ALR
            
            for GM,row in SeqDB.iterrows():
                i=-1
                GM_fn = row['horizontal_1_filename']
                GM_dt = row['dt_sequence_horizontal1']
                GM_npt = row['npt_sequence_horizontal1']
                print('GM = ',GM_fn)
                GM_file=GM_Path+'\\'+GM_fn
                for cover in icover:
                    i=i+1
                    for Time in iTime:
                        for wcr in iwcr:
                            #set Functions for Fiber Model and NLTHA

                            if (counter % np) == pid :
                            
                                print ('cover is: ', cover,' and Time is:', Time,'and w/c: ',wcr)
                                Tcorr=iTcorr[i]
                                
                                Abi   = 0.25*math.pi*(dbi)**2
                                dblc  = dbi*25.4-(((1.0508*(1-wcr)**(-1.64))/(cover*10))*(Time-Tcorr)**0.71)
                                Ablc  = 0.25*math.pi*dblc**2
                                Ablcm = Ablc/(1000.**2)
                                Mcorr = Ablcm*7800.
                                CLl   = (1-Ablcm/(Abi*0.0254**2))*100
                                
                                
                                Ati   =  0.25*math.pi*(dti)**2
                                dbtc  = dti*25.4-(((1.0508*(1-wcr)**(-1.64))/(cover*10))*(Time-Tcorr)**0.71)
                                Atc  = 0.25*math.pi*dbtc**2
                                Atcm = Atc/(1000.**2)
                                CLt   = (1-Atcm/(Ati*0.0254**2))*100
                                datadir=rootdir+"\\"+"data"+"\\"+GM_fn+"\\"+str(cover)+"\\"+str(wcr)+"\\"+str(Time)+"\\D"+str(D)+"\\SF"+str(shapefactor)+"\\ALR"+str(ALR)+"\\RhoL"+str(rhol)+"\\Rhov"+str(rhov)
                                
                                
                                if not os.path.exists(datadir):
                                    os.makedirs(datadir)
                
    
                                Build_RC_Column.Build_RC_Column(D,Height_of_Column, compressive_strength_concrete,dbi, dti, CLl, dblc, nbi, cover, Ablc, CLt, Atc, dbtc, sti, datadir, AxialLoad, GM_file, GM_dt, GM_npt)
                                with open(datadir+"\\Conditions.out", 'w') as f:
                                    f.write("%s %s %s %s %s \n" %(cover,Time,wcr,CLl,CLt) )
                                f.close
    
                                Postprocessor_of_data.Postprocessor_of_data(GM_fn, cover, wcr, Time,D,shapefactor,ALR,rhol,rhov)
                                os.remove(datadir+"\\StressStrain.out")
                                os.remove(datadir+"\\StressStrain2.out")
                                os.remove(datadir+"\\StressStrain3.out")
                                os.remove(datadir+"\\StressStrain4.out")
                                os.remove(datadir+"\\Conditions.out")
                                os.remove(datadir+"\\DFree.out")
                                os.remove(datadir+"\\mat.out")
                                os.remove(datadir+"\\Period.out")
                                os.remove(datadir+"\\PGA.out")
                                os.remove(datadir+"\\RBase.out")
                                
                                counter += 1
    
                
                
    #                             print('-------------------------------------------------------------------------------------------')
    #                             print("OUTPUT FILES DELETED")
    #                             print('-------------------------------------------------------------------------------------------')
# print("ALL ANALYSIS COMPLETE")
print("--- %s minutes ---" % ((time.time() - start_time)/60))
