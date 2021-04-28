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
import pandas as pd


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



icover= [4] #[4.,5.,7.5] #
iTcorr=  [1.1307] #[1.1307,1.7667,3.975] #
iTime= [5.]#[5.,10.,15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75.] #
iwcr= [0.40]#[0.40, 0.45, 0.50, 0.55, 0.60] #
dbi= 0.75   
dti= 0.375  
GM_Path=r'C:\ConditionDependentPBEE\GroundMotion_Squences_Records'
GMListing = os.listdir(GM_Path)
rootdir=r'C:\ConditionDependentPBEE\NLTHA_Sequences'
PCol =225.0*kip
SeqDB = pd.read_csv('"C:\ConditionDependentPBEE\GroundmotionSelection\mainshock_aftershock_file_database.csv"')
# ----------------------------------------------------------------------------
#|                             BATCH RUN
# ----------------------------------------------------------------------------



for GM in GMListing:
    i=-1
    print('GM = ',GM)
    for cover in icover:
        i=i+1
        for Time in iTime:
            for wcr in iwcr:
                #set Functions for Fiber Model and NLTHA
                
                print ('cover is: ', cover,' and Time is:', Time,'and w/c: ',wcr)
                Tcorr=iTcorr[i]
                
                
                dblc  = dbi*25.4-(((1.0508*(1-wcr)**(-1.64))/(cover*10))*(Time-Tcorr)**0.71)
                Ablc  = 0.25*math.pi*dblc**2
                Ablcm = Ablc/(1000.**2)
                Mcorr = Ablcm*7800.
                CLl   = (1-Ablcm*7800./2.223179)*100
                
                
                dbtc  = dti*25.4-(((1.0508*(1-wcr)**(-1.64))/(cover*10))*(Time-Tcorr)**0.71)
                Atc  = 0.25*math.pi*dbtc**2
                Atcm = Atc/(1000.**2)
                CLt   = ((0.55795-Atcm*7800.)/0.55795)*100
                datadir=rootdir+"\\"+"data"+"\\"+GM+"\\"+str(cover)+"\\"+str(wcr)+"\\"+str(Time)
                
                
                
#                print('Tcorr = ',Tcorr)
#                print('Ablc  = ',Ablc)
#                print('Ablcm = ',Ablcm)
#                print('Mcorr = ',Mcorr)
#                print('CLl   = ',CLl )
#                print('dbtc  = ',dbtc)
#                print('Atc  = ',Atc)
#                print('Atcm = ',Atcm)
#                print('CLt   = ',CLt)
                
                if not os.path.exists(datadir):
                    os.makedirs(datadir)

                # Obtain GM dt and npt

                Build_RC_Column.Build_RC_Column(dbi,dti,CLl,dblc,cover,Ablc,CLt,Atc,dbtc,datadir,PCol,GM_file, GM_dt, GM_npt)
                with open(datadir+"\\Conditions.out", 'w') as f:
                    f.write("%s %s %s %s %s \n" %(cover,Time,wcr,CLl,CLt) )
                f.close
                    
                    
print("ALL ANALYSIS COMPLETE")
print("--- %s minutes ---" % ((time.time() - start_time)/60))