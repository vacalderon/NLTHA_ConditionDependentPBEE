# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:45:37 2021

@author: VACALDER
"""
# PROGRAM TO ANALYZE DATA FROM RESPONSE SPECTRUM
#   Victor A Calderon
#   PhD Student/ Research Assistant
#   NC STATE UNIVERSITY 
#   2021 (c)
#
#
# ----------------------------------------------------------------------------
#|                             IMPORTS
# ----------------------------------------------------------------------------
import numpy as np
import pandas as pd
#-----------------------------------------------------------------------------

rootdir=r'C:\ConditionDependentPBEE\NLTHA_ConditionDependentPBEE'
SpectrumDir='C:\ConditionDependentPBEE\GroundmotionSelection\Response Spectrum Analysis\Mainshocks_RS'
ResultsDir='Mainshocks'
results_table=pd.read_csv(rootdir+'\\'+ResultsDir+'\\'+'PostprocData_MS.csv')
SpectralDisplacement_Results=[]
PGD_Results=[]

for GM,row in results_table.iterrows():
    groundmotion=row['earthquake']
    EffectivePeriod=row['Effective period, Teff']
    SpectrumFile = open(SpectrumDir+'\\'+groundmotion+'.csv')
    SpectrumContent=SpectrumFile.readlines()
    SDC=SpectrumContent[12:109]
    SDC_cols=['Period','SD','PSV','PSA']
    SDC_Data=[line.split(',') for line in SDC[:]]
    SDC_DF=pd.DataFrame(columns=SDC_cols,data=SDC_Data)
    PeriodStringList=list(SDC_DF['Period'])
    SpectralDisplacementStringList=list(SDC_DF['SD'])
    PGD=float(SpectralDisplacementStringList[-1])
    T=[float(i) for i in PeriodStringList]
    SpectralDisplacementList=list(SDC_DF['SD'])
    SD_Float=[float(i) for i in SpectralDisplacementList]
    SD_at_Teff=np.interp(EffectivePeriod,T,SD_Float)
    SpectralDisplacement_Results.append(SD_at_Teff)  
    PGD_Results.append(PGD)
    
dataDict = {'SD_at_Teff': SpectralDisplacement_Results,
            'PGD_MS': PGD_Results}
    
    