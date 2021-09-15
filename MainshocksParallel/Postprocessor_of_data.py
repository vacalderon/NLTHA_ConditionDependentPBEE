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
# |                             IMPORTS
# ----------------------------------------------------------------------------

import time

start_time = time.time()
import os
import math
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import shutil
from LibUnitsMUS import *

# -----------------------------------------------------------------------------

def Postprocessor_of_data(GM_fn, cover, wcr, Time,D,SF,ALR,rhol,rhov):
    # 1. Opening folder to acces data

    SpectrumDir='/share/kowalsky/vacalder/MainshocksParallel/ResponseSpectrumAnalysis'
    rootdir = '/share/kowalsky/vacalder/MainshocksParallel/data'

    Es = 29000
    earthquake = []
    PGA_MS = []
    covers = []
    times = []
    WaterCement_Ratios = []
    CorrosionLvls_Long = []
    CorrosionLvls_Trans = []
    Steel_Strains = []
    CConc_Strains = []
    UConc_Strains = []
    YieldStresses = []
    YielStressesTrans = []
    AreaOfSteels = []
    spacings = []
    CoreDiameters = []
    AxialLoads = []
    Diameters = []
    AreaRebars = []
    CompStrength = []
    LS_ConcCover = []
    LS_SteelBB = []
    LS_ConfYield = []
    FirstPeriods = []
    EffectivePeriods = []
    Forces = []
    Displacements = []
    SpectralDisplacement_Results=[]
    PGD_Results=[]
    Rho_ls=[]
    Rho_vs=[]
    Heights=[]

    datadir = rootdir + '/' + GM_fn + "/" + str(cover) + "/" + str(wcr) + "/" + str(Time)+"/D"+str(D)+"/SF"+str(SF)+"/ALR"+str(ALR)+"/RhoL"+str(rhol)+"/Rhov"+str(rhov)

    # 2. Read Conditions
    groundmotion = GM_fn
    with open(datadir + "/PGA.out") as pgafile:
        linespgafile = pgafile.readline()
    pga = float(linespgafile.split()[0])
    with open(datadir + "/Conditions.out") as  conditions:
        linesconditions = conditions.readline()


    cov = float(linesconditions.split()[0])
    t = float(linesconditions.split()[1])
    wc = float(linesconditions.split()[2])
    CLl = float(linesconditions.split()[3])
    CLt = float(linesconditions.split()[4])


    # 3. Read Period of the Structure
    with open(datadir + "/Period.out") as  Period_01:
        lines_Period_01 = Period_01.readline()
    T1 = float(lines_Period_01.split()[0])

    
    # 4. Read Material Properties for run

    with open(datadir + "/mat.out") as material_prop:
        lines_material_prop = material_prop.readline()

    YieldStress_Long = float(lines_material_prop.split()[0])
    YieldStress_Trans = float(lines_material_prop.split()[1])
    AreaOfSteel = float(lines_material_prop.split()[2])
    spacing_of_steel = float(lines_material_prop.split()[3])
    CoreDiameter = float(lines_material_prop.split()[4])
    AxialLoad = float(lines_material_prop.split()[5])
    Diameter = float(lines_material_prop.split()[6])
    Height = float(lines_material_prop.split()[7])
    AreaRebar = float(lines_material_prop.split()[8])
    CompStrengths = float(lines_material_prop.split()[9])

    ros = (4 * AreaOfSteel) / (CoreDiameter * spacing_of_steel)
    Ag = 0.25 * math.pi * Diameter ** 2

    e_ccc = 0.004
    e_bb = 0.03 + 700 * ros * YieldStress_Trans / Es - 0.1 * AxialLoad / (CompStrengths * Ag)
    e_csy = 0.009 - 0.3 * AreaRebar / Ag + 3.9 * YieldStress_Trans / Es


    # 5. Force Displacement Plot

    with open(datadir + "/DFree.out") as d:
        linesd = d.readlines()
    with open(datadir + "/RBase.out") as F:
        linesf = F.readlines()

    x = [line.split()[1] for line in linesd[:-1]]
    y = [line.split()[1] for line in linesf[:-1]]

    X = [float(i) for i in x]
    Y = [-float(i) for i in y]
    maxDisp = max(X)
    minDisp = min(X)
    if maxDisp>abs(minDisp) :
        AbsMaxDisp=maxDisp
    elif maxDisp<abs(minDisp) :
        AbsMaxDisp=minDisp
    maxDispPoss = X.index(AbsMaxDisp)
    maxForce_at_maxDisp = Y[maxDispPoss]
    Keff = abs(maxForce_at_maxDisp)/abs(AbsMaxDisp)
    meff = 225.0*kip/g
    Teff = (2*math.pi)*(math.sqrt((meff/Keff)))


    # 6. Steel Stress Strain Analysis
    with open(datadir + "/StressStrain.out") as SteelStressStrain:
        linesSteelStressStrain = SteelStressStrain.readlines()
    StlStress = [line.split()[1] for line in linesSteelStressStrain]
    StlStrain = [line.split()[-1] for line in linesSteelStressStrain]
    siGM_fnaStl = [float(i) for i in StlStress]
    epsilonStl = [float(i) for i in StlStrain[:-1]]

    # 7. Confined Concrete Stress Strain Analysis
    with open(datadir + "/StressStrain2.out") as CConcStressStrain:
        linesCConcStressStrain = CConcStressStrain.readlines()
    CConcStress = [line.split()[1] for line in linesCConcStressStrain]
    CConcStrain = [line.split()[2] for line in linesCConcStressStrain]
    siGM_fnaCConc = [float(i) for i in CConcStress]
    epsilonCConc = [float(i) for i in CConcStrain[:-1]]

    # 8. UncConfined Concrete Stress Strain Analysis
    with open(datadir + "/StressStrain3.out") as UnConcStressStrain:
        linesUnConcStressStrain = UnConcStressStrain.readlines()
    UnConcStress = [line.split()[1] for line in linesUnConcStressStrain]
    UnConcStrain = [line.split()[2] for line in linesUnConcStressStrain]
    siGM_fnaUnConc = [float(i) for i in UnConcStress]
    epsilonUnConc = [float(i) for i in UnConcStrain[:-1]]
    
    # 9. Writing SD_teff
    SpectrumFile = open(SpectrumDir+'/'+groundmotion+'.csv')
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
    SD_at_Teff=np.interp(Teff,T,SD_Float)
    
    
    # 10. Writing data to variables 
    
    earthquake.append(groundmotion)
    PGA_MS.append(pga)
    covers.append(cov)
    times.append(Time)
    WaterCement_Ratios.append(wcr)
    CorrosionLvls_Long.append(CLl)
    CorrosionLvls_Trans.append(CLt)
    Steel_Strains.append(max(max(epsilonStl), abs(min(epsilonStl))))
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
    EffectivePeriods.append(Teff)
    Forces.append(maxForce_at_maxDisp)
    Displacements.append(AbsMaxDisp)
    SpectralDisplacement_Results.append(SD_at_Teff)  
    PGD_Results.append(PGD)
    Rho_ls.append(rhol)
    Rho_vs.append(rhov)
    Heights.append(Height)
    
    # 10. Preparing dictionary to wirte output database
    
    dataDict = {'earthquake': earthquake,
                'pga_(g)': PGA_MS,
                'cover_cm': covers,
                'water_cement_ratio': WaterCement_Ratios,
                'time_yrs': times,
                'CorrosionLvl_Long': CorrosionLvls_Long,
                'CorrosionLvl_Transv': CorrosionLvls_Trans,
                'First_Period_s': FirstPeriods,
                'Steel_Strain': Steel_Strains,
                'Conf_Conc_Strain': CConc_Strains,
                'Unc_Conc_srain': UConc_Strains,
                'Fy_ksi': YieldStresses,
                'fyt_ksi': YielStressesTrans,
                'Ast_in2': AreaOfSteels,
                'st_in': spacings,
                'Dprime_in': CoreDiameters,
                'PCol_kip': AxialLoads,
                'DCol_in': Diameters,
                'barAreaSec_in2': AreaRebars,
                'fc_ksi': CompStrength,
                'LimitState_ConcreteCoverCrushing': LS_ConcCover,
                'ConfinementSteelYielding': LS_ConfYield,
                'LongitudinlSteelBuckling': LS_SteelBB,
                'Effective period, Teff' : EffectivePeriods,
                'Force' : Forces,
                'MaxDisplacement at MaxForce' : Displacements,
                'SD_at_Teff': SpectralDisplacement_Results,
                'rhol' :Rho_ls,
                'rhov' : Rho_vs,
                'height_of_col': Heights}
    
    # 11. Generating data frame to write data to csv file
    
    DataFrame_Out = pd.DataFrame(dataDict)
    
    # 12. Writing CSV File
    DataFrame_Out.to_csv('/share/kowalsky/vacalder/MainshocksParallel/results/PosprocData.csv', mode='a', header=False)
    
    # Output to show in console
    print('-----------------------------------------------------------------------------------------------------------')
    print("POSTPROCESSING COMPLETE")
    print('-----------------------------------------------------------------------------------------------------------')
    print("--- %s minutes ---" % ((time.time() - start_time) / 60))
