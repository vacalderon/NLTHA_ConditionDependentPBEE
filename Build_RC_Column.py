# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:12:06 2019

@author: pchi893
"""

from openseespy.opensees import *
# import the os module
#    import os
import math
import numpy as np
#    import matplotlib.pyplot as plt
import ReadRecord
from LibUnitsMUS import *
import ManderCC
from __main__ import *


def Build_RC_Column(dbi,dti,CLl,dblc, cover, Ablc, CLt, Atc, dbtc, datadir,PCol,MS_path,GM):
    # -----------------------------------------------------------------------------
    #
    #    ^Y
    #    |
    #    3       __ 
    #    |          | 
    #    |          |
    #    |          |
    #  (2)       LCol
    #    |          |
    #    |          |
    #    |          |
    #  =2=         _|_  -------->X
    #  =1=      ZeroLength

    # ------------------------------------------------------------------------------
    # |                      IMPORTS
    # ------------------------------------------------------------------------------
    wipe()

    # ------------------------------------------------------------------------------
    #                           GENERATE GEOMETRY
    # ------------------------------------------------------------------------------

    model('basic', '-ndm', 2, '-ndf', 3)
    LCol = 8.0 * ft  # column length
    Weight = PCol  # superstructure weight

    # define section geometry
    DCol = 24.0 * inch  # Column Diameterepth
    

    #PCol = Weight  # nodal dead-load weight per column
    Mass = PCol / g

    ACol = 0.25 * math.pi * DCol ** 2  # cross-sectional area, make stiff
    IzCol = 0.25 * math.pi * DCol ** 4  # Column moment of inertia
    
    node(1, 0.0, 0.0)
    node(2, 0.0, 0.0)
    node(3, 0.0, LCol)
    #  Node, Dx, Dy, Rz
    fix(1, 1, 1, 1)
    fix(2, 1, 0, 0)
    
    mass(3, Mass, 1e-9, 0.0)
    

    
    # MATERIAL parameters
    IDconcC = 1  # material ID tag -- confined cover concrete
    IDconcU = 2  # material ID tag -- unconfined cover concrete
    IDreinf = 3  # material ID tag -- reinforcement
    IDSP    = 4  # material ID tag -- Strain Penetration
    # Define materials for nonlinear columns
    # ------------------------------------------
    # Longitudinal steel properties
    Fy = 60.0 * ksi * (1 - 0.021 * CLl)  # STEEL yield stress
    Fu = 1.375*Fy #Steel Ultimate Stress
    Es = 29000.0 * ksi  # modulus of steel
    Bs = 0.012  # strain-hardening ratio
    R0 = 20.0  # control the transition from elastic to plastic branches
    cR1 = 0.90  # control the transition from elastic to plastic branches
    cR2 = 0.08  # control the transition from elastic to plastic branches
    a1=0.039
    a2=1
    a3=0.029
    a4=1.0
    c = cover * cm  # Column cover to reinforcing steel NA.
    numBarsSec = 16  # number of uniformly-distributed longitudinal-reinforcement bars
    barAreaSec = Ablc * mm2  # area of longitudinal-reinforcement bars
    dbl = dblc * mm
     

    # Transverse Steel Properties
    fyt = 60.0 * ksi * (1 - 0.021 * CLt)  # Yield Stress of Transverse Steel
    Ast = Atc * mm2  # Area of transverse steel
    dbt = dbtc * mm  # Diameter of transverse steel
    st = 2.0 * inch  # Spacing of spiral
    Dprime = DCol - 2 * c - dti*0.5  # Inner core diameter
    #print('Dprime= ',Dprime)
    Rbl = Dprime * 0.5 - dti*0.5 - dbi * 0.5  # Location of longitudinal bar
    #print('Rbl =', Rbl)

    # nominal concrete compressive strength
    fpc = 5.0 * ksi  # CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
    Ec = 57.0 * ksi * math.sqrt(fpc / psi)  # Concrete Elastic Modulus

    # unconfined concrete
    fc1U = -fpc;  # UNCONFINED concrete (todeschini parabolic model), maximum stress
    eps1U = -0.003  # strain at maximum strength of unconfined concrete
    fc2U = 0.2 * fc1U  # ultimate stress
    eps2U = -0.01  # strain at ultimate stress
    lambdac = 0.1  # ratio between unloading slope at $eps2 and initial slope $Ec

    mand = ManderCC.ManderCC(fpc, Ast, fyt, Dprime, st)

    fc = mand[0]
    eps1 = mand[1]
    fc2 = mand[2]
    eps2 = mand[3]

    # CONCRETE                  tag   f'c        ec0   f'cu        ecu
    # Core concrete (confined)
    uniaxialMaterial('Concrete01', IDconcC, fc, eps1, fc2, eps2)

    # Cover concrete (unconfined)
    uniaxialMaterial('Concrete01', IDconcU, fc1U, eps1U, fc2U, eps2U)

    # STEEL
    # Reinforcing steel 
    params = [R0, cR1, cR2]
    #                        tag  fy E0    b
    uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs,R0, cR1, cR2)
    
    
    # STRAIN PENETRATION MATERIAL
    SPalpha=0.4
    SPsy=0.1*((dbl*Fy)*(2*SPalpha+1)/(4000*((-fc)**0.5)))**(1/SPalpha)+0.013
    SPsu=35*SPsy
    SPb=0.45
    SPR=1.01
    
    
    #uniaxialMaterial StrPen01   Tag  fy  sy fu su b  R  
    uniaxialMaterial('Bond_SP01',IDSP, Fy,SPsy,Fu,SPsu,SPb,SPR)
    
    # Writing Material data to file
    with open(datadir+"\\mat.out", 'w') as matfile:
        matfile.write("%s %s %s %s %s %s %s %s %s %s %s %s %s \n" %(Fy, fyt, Ast, st, Dprime, PCol, DCol, barAreaSec, fc,SPsy,SPsu,SPb,SPR))
    matfile.close
    
    #-------------------------------------------------------------------------
    #               DEFINE PLASTICE HIGE PROPERTIES
    #-------------------------------------------------------------------------    
    
    k=0.2*(Fu/Fy - 1)
    if k > 0.08:
        k=0.08
    Leff=LCol
    Lpc=k*Leff + 0.4*DCol
    gamma=0.33 #Assuming unidirectional action
    Lpt=Lpc+gamma*DCol
    # FIBER SECTION properties -------------------------------------------------------------
    # Define cross-section for nonlinear columns
    # ------------------------------------------

    # set some paramaters Section 1
    ColSecTag = 1
    ri = 0.0
    ro = DCol / 2.0
    nfCoreR = 8
    nfCoreT = 8
    nfCoverR = 2
    nfCoverT = 8
    rc = ro - c
    theta = 360.0 / numBarsSec

    section('Fiber', ColSecTag,'-GJ',1e+10)

    # Create the concrete fibers
    patch('circ', 1, nfCoreT, nfCoreR, 0.0, 0.0, ri, rc, 0.0, 360.0)  # Define the core patch
    patch('circ', 2, nfCoverT, nfCoverR, 0.0, 0.0, rc, ro, 0.0, 360.0)  # Define Cover Patch

    # Create the reinforcing fibers
    layer('circ', 3, numBarsSec, barAreaSec, 0.0, 0.0, Rbl, theta, 360.0)
    
    #Set parameters for ZeroLength Element
    
    SecTag2 = 2
    section('Fiber',SecTag2,'-GJ',1e+10)

    # Create the concrete fibers
    patch('circ', 1, nfCoreT, nfCoreR, 0.0, 0.0, ri, rc, 0.0, 360.0)  # Define the core patch
    patch('circ', 2, nfCoverT, nfCoverR, 0.0, 0.0, rc, ro, 0.0, 360.0)  # Define Cover Patch

    # Create the reinforcing fibers
    layer('circ', IDSP, numBarsSec, barAreaSec, 0.0, 0.0, Rbl, theta, 360.0)
    
    # Creating Elements
    
    ColTransfTag = 1
    geomTransf('Linear', ColTransfTag)
    
    ZL_eleTag=1
    element('zeroLengthSection',ZL_eleTag,1,2,SecTag2,'-orient', 0., 1., 0., 1., 0.,0.)
    
    ColeleTag = 2
    
    # Defining Fiber Elements as ForceBeamColumn
    #element('nonlinearBeamColumn', eleTag, 1, 2, numIntgrPts, ColSecTag, ColTransfTag)
    ColIntTag=1
    # beamIntegration('Lobatto',ColIntTag,ColSecTag,numIntgrPts)
    beamIntegration('HingeRadau',ColIntTag,ColSecTag,Lpt,ColSecTag,1e-10,ColSecTag)
    element('forceBeamColumn', ColeleTag, 2, 3, ColTransfTag,ColIntTag,'-mass',0.0)
    
    # Setting Recorders
    
    
    recorder('Node', '-file', datadir + '/DFree.out', '-time','-node', 3, '-dof', 1, 2, 3, 'disp')
    # recorder('Node', '-file', datadir + '/DBase.out', '-time', '-node', 1, '-dof', 1, 2, 3, 'disp')
    recorder('Node', '-file', datadir + '/RBase.out', '-time', '-node', 2, '-dof', 1, 2, 3, 'reaction')
    # recorder('Drift', '-file', datadir+'Data-2c/Drift.out','-time', '-node', 1, '-dof', 1,2,3, 'disp')
    # recorder('Element', '-file', datadir + '/FCol.out', '-time', '-ele', 1, 'globalForce')
    # recorder('Element', '-file', datadir + '/ForceColSec1.out', '-time', '-ele', 1, 'section', 1, 'force')
    recorder('Element', '-file', datadir + '/StressStrain.out', '-time','-ele', 2, 'section', '1', 'fiber', str(Rbl)+', 0.0','mat','3','stressStrain')  #Rbl,0, IDreinf
    recorder('Element', '-file', datadir + '/StressStrain2.out','-time','-ele', 2, 'section', '1', 'fiber', str(-Dprime)+', 0.0','mat','1','stressStrain')  #Rbl,0, IDreinf
    recorder('Element', '-file', datadir + '/StressStrain3.out','-time','-ele', 2, 'section', '1', 'fiber', str(-DCol)+', 0.0','mat','2','stressStrain')
    # recorder('Element', '-file', datadir+'Data-2c/DCol.out','-time', '-ele', 1, 'deformations')
    
    #------------------------------------------------------------------------------ 
    #|                      NLTHA RUN
    #------------------------------------------------------------------------------

    infile=MS_path+"\\"+GM
    outfile=r'C:\ConditionDependent_PBD\EarthquakeSelection\Mainshock_Test_g3files'+"\\"+GM+".g3"
    dt,npt = ReadRecord.ReadRecord(infile,outfile)
    with open(datadir+"\\PGA.out", 'w') as PGAfile:
        accelerations=open(outfile)     
        linesacc = accelerations.readlines()
        acc= [line.split() for line in linesacc]
        flat_list = []
        for sublist in acc:
            for item in sublist:
                flat_list.append(item)
        
        ACC=[float(i) for i in flat_list]
        PGA=max(abs(max(ACC)),abs(min(ACC)))
        PGAfile.write("%s \n" %(PGA))
        PGAfile.close

    #('DCol =',DCol)
    ##import the os module
    #import os
    #import math
    #import numpy as np
    #import matplotlib.pyplot as plt
    #wipe()
    #from LibUnitsMUS import *
    #import ManderCC
    
    #defining gravity loads
    timeSeries('Linear', 1)
    pattern('Plain', 1, 1)
    load(3, 0.0, -PCol, 0.0)
    
    Tol = 1e-8 # convergence tolerance for test
    NstepGravity = 10
    DGravity = 1/NstepGravity
    integrator('LoadControl', DGravity) # determine the next time step for an analysis
    numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
    system('BandGeneral') # how to store and solve the system of equations in the analysis
    constraints('Plain') # how it handles boundary conditions
    test('NormDispIncr', Tol, 6) # determine if convergence has been achieved at the end of an iteration step
    algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
    analysis('Static') # define type of analysis static or transient
    analyze(NstepGravity) # apply gravity
    
    loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
     
    #applying Dynamic Ground motion analysis
    GMdirection = 1
    GMfile = outfile
    GMfact = 1.0
    
    
    
    Lambda = eigen('-fullGenLapack', 2) # eigenvalue mode 1
    Omega = math.pow(Lambda[0], 0.5)
    T1 = 2*np.pi/Omega
    
    
    with open(datadir+"\\Period.out", 'w') as Periodfile:
        Periodfile.write("%s\n" %(T1))
    Periodfile.close
    
    xDamp = 0.05				# 5% damping ratio
    betaKcomm = 2 * (xDamp/Omega)
    alphaM = 0.0				# M-pr damping; D = alphaM*M	
    betaKcurr = 0.0		# K-proportional damping;      +beatKcurr*KCurrent
    betaKinit = 0.0 # initial-stiffness proportional damping      +beatKinit*Kini
    
    rayleigh(alphaM,betaKcurr, betaKinit, betaKcomm) # RAYLEIGH damping
    
    # Uniform EXCITATION: acceleration input
    IDloadTag = 400			# load tag
    Dt = dt			# time step for input ground motion
    GMfatt = GMfact*g			# data in input file is in g Unifts -- ACCELERATION TH
    maxNumIter = 50
    timeSeries('Path', 2, '-dt', Dt, '-filePath', GMfile, '-factor', GMfatt)
    pattern('UniformExcitation', IDloadTag, GMdirection, '-accel', 2) 
    
    wipeAnalysis()
    constraints('Transformation')
    numberer('Plain')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('ModifiedNewton')
    
    NewmarkGamma = 0.5
    NewmarkBeta = 0.25
    integrator('Newmark', NewmarkGamma, NewmarkBeta)
    analysis('Transient')
    analysis_substeps=100
    
    DtAnalysis = dt/analysis_substeps
    TmaxAnalysis = DtAnalysis*analysis_substeps*npt
    
    Nsteps =  int(TmaxAnalysis/ DtAnalysis)
    
    ok = analyze(Nsteps, DtAnalysis)
    
    tCurrent = getTime()
    
    # for gravity analysis, load control is fine, 0.1 is the load factor increment (http://opensees.berkeley.edu/wiki/index.php/Load_Control)
    
    Atest = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    Algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    
    
    for i in Atest:
        
        for j in Algorithm:
            
    
            if ok != 0:
                if j < 4:
                    algorithm(Algorithm[j], '-initial')
                    
                else:
                    algorithm(Algorithm[j])
                    
                test(Atest[i], Tol, 1000)
                ok = analyze(Nsteps, DtAnalysis)                            
                algorithm('ModifiedNewton')
                if ok == 0:
                    print('Analysis succesful: ',Atest[i], ' ', Algorithm[j],' OK = ', ok)
                    
                    break
            else:
                continue
    
    u3 = nodeDisp(3, 1)
    print("GroundMotion Done ", getTime())
    

    
