# *****************************************
# **  author: Stefano Magni              **
# **  email: magnistefano01[AT]gmail.com **
# **                                     **
# **  first created: 11.01.2016          **
# *****************************************

# This code implements the dynamical model of the Heat Shock Response in Chlamydomonas reinhardtii
# developed by Stefano Magni, Antonella Succurro, Alexander Skupin, and Oliver Ebenhoeh. 
# This model and the associated scientific results are discussed in the following paper: XXX
# This code can be downloaded from: XXX

import os

from HSM_SimulateClass import *
from HSM_StudyEquilibrium import *
from HSM_StudyHPproduction import *
from HSM_VaryParamsRMSvsData import *
from HSM_StudyMaxPdFuncOfTau import *
from HSM_calibrationRMSmainFunctions import *


if __name__ == '__main__':

    FigureExtension = ".pdf"  # Choose among .eps, .pdf, .png,... the extension to save all figures

    # Select which simulations you want to launch
    # PART I: Calibration and RMS stuff
    ComputeRMSfiducialParamSet = "No"
    MyNumberOfRandomSets = 2
    MCrandomScanParamSpaceRMS = "Yes"
    Part1_GenerateMCsamples = "No"
    Part2_PlotRMSfromMC = "Yes"
    Part3_ComputeRMSvsDoubleHS = "No"
    GradientSearchRMS = "No"
    ComputeRMSfinalParamSet = "Yes"
    # PART II: Simulations of HSR and Data comparison (main figures paper)
    HS_Response = "Yes"
    Time_Course_Data = "Yes"
    Schmollingen_Data = "Yes"
    Plots3Dfeeding = "No"
    Double_Heat = "Yes"
    Double_Heat_ARS = "Yes"
    Early_Late_24h_HS = "Yes"
    Hot_Day_Behaviour = "Yes"
    HPproductionForVariousTempDurationsHS = "Yes"
    SteadyStateSystematicStudy = "Yes"
    Study_Stability_Short = "Yes"
    UnfoldedPfuncOfTau = "Yes"

    # Here specify the names of the folders to be used after, mainly for the calibration
    FolderContainingCsvFiles = "CsvFilesWithRMSandParamSets/"
    FolderForAllCalibrationPlots = "PlotsCalibrationRMS/"
    FolderContainingDataVsSimuCalibration = FolderForAllCalibrationPlots + "PlotsDataVsSimuCalibration/"
    FolderContaining1ParametrsRMSplots = FolderForAllCalibrationPlots + "PlotsRMS1Parameter/"
    FolderContaining2ParametrsRMSplots = FolderForAllCalibrationPlots + "PlotsRMS2Parameters/"
    FolderContainingGradientSearchPlots = FolderForAllCalibrationPlots + "PlotsGradientSearch/"
    FolderRMS1vs2 = FolderForAllCalibrationPlots + "PlotRMS1vs2/"

    if not os.path.exists(FolderContainingCsvFiles):
        os.makedirs(FolderContainingCsvFiles)
    if not os.path.exists(FolderContainingDataVsSimuCalibration):
        os.makedirs(FolderContainingDataVsSimuCalibration)
    if not os.path.exists(FolderContaining1ParametrsRMSplots):
        os.makedirs(FolderContaining1ParametrsRMSplots)
    if not os.path.exists(FolderContaining2ParametrsRMSplots):
        os.makedirs(FolderContaining2ParametrsRMSplots)
    if not os.path.exists(FolderContainingGradientSearchPlots):
        os.makedirs(FolderContainingGradientSearchPlots)
    if not os.path.exists(FolderRMS1vs2):
        os.makedirs(FolderRMS1vs2)

    DefaultParamSetInitCond = {  # Set initial conditions of ODEs system
        "Pin": 100000.,          # (microM) protein P
        "Phin": 1.,              # (microM) deg. protein P#
        "Sin": 0.1,              # (microM) stresskinease S
        "Ssin": .1 / 20.,        # (microM) phosphor stresskinease S*
        "Fin": 10.5,             # (microM) HSF F
        "Fsin": 1.,              # (microM) phosphor. HSF F*
        "Gin": 0.0012,           # (microM) free Gene G
        "FsGin": 0.0002,         # (microM) active Gene F*G
        "FGin": 0.0008,          # (microM) inactive Gene FG
        "RFin": 0.0036,          # (microM) mRNA_F
        "RHPin": 0.00360,        # (microM) mRNA_HP
        "HPin": 1.}              # (microM) heatshock protein HP

    # IMPORTANT: THE FOLLOWING PARAMETER SET IS THE ONE USED AS A STARTING POINT FOR THE CALIBRATION PROCEDURE. 
    # THE SIMULATIONS WILL USE THE PARAMETER SET PROVIDED BY THE TXT FILE
    DefaultParamSetRATES = {  # Set values for the rates
        "kP0": 10.,  # //0.5canav!         # P#-->P       ((microM s)^-1)
        "kP0p": 100.,  # 70.,              # P-->P#       (s^-1)
        "kS": 100.,                        # S*-->S       (s^-1)
        "kSp0": 500.,                      # S-->S*       (s^-1)
        "kFp0": 1.,                        # F*-->F       (s^-1)          
        "kF0": 1.,                         # F-->F*       (s^-1) 
        "kFpi0": 0.02,  # //0.01           # mRF: F       (s^-1)          
        "kFGp": 0.10,                      # FG-> F +G    (s^-1)          
        "kFG": 0.0050,                     # F +G -> FG   ((microM s)^-1) 
        "ketaF": 0.001,                    # F-->         (s^-1)          
        "kFsG": 1.0,                       # F* + G->F*G  ((microM s)^-1) 
        "kFsGp": 0.50,                     # F*G->F* + G  (s^-1)         
        "kFsp": 0.010,                     # F*G->FG      (s^-1)         
        "kFs": 0.010,                      # FG->F*G      (s^-1)         
        "kpiRF": 2.0*8.,                   # F*G: mRF     (s^-1)
        "kpiRH": 9./2.,                    # F*G: mRHP    (s^-1)
        "kpiHP": 0.5,                      # mRHP: HP     (s^-1)
        "ketaHP": 1.72*0.00005, # //0.0001 # HP-->        (s^-1)          
        "ketaRF": 0.006/4.,                # mRF-->       (s^-1)         
        "ketaRHP": 0.006/5.}               # mRHP-->      (s^-1)         

    DefaultParamSetForREACTIONS = {
        "n1": 10,                     # (adimensional)
        "T0const": 36.,               # (deg C) Threshold Temperature
        "n2": 10./2.,                 # (adimensional) n2 = 10 in Alex's code, 1 in his draft... Nonlinear protein degradation
        "P0const": 600.,              # (microM)
        "I": 1.,                      # (microM) iso...kinease # IS I NOT VARYING?
        "piRFconst": DefaultParamSetRATES["ketaRF"]*DefaultParamSetRATES["ketaF"]/DefaultParamSetRATES["kFpi0"]*0.02125,
                                      # Basal production of RF, added to Alex's model [CHANGE IT IF n2 CHANGES!!!]
        "piRHPconst": DefaultParamSetRATES["ketaRHP"]*DefaultParamSetRATES["ketaHP"]/DefaultParamSetRATES["kpiHP"]*17.5,}
                                      # Basal production of RHP, added to Alex's model [CHANGE IT IF n2 CHANGES!!!]

    ################################################################################################################
    ################################################################################################################
    ######################### PART I: Play around to chose the parameter set for rates #############################
    ################################################################################################################
    ################################################################################################################

    ########################################################################
    #####################    BIG 0: Preliminaries     ######################
    ########################################################################

    StartingParamSetRATES = DefaultParamSetRATES

    ############ Fix REACTIONS PARAMETERS SET
    TestParamSetForREACTIONS = deepcopy(DefaultParamSetForREACTIONS) 

    ############# 1) ########### EXTRACT EXPERIMENTAL DATA USED FOR RMS Feeding and RMS Double HS
    ##### FEEDING EXPERIMENTS - FROM FILES ALEXANDER SKUPIN,
    AllDataControlsFeeding = ExtractDataControlsFeedingExperimentsFromFilesIntoListOfDictionaries()
    ##### Double HS from files from figure 7 Schroda et al. 2000 (for later use)
    AllDataControlsDoubleHS = ExtractDataControlsDoubleHSExperimentFromFilesIntoListOfDictionaries()

    ############ 2) ############ COMPUTE RMS w.r.t. double HS for the FIDUCIAL PARMETER SET
    if ComputeRMSfiducialParamSet == "Yes":
        
        RMSFeeding = ComputeRMSfeedingForGivenParameterSet(DefaultParamSetRATES, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsFeeding)
        RMSdoubleHS = ComputeRMSdoubleHSforGivenParameterSet(DefaultParamSetRATES, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsDoubleHS)
        print("\nFor the FIDUCIAL PARAMETER SET, RMS w.r.t. Feeding is " + str(RMSFeeding) + " and RMS w.r.t. Double HS is " + str(RMSdoubleHS) + "\n")


    if MCrandomScanParamSpaceRMS == "Yes":

        #################################################################################################################################
        ################ BIG 1: Generate random (MC) params sets or change params 1by1, and compute RMS Feeding #########################
        #################################################################################################################################
        # A switch to switch between "random" and "1by1" for changing the values of the parameters !!!
        SwitchRandomSetsOrParametersK1by1Sets = "RandomSets" # "ParametersK1by1" "RandomSets" DO BOTH!!!!!!!

        if Part1_GenerateMCsamples == "Yes":
            print("\nSTARTING MC RANDOM SCAN METER SPACE...\n")
        
            # If "RandomSets" then the following parameters will be used:
            NumberOfRandomSets = MyNumberOfRandomSets # we used up to 100000
            FactorOfRandom = 0.5 # 0.5 means 50% variation of the parameter
            # If "ParametersK1by1" then the following parameters will be used:
            NumberOfValuesForEachParameterk = 3 # we used 100 
            FactorOfK1by1 = 0.5 # 0.5 means 50% variation of the parameter

            NameOfOutputFileRMSmanyParamsSets = 'OutputFileRMSmanyParamsSets.csv'
            NameOfOutputFileKeys = 'OutputFileKeys.csv'
    
            GenerateMCRandomOrNotParSetsAndComputeRMSFeeding(SwitchRandomSetsOrParametersK1by1Sets, NumberOfRandomSets, NumberOfValuesForEachParameterk, FactorOfRandom, FactorOfK1by1, FolderContainingCsvFiles, FolderContainingDataVsSimuCalibration, NameOfOutputFileRMSmanyParamsSets, NameOfOutputFileKeys, StartingParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, FigureExtension)

        
        ##############################################################################################################
        #####################    BIG 2: Plot RMS values as function of parameters, from file     #####################
        ##############################################################################################################
        if Part2_PlotRMSfromMC == "Yes":
            print("\nSTARTING TO PLOT RMS VALUES...\n")

            #FileNameManyParamsSetsRMS = FolderContainingCsvFiles + 'OutputFileRMSmanyParamsSets.csv'
            FileNameManyParamsSetsRMS = FolderContainingCsvFiles + 'InterestingListsOfParamsSets/OutputFile100000.csv'#'aaa.csv'
            # FolderContainingCsvFiles + 'OutputFileRMSmanyParamsSets.csv'                           
            # 'InterestingListsOfParamsSets/OutputFile100000.csv'                                    
            # 'InterestingListsOfParamsSets/OutputFileKparams1by1.csv'                               

            FileNameKeysNamesParamsSets = FolderContainingCsvFiles + 'InterestingListsOfParamsSets/OutputFileKeys100000.csv'#'OutputFileKeys.csv'          
            # FolderContainingCsvFiles + 'OutputFileKeys.csv'                                      
            # 'InterestingListsOfParamsSets/OutputFileKeys100000.csv'                               
            # 'InterestingListsOfParamsSets/OutputFileKeys1by1.csv'    

            NumberOfBestRMSparamsSetsPlotted = 300#MyNumberOfRandomSets 
                             
            PlotRMSvaluesAsFunctionOfParametersFromFile(FolderContainingCsvFiles, FolderContaining1ParametrsRMSplots, FolderContaining2ParametrsRMSplots, FileNameManyParamsSetsRMS, FileNameKeysNamesParamsSets, NumberOfBestRMSparamsSetsPlotted, StartingParamSetRATES, SwitchRandomSetsOrParametersK1by1Sets, FigureExtension, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding)

        ######################################################################################################################################
        ##############   BIG 3: For the best 5000 points wrt RSM Feeding, compute the corresponding RMS w.r.t. Double HS   ###################
        ######################################################################################################################################
        if Part3_ComputeRMSvsDoubleHS == "Yes":
            print("\nSTART TO COMPUTE RMS DOUBLE HS FOR ALL BEST 5000 POINTS FEEDING...\n")

            NameFigureRMS1vs2 = "RMS1vs2"

            FileNameManyParamsSetsRMSDoubleHS = FolderContainingCsvFiles + 'OutputFileRMSmanyParamsSets.csv'                       
            # FileNameManyParamsSetsRMSDoubleHS = FolderContainingCsvFiles + 'OutputFileRMSmanyParamsSets.csv'                     
            # FileNameManyParamsSetsRMSDoubleHS = 'InterestingListsOfParamsSets/OutputFileKparams1by1.csv'                         
            # FileNameManyParamsSetsRMSDoubleHS = FolderContainingCsvFiles + 'CutORDEREDOutputFileRMSmanyParamsSetsMOD5000.csv'    
            # FileNameManyParamsSetsRMSDoubleHS = 'InterestingListsOfParamsSets/CutORDEREDOutputFileRMSmanyParamsSetsMOD5000.csv'  

            ForBestRMSfeedingPointspointsComputeRMSDoubleHS(FolderRMS1vs2, FolderContainingCsvFiles, NameFigureRMS1vs2, FileNameManyParamsSetsRMSDoubleHS, StartingParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsDoubleHS, FigureExtension)


    if GradientSearchRMS == "Yes":

        #####################################################################
        ############## BIG 4: Launch the gradient search ####################
        #####################################################################

        print("\nSTARTING GRADIENT SEARCH...\n")

        ### Choose if for the minimization you want to use the RMS w.r.t. the controls of the feeding, or RMSfeedings + RMSdoubleHS
        UseRMSForFeedingOrTotal = "Feeding" # "FeedingPlusDouble"

        MaxNumberOfIterations = 150
        NumberOfIterationsForAverage = 10
        ThresholdAverageRMSdecrease = 0.00003 # 0.00001
        IncrementInComputingDerivative = 1.e-6

        GammaMin = 0.
        GammaMax = 5.
        GammaBisectionStep = 0.01

        NameOutputFileBestParametersSet = 'OutputFileBestParametersSet.csv'

        ExecuteGradientSearch(FolderContainingGradientSearchPlots, UseRMSForFeedingOrTotal, MaxNumberOfIterations, NumberOfIterationsForAverage, ThresholdAverageRMSdecrease, IncrementInComputingDerivative, NameOutputFileBestParametersSet, StartingParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, GammaMin, GammaMax, GammaBisectionStep, FigureExtension)


    ##########################################################################################################################################
    ############## BIG 5: Charge the parameter set for rates from the file containing the one selected by the gradient search ################
    ##########################################################################################################################################
    print("\nSTARTING TO CHARGE PARAMTER SET FROM FILE..\n")
    # Read the whole file into a variable which is a list of every row of the file.
    Datafile = open('OutputFileBestParametersSet.csv', 'r')
    DataLines = Datafile.readlines()
    Datafile.close()

    # Initialize the lists which will contain the data:
    BestParameterSetFromGradientSearchFromFile = {}

        # Scan the rows of the file stored in lines, and put the values into some variables:
    for line in DataLines:
        SplittedLine = line.split()
        Key = str(SplittedLine[0])
        Value = float(SplittedLine[1])
        BestParameterSetFromGradientSearchFromFile.update({ Key : Value })

    print("The parameter set charged from the txt file and used for all the simulations from now on is:")
    print(BestParameterSetFromGradientSearchFromFile)

    
    ################################################################################################################
    ################################################################################################################
    ####################################### PART II: Now do the real stuff!!! ######################################
    ################################################################################################################
    ################################################################################################################
    print("\nSTART RUNNING SIMULATIONS OF PAPER...\n")
    ############ 1) ############ Plot Controls data / Double HS data vs model to see if CALIBRATION was ok
    ############ 2) ############ Plot other data vs model to see if VALIDATION is ok
    ############ 3) ############ See the TYPICAL BEHAVIOUR of the system
    ############ 4) ############ Use the model to SIMULATE FURTHER situations
    ############ 5) ############ Verify stability of the system and study steady state

    # Set default parameter values
    #(BestRMSDictionaryParameters)
    #(DefaultParamSetRATES)
    #(FiducialParameterSetRates)
    #(StartingParameterSetKs)
    #(FinalParamSetMinimizingRMS)
    #(ORIGINALFinalParamSetMinimizingRMS)
    ThisParametrSet = deepcopy(BestParameterSetFromGradientSearchFromFile)

    if ComputeRMSfinalParamSet == "Yes":
        RMSFeedingList = ComputeRMSfeedingForGivenParameterSet(ThisParametrSet, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsFeeding)
        RMSFeeding = RMSFeedingList[0]
        RMSdoubleHS = ComputeRMSdoubleHSforGivenParameterSet(ThisParametrSet, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsDoubleHS)
        print("\nFor the CURRENT PARAMETER SET, RMS w.r.t. Feeding is " + str(round(RMSFeeding,5)) + " and RMS w.r.t. Double HS is " + str(round(RMSdoubleHS,5)) + "\n")
        PlotResultOfBestFitToData(FolderContainingDataVsSimuCalibration, ThisParametrSet, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, FigureExtension)

    MyParamSetRATES = ParametersSet(ThisParametrSet)
    MyParamSetIC = ParametersSet(DefaultParamSetInitCond)
    MyParamSetForREACTIONS = ParametersSet(DefaultParamSetForREACTIONS)

    # Create an object of the class "Heat shock models" with these parameters
    MyHSM = HeatShockModel(MyParamSetIC, MyParamSetRATES, MyParamSetForREACTIONS)

    #################################################################
    ###############  Simulate HEAT SHOCK RESPONSE  ##################
    #################################################################

    if HS_Response == "Yes":
        print("\nSTARTING TO SIMULATE HEAT SHOCK RESPONSE...\n")
        # Set parameters for the Temperature
        TsetHSR = ParametersSet({"Ttype": 1, "Tin": 25., "Tup": 42., "tau": 5., "ta": 20. * 60. + vorl})
        # Set parameters for the Time
        timesetHSR = ParametersSet({"t_start": 0., "t_stop": 140. * 60. + vorl, "delta_t": 5.0})
        # Run a simulation on the model, using Temperature and time settings defined above
        # and giving a name to label the output plots
        SimulationHSR = Simulate(MyHSM, timesetHSR, TsetHSR, "SimulationHSResponse" + FigureExtension)
        SimulationHSR.TimeRun(ZoomInPanelA="Yes3")


    ################################################################################################
    ############# Simulate time course experiment Muehlhaus 2011 and compare with DATA #############
    ################################################################################################

    if Time_Course_Data == "Yes":
        print("\nSTARTING TO Simulate time course experiment Muehlhaus 2011 and compare with DATA...\n")
        TsetHSRdata = ParametersSet({"Ttype": 1, "Tin": 25., "Tup": 42., "tau": 5., "ta": 0. * 60. + vorl})
        timesetHSRdata = ParametersSet({"t_start": 0., "t_stop": 185. * 60. + vorl, "delta_t": 5.0})
        SimulationHSRdata = Simulate(MyHSM, timesetHSRdata, TsetHSRdata, "SimulationHSResponseVsData" + FigureExtension)

        SimulationHSRdata.TimeCourseVsDataPlot()

        SimulationHSRdata.TimeCourseVsDataPlotAllInOne()

    ##############################################################################################
    ############   Simulate experiments Schmollingen 2013 and compare with DATA    ###############
    ##############################################################################################

    if Schmollingen_Data == "Yes":
        print("\nSTARTING TO Simulate experiments Schmollingen 2013 and compare with DATA...\n")
        TsetSchmol2013data = ParametersSet({"Ttype": 1, "Tin": 25., "Tup": 40., "tau": 5., "ta": 0. * 60. + vorl})
        Ylegend = r'Concentration of mRNA$_{F}$ ($\%$ of maximum)'
        NameOfFigure = "RF"

        ############ Simulate feeding with STAUROSPORINE experiment and compare with DATA  ###########
        # Settings
        timesetSTAURdata = ParametersSet({"t_start": 0., "t_stop": 60. * 60. + vorl, "delta_t": 5.0})
        SimulationSTAURdata = Simulate(MyHSM, timesetSTAURdata, TsetSchmol2013data,
                                       "SimulationSTAUR2data" + FigureExtension)
        ListOfKvaluesSTAUR = [ThisParametrSet["kFp0"], ThisParametrSet["kFp0"]/100.*60., ThisParametrSet["kFp0"]/100.*10.]
        ListOfKvaluesSTAURMOD = []
        for i in range(len(ListOfKvaluesSTAUR)):
            #ListOfKvaluesSTAURMOD.append(round(ListOfKvaluesSTAUR[i] * 100. / max(ListOfKvaluesSTAUR)))
            ListOfKvaluesSTAURMOD.append(round(ListOfKvaluesSTAUR[i] * 100. / ThisParametrSet["kFp0"]))
        DataFileName = "DataFiles/DataSchmol2013StaurFig1B1.dat"

        # Plot settings
        ModelLegend = [r"k$_F$' = %r" % ListOfKvaluesSTAURMOD[0] + "% of k$_F$' nominal", r"k$_F$' = %r"
                       % ListOfKvaluesSTAURMOD[1] + "% of k$_F$' nominal",
                       r"k$_F$' = %r" % ListOfKvaluesSTAURMOD[2] + "% of k$_F$' nominal"]
        DataLegend = [r"Control", r"+ Staur $20$ nM", r"+ Staur $1$ $\mu$M"]
        LegendPosition = "upper right"

        # Simulate
        SimulationSTAURdata.FeedingExperimentPlotsVsData("kFp0", ListOfKvaluesSTAUR,
                                                         ModelLegend, SimulationSTAURdata.RF, DataFileName, DataLegend,
                                                         Ylegend, NameOfFigure, LegendPosition, 4)
        SimulationSTAURdata.FeedingExperimentPlotsVsDataALLinONE("kFp0", ListOfKvaluesSTAUR,
                                                         ModelLegend, SimulationSTAURdata.RF, DataFileName, DataLegend,
                                                         Ylegend, NameOfFigure, LegendPosition, 4)
        if Plots3Dfeeding == "Yes":
            SimulationSTAURdata.FeedingExperiment3Dplots(1, "kFp0", r"k$_F$' (s$^{-1}$)", 0.1, 1., 20, [23., 283.])


        #############  Simulate feeding with RADICICOL experiment and compare with DATA   ############
        # Settings
        timesetRADdata = ParametersSet({"t_start": 0., "t_stop": 180. * 60. + vorl, "delta_t": 5.0})
        SimulationRADdata = Simulate(MyHSM, timesetRADdata, TsetSchmol2013data, "SimulationRADdata" + FigureExtension)
        ListOfKvaluesRAD = [ThisParametrSet["kP0"], ThisParametrSet["kP0"]/100.*60., ThisParametrSet["kP0"]/100.*30.]
        ListOfKvaluesRADMOD = []
        for i in range(len(ListOfKvaluesRAD)):
            ListOfKvaluesRADMOD.append(round(ListOfKvaluesRAD[i] * 100. / ThisParametrSet["kP0"]))
        DataFileName = "DataFiles/DataSchmol2013StaurFig4B1.dat"

        # Plot settings
        ModelLegend = [r"$k_{P}$ = %r"
                       % ListOfKvaluesRADMOD[0] + "% of $k_{P}$ nominal", r"$k_{P}$ = %r"
                       % ListOfKvaluesRADMOD[1] + "% of $k_{P}$ nominal",
                       r"$k_{P}$ = %r" % ListOfKvaluesRADMOD[2] + "% of $k_{P}$ nominal"]
        DataLegend = [r"Control", r"+ Radicicol $10$ $\mu M$", r"+ Radicicol $100$ $\mu M$"]
        LegendPosition = "upper right"

        # Simulate
        SimulationRADdata.FeedingExperimentPlotsVsData("kP0", ListOfKvaluesRAD,
                                                       ModelLegend, SimulationRADdata.RF, DataFileName, DataLegend,
                                                       Ylegend, NameOfFigure, LegendPosition, 4)
        SimulationRADdata.FeedingExperimentPlotsVsDataALLinONE("kP0", ListOfKvaluesRAD,
                                                       ModelLegend, SimulationRADdata.RF, DataFileName, DataLegend,
                                                       Ylegend, NameOfFigure, LegendPosition, 4)
        if Plots3Dfeeding == "Yes":
            SimulationRADdata.FeedingExperiment3Dplots(1, "kP0", r"k$_P$' (($\mu$Ms)$^{-1}$)", 0.5, 30, 20, [18., 100.])

    #######################################################################
    #############  Simulate DOUBLE HEAT SHOCK EXPERIMENT  #################
    #######################################################################

    if Double_Heat == "Yes":
        print("\nSTARTING TO Simulate DOUBLE HEAT SHOCK EXPERIMENT...\n")

        Tin = 25.
        Tup = 42.

        ### 30' ###
        Tset30 = ParametersSet(
                {"Ttype": 3, "Tin": Tin, "Tup": Tup, "tau": 5., "ta": 40. * 60. + vorl, "tb": 160. * 60. + vorl,
                 "tc": 190. * 60. + vorl})
        timeset30 = ParametersSet({"t_start": 0., "t_stop": 1000. * 60. + vorl, "delta_t": 5.0})
        Simulation30 = Simulate(MyHSM, timeset30, Tset30, "SimulationDoubleHeat30min" + FigureExtension)
        Simulation30.TimeRun(ZoomInPanelA="Yes3")

        ### 2h ###
        Tset2h = ParametersSet(
                {"Ttype": 3, "Tin": Tin, "Tup": Tup, "tau": 5., "ta": 40. * 60. + vorl, "tb": 160. * 60. + vorl,
                 "tc": 280. * 60. + vorl})
        timeset2h = ParametersSet({"t_start": 0., "t_stop": 1000. * 60. + vorl, "delta_t": 5.0})
        Simulation2h = Simulate(MyHSM, timeset2h, Tset2h, "SimulationDoubleHeat2h" + FigureExtension)
        Simulation2h.TimeRun(ZoomInPanelA="Yes3")

        ### 3.5h ###
        Tset35h = ParametersSet(
                {"Ttype": 3, "Tin": Tin, "Tup": Tup, "tau": 5., "ta": 40. * 60. + vorl, "tb": 160. * 60. + vorl,
                 "tc": 370. * 60. + vorl})
        timeset35h = ParametersSet({"t_start": 0., "t_stop": 1000. * 60. + vorl, "delta_t": 5.0})
        Simulation35h = Simulate(MyHSM, timeset35h, Tset35h, "SimulationDoubleHeat3h30min" + FigureExtension)
        Simulation35h.TimeRun(ZoomInPanelA="Yes3")

        ### 5h ###
        Tset5h = ParametersSet(
                {"Ttype": 3, "Tin": Tin, "Tup": Tup, "tau": 5., "ta": 40. * 60. + vorl, "tb": 160. * 60. + vorl,
                 "tc": 460. * 60. + vorl})
        timeset5h = ParametersSet({"t_start": 0., "t_stop": 1000. * 60. + vorl, "delta_t": 5.0})
        Simulation5h = Simulate(MyHSM, timeset5h, Tset5h, "SimulationDoubleHeat5h" + FigureExtension)
        Simulation5h.TimeRun(ZoomInPanelA="Yes3")


    ###################################################################################
    #############  Simulated Schroda et al. 2000 exp with ARS enzyme  #################
    ###################################################################################

    if Double_Heat_ARS == "Yes":
        print("\nSTARTING TO Simulate DOUBLE HEAT SHOCK from Schroda et al. 2000 with ARS enzyme...\n")

        # Simulate mRNA and ARS after 1h HS (Fig 6b)
        TsetARS = ParametersSet(
                {"Ttype": 2, "Tin": 23., "Tup": 40., "tau": 5., "ta": 0. * 60. + vorl, "tb": 60. * 60. + vorl})
        timesetARS = ParametersSet({"t_start": 0., "t_stop": 360. * 60. + vorl, "delta_t": 5.0})
        SimulationARS = Simulate(MyHSM, timesetARS, TsetARS, "SimulationARSexperiment" + FigureExtension)
        SimulationARS.TimeRunPlusARS()

        # Simulate mRNA and ARS after double HS (Fig 7b)
        HSduration = 30.  # (min)
        TsetARSdoubleHS = ParametersSet(
                {"Ttype": 2, "Tin": 23., "Tup": 40., "tau": 5., "ta": 0. * 60. + vorl, "tb": HSduration * 60. + vorl})
        timesetARSdoubleHS = ParametersSet(
                {"t_start": 0., "t_stop": (2. * HSduration + 5 * 60. + 60.) * 60 + vorl, "delta_t": 5.0})
        SimulationARSdoubleHS = Simulate(MyHSM, timesetARSdoubleHS, TsetARSdoubleHS,
                                         "SimulationARSdoubleHSshort" + FigureExtension)
        EmptyListToBeFilled = []
        SimulationARSdoubleHS.TimeRunPlusARSdoubleHSMOD(EmptyListToBeFilled)


    #################################################################
    #############          24h behaviour          ###################
    #################################################################

    if Early_Late_24h_HS == "Yes":
        print("\nSTARTING TO Simulate Early HS, Late HS, 8h Recovery...\n")
        ### Early HS, Late HS, 8h Recovery (cf Hemme et al. 2014, fig. 8) ###
        TsetEarlyHSLateHSRecovery = ParametersSet(
                {"Ttype": 2, "Tin": 25., "Tup": 42., "tau": 5., "ta": 0. * 60. + vorl, "tb": 24. * 60. * 60. + vorl})
        timesetEarlyHSLateHSRecovery = ParametersSet(
                {"t_start": 0., "t_stop": (24. + 8.) * 60. * 60. + vorl, "delta_t": 5.0})
        SimulationEarlyHSLateHSRecovery = Simulate(MyHSM, timesetEarlyHSLateHSRecovery, TsetEarlyHSLateHSRecovery,
                                                   "SimulationEarlyHSLateHSRecovery" + FigureExtension)
        SimulationEarlyHSLateHSRecovery.TimeRun(ZoomInPanelA="Yes2")

    if Hot_Day_Behaviour == "Yes":
        print("\nSTARTING TO Simulate Daily (sinusoidal) T variation...\n")
        ### Daily (sinusoidal) T variation ###
        NperiodsBeforePlotting = 4
        timesetWEIRDTsin = ParametersSet({"t_start": 0., "t_stop": (NperiodsBeforePlotting+1)*24. * 60. * 60. + vorl, "delta_t": 5.0})
        TsetWEIRDTsin = ParametersSet({"Ttype": 7, "Tin": 22., "Tup": 40., "Period": 24 * 60 * 60})
        SimulationWEIRDTsin = Simulate(MyHSM, timesetWEIRDTsin, TsetWEIRDTsin, "SimulationWEIRDTsin" + FigureExtension)
        SimulationWEIRDTsin.TimeRun(ZoomInPanelA="Yes", tminMANUAL=("Yes", NperiodsBeforePlotting * 24. * 60. * 60.))


    #################################################################################################################
    ###############  Sistematic Study of HP production for different Temperatures/Durations of HS  ##################
    #################################################################################################################

    if HPproductionForVariousTempDurationsHS == "Yes":
        print("\nSTARTING TO Simulate Sistematic Study of HP production for different Temperatures/Durations of HS...\n")

        MyTEMPERATUREstart = 20.   # (deg C)
        MyTEMPERATUREstop = 45.    # (deg C)
        MyNstepsTEMPERATURE = 3 # put at least 30!!!  # (adimensional)
        MyDURATIONstart = 10*60.      # (seconds)
        MyDURATIONstop = 1000.*60. # (seconds)
        MyNstepsDURATION = 4 # put at least 30!!!     # (adimensional)

        StudyHPproductionForDifferentTemeperaturesDurationsHS(MyHSM,
                                                              MyTEMPERATUREstart,
                                                              MyTEMPERATUREstop,
                                                              MyNstepsTEMPERATURE,
                                                              MyDURATIONstart,
                                                              MyDURATIONstop,  
                                                              MyNstepsDURATION,
                                                              "ExploreTemperatureDurationHS" + FigureExtension)

    ####################################################################################
    #############  Study evolution of steady state (Oliver's suggestion)  ##############
    ####################################################################################

    if SteadyStateSystematicStudy == "Yes":
        print("\nSTARTING TO Simulate evolution of steady state...\n")

        InitialTemperature = 0  # (deg C)
        FinalTemperature = 100  # (deg C)
        #StudyUnfoldingRateFuncOfT(ThisParametrSet["kP0p"], InitialTemperature, FinalTemperature, "lower right", "UnfoldingRate" + FigureExtension) 
        # Units of measure are: (s^-1, deg C, deg C, none, none) 

        #MaximalNuPpOverPTested = 150.  # (s^-1)
        NumberOfPointsInNuPpRange = 51 # 21 # Use 51
        
        NumbEquations9or10or12 = 9
        IC_PplusPp = MyHSM.ParamsSetIC.CurrentParams["Pin"] + MyHSM.ParamsSetIC.CurrentParams["Phin"]
        IC_SplusSs = MyHSM.ParamsSetIC.CurrentParams["Sin"] + MyHSM.ParamsSetIC.CurrentParams["Ssin"]
        IC_GplusFsGplusFG = MyHSM.ParamsSetIC.CurrentParams["Gin"] + MyHSM.ParamsSetIC.CurrentParams["FsGin"] + MyHSM.ParamsSetIC.CurrentParams["FGin"] 
        print()
        print("LOOK HERE!!!")
        print()
        print(IC_PplusPp)
        print()
        print(IC_SplusSs)
        print()
        print(IC_GplusFsGplusFG)
        print()
        #StudySteadyStateAndEquilibriumEvolution(MyHSM, FigureExtension, ThisParametrSet, 
        #                                        DefaultParamSetForREACTIONS, MaximalNuPpOverPTested, NumberOfPointsInNuPpRange, 
        #                                        "EvolutionOfEquilibriumPoint9eqs" + FigureExtension, 
        #                                        "EvolutionOfEigenvaluesOfJacobianAtEquilibriumPoint9eqs" + FigureExtension,
        #                                        NumbEquations9or10or12,
        #                                        IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)

        StudySteadyStateAndEquilibriumEvolutionFuncOfTEMPERATURE(MyHSM, FigureExtension, ThisParametrSet, 
                                                                 DefaultParamSetForREACTIONS, InitialTemperature, FinalTemperature, NumberOfPointsInNuPpRange, 
                                                                 "EvolutionOfEquilibriumPoint9eqsFuncOfnuPp" + FigureExtension, 
                                                                 "EvolutionOfEquilibriumPoint9eqsFuncOfTEMPERATURE" + FigureExtension, 
                                                                 "EvolutionOfEigenvaluesOfJacobianAtEquilibriumPoint9eqsFuncOfnuPp" + FigureExtension,
                                                                 "EvolutionOfEigenvaluesOfJacobianAtEquilibriumPoint9eqsFuncOfTEMPERATURE" + FigureExtension,
                                                                 NumbEquations9or10or12,
                                                                 IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)
    
    ###########################################################################################################################
    #############  Run for vorl and see the steady state concentrations at t = vorl + 0 and tis stability  ###################
    ###########################################################################################################################

    if Study_Stability_Short == "Yes":
        print("\nSTARTING TO run for a very long time without HS...\n")

        FixedTemp = 20.  # (deg C)
        ShortTime = vorl

        ### 0) See which are the values at which the system settle on a short timescale ###
        TsetSteadyStateShort = ParametersSet({"Ttype": 0, "Tin": FixedTemp})
        timesetSteadyStateShort = ParametersSet({"t_start": 0., "t_stop": ShortTime, "delta_t": 5.0})
        SimulationSteadyStateShort = Simulate(MyHSM, timesetSteadyStateShort, TsetSteadyStateShort,
                                              "SimulationSteadyStateShort" + FigureExtension)
        SimulationSteadyStateShort.TimeRun()
        ### 1) Infer a guess for an equilibrium point by taking the variables' values at the final time of the simulation above ###
        yEquilibriumGuessShort = FinalValuesOf(SimulationSteadyStateShort)
        ### 2) Find a root of the system f(y), i.e. a point of equilibrium, starting from the guess for y previously found ###
        yEquilibriumPointShort = FindRootOfFuncAndPrint(ODEsSysthAsFunction, yEquilibriumGuessShort,
                                                        (ThisParametrSet, DefaultParamSetForREACTIONS, FixedTemp, ["No",0.],
                                                        12, 0, 0, 0))
        ### 3) Compute the Jacobian of f(y) and find its eigenvalues to determine the stability of the equilibrium point ###
        eigenvaluesShort = EigenvaluesOfJacobianAtEquilibrium(yEquilibriumPointShort, ThisParametrSet,
                                                              DefaultParamSetForREACTIONS, FixedTemp, ["No",0.],
                                                              12, 0, 0, 0)

    #####################################################################################################
    ###############  Simulate Unfolded Proteins as func Of time taken by Temp to go up  #################
    #####################################################################################################

    if UnfoldedPfuncOfTau == "Yes":
        print("\nSTARTING TO SIMULATE How Unfolded Proteins change as func Of time taken by Temp to go up...\n")

        Tup = 42.
        Tdown = 25.
        PlotTemperatureManyTau(Tup, Tdown, FigureExtension)

        TauMin = 10. # (s)
        TauMax = 24. * 60. * 60. # (s) = 24h
        NumberOfSteps = 16
        SimulationName = "SimulationUnfoldedPfuncOfTau" + FigureExtension
        FigureName = "UnfoldedPfuncOfTau" + FigureExtension

        ComputeMaxUnfoldedProteinsAsFunctionOfTimeToIncreaseTemperature(MyHSM, TauMin, TauMax, NumberOfSteps, SimulationName, FigureName, Tup, Tdown)


