
from copy import deepcopy
import math

#from HSM_ParametersClass import *
from HSM_SimulateClass import *



def GenerateRandomParametersSets(NumberOfRandomSets, FactorOf, DefaultParamSetRATES):

    """GENERATE SETS OF RANDOM PARAMETERS FROM A FLAT DISTRIBUTION CENTERED AROUND FIDUCIAL VALUES AND WITHIN A FACTOR OF FactorOf"""

    DictionaryOfMaxs = {}
    DictionaryOfMins = {}

    for key in DefaultParamSetRATES:
        #DictionaryOfMaxs[key] = deepcopy(DefaultParamSetRATES[key]*FactorOf)
        DictionaryOfMaxs[key] = deepcopy(DefaultParamSetRATES[key] - DefaultParamSetRATES[key]*FactorOf)
        #DictionaryOfMins[key] = deepcopy(DefaultParamSetRATES[key]/FactorOf)
        DictionaryOfMins[key] = deepcopy(DefaultParamSetRATES[key] + DefaultParamSetRATES[key]*FactorOf)

    import random

    ListOfManyDictionaryOfRandomParameters = []

    for i in range(NumberOfRandomSets):
        DictionaryOfRandomParameters = {}
        for key in DefaultParamSetRATES:
            RandomNumber = random.random()
            NewValue = deepcopy(DictionaryOfMins[key] + RandomNumber*(DictionaryOfMaxs[key]-DictionaryOfMins[key]))
            DictionaryOfRandomParameters[key] = NewValue
        ListOfManyDictionaryOfRandomParameters.append(deepcopy((DictionaryOfRandomParameters)))

    return ListOfManyDictionaryOfRandomParameters



def GenerateParametersSetsChangingOneParameter(NumberOfValuesForEachParameterk, FactorOfII, DefaultParamSetRATES):
    
    """GENERATE SETS OF PARAMETERS By Changing Only 1 PARAMETER AT A TIME"""

    ListOfManyDictionariesOfParametersVarying1by1 = []

    for key in DefaultParamSetRATES:
        DictionaryOfTestParameters = deepcopy(DefaultParamSetRATES)
        for j in range(NumberOfValuesForEachParameterk+1):
            NewValue = deepcopy(DefaultParamSetRATES[key] + FactorOfII * DefaultParamSetRATES[key] * ( 2*j - NumberOfValuesForEachParameterk)/NumberOfValuesForEachParameterk)
            DictionaryOfTestParameters[key] = deepcopy(NewValue)
            ListOfManyDictionariesOfParametersVarying1by1.append(deepcopy(DictionaryOfTestParameters))
    
    return ListOfManyDictionariesOfParametersVarying1by1



def ExtractDataControlFeedExpHSFandHSP90aFromFiles(DataFileNameHSFcontrol, DataFileNameHSP90Acontrol):

    ListToBeFilledWithResults = []    

    # CONTROL HSF 
    ColumnNumber = 2  # Only the time and the row data are used
    ListOfDataArraysHSF = []
    FromDataFileToArrays(DataFileNameHSFcontrol, ColumnNumber, ListOfDataArraysHSF)  # Read data file, put in list of arrays
    ListToBeFilledWithResults.append(ListOfDataArraysHSF[1])

    # Times
    ListToBeFilledWithResults.append(ListOfDataArraysHSF[0])

    # CONTROL HSP90a 
    ColumnNumber = 2  # Only the time and the row data are used
    ListOfDataArraysHSP90a = []
    FromDataFileToArrays(DataFileNameHSP90Acontrol, ColumnNumber, ListOfDataArraysHSP90a)  # Read data file, put in list of arrays
    ListToBeFilledWithResults.append(ListOfDataArraysHSP90a[1])

    return ListToBeFilledWithResults



def ComputePartOfRMSSimulationVsData(ListOfDataTimes, timesetDataRMS, SimulationExperimentDataRMS, ListOfDataHSF, ListOfDataHSP90a):
    """Compute the Sum over all datapoints of ( Xth - Xdata )^2, for 1 feeding experiment, taking into caccount HSF + HSP """

    ##### 2-C1.1: EXTRACT SIMULATION VALUES AT PROPER TIME POINTS FOR COMPARISON WITH DATA
    ListOfmRNAHSFsimulation = []
    ListOfmRNAHSPsimulation = []
    ListOfTimesForDatapoints = ListOfDataTimes #= [0., 15., 30., 45., 60., 120.]

    for val in ListOfTimesForDatapoints:
        j = ( val * 60. + vorl ) / timesetDataRMS.CurrentParams["delta_t"] # ( seconds/seconds = adimensional )
        ListOfmRNAHSFsimulation.append(SimulationExperimentDataRMS.RF[j])
        ListOfmRNAHSPsimulation.append(SimulationExperimentDataRMS.RHP[j])
    ArrayOfmRNAHSFsimulation = np.array([val for sublist in ListOfmRNAHSFsimulation for val in sublist])
    ArrayOfmRNAHSPsimulation = np.array([val for sublist in ListOfmRNAHSPsimulation for val in sublist])

    ##### 2-C1.2: COMPUTE pieces of LS TH VS DATA - STAUROSPORINE

    # print("We now wants to compare these...")
    # print(ListOfDataHSF_stau)
    # print(ListOfDataHSP90a_stau)
    # print(ArrayOfmRNAHSFsimulation)
    # print(ArrayOfmRNAHSPsimulation)

    k = 0
    SumOverDataPointsHSF = 0.
    for val in ListOfDataHSF:
         DeviationHSF = ArrayOfmRNAHSFsimulation[k]/max(ArrayOfmRNAHSFsimulation) - ListOfDataHSF[k]/max(ListOfDataHSF)
         SumOverDataPointsHSF = SumOverDataPointsHSF + pow(DeviationHSF, 2)
         k = k + 1

    l = 0
    SumOverDataPointsHSP90a = 0.

    for val in ListOfDataHSP90a:
        DeviationHSP90a = ArrayOfmRNAHSPsimulation[l]/max(ArrayOfmRNAHSPsimulation) - ListOfDataHSP90a[l]/max(ListOfDataHSP90a)
        SumOverDataPointsHSP90a = SumOverDataPointsHSP90a + pow(DeviationHSP90a, 2)
        l = l + 1

    SumOverDatapoints = SumOverDataPointsHSF + SumOverDataPointsHSP90a

    return SumOverDatapoints



def PlotSimulationVsDataFeeding(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, timesetDataRMS, ListOfDataHSF, ListOfDataHSP90a, ListOfDataTimes, FigureName, FigureExtension, FolderContainingDataVsSimuCalibration):
    """Plot mRNAs for HSF and HSP90a feeding experiments data VS simulations to see if it makes sense"""

    fig = figure()

    ############ Simulations  

    ax1 = plt.subplot(121)

    SubPlot(ax1, SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, 'Time (min)', 0.,
            (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., " ", 0, 0, "upper right", "A",
            Legendfontsize="small", Legendfancybox=True)

    ax2 = plt.subplot(122)

    SubPlot(ax2, SimulationFeedingControlsDataRMStimes, ListForPlottingHSP, 'Time (min)', 0.,
            (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., " ", 0, 0, "upper right", "B",
            Legendfontsize="small", Legendfancybox=True)

    ############ and Data Points   

    ListOfDataHSFNORM = []
    k = 0
    for val in ListOfDataHSF:
         ListOfDataHSFNORM.append(ListOfDataHSF[k]/max(ListOfDataHSF))
         k = k + 1 

    ListOfDataHSP90aNORM = []
    k = 0
    for val in ListOfDataHSP90a:
         ListOfDataHSP90aNORM.append(ListOfDataHSP90a[k]/max(ListOfDataHSP90a))
         k = k + 1 

    DataLegend = [r"Data" + str(FigureName)]
    DataSubPlot(ax1, ListOfDataTimes, [ListOfDataTimes, np.asarray(ListOfDataHSFNORM)], 'Time (min)', 0.,
                (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., r"mRNA$_{HSF}$ (Normalizet to Max)", 0., 1.,  "upper right", DataLegend, "",
                Legendfontsize="small", Legendfancybox=True, Black = "Yes")

    DataLegend = [r"Data" + str(FigureName)]
    DataSubPlot(ax2, ListOfDataTimes, [ListOfDataTimes, np.asarray(ListOfDataHSP90aNORM)], 'Time (min)', 0.,
                (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., r"mRNA$_{HSP}$ (Normalizet to Max)", 0., 1.,  "upper right", DataLegend, "",
                Legendfontsize="small", Legendfancybox=True, Black = "Yes")

    PlotAndSave(fig, FolderContainingDataVsSimuCalibration + "FittingToData" + str(FigureName) + FigureExtension, "PS", 1, 0)


##################################################### LOOK HERE!!!!!!!!!!!!!!!!!!!!! ##########################################################
def PlotSimulationVsDataFeedingModelVSFittedData(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, timesetDataRMS, ListOfDataHSF, ListOfDataHSP90a, ListOfDataTimes, FigureName, FigureExtension, FolderContainingDataVsSimuCalibration):
    """Plot mRNAs for HSF and HSP90a feeding experiments data VS simulation best fit for paper"""

    ##### for key in AllDataControlsFeeding[5]:#ListOfFeedingKeys:

    fig = figure()

    ############ Simulations  

    ax1 = plt.subplot(121)

    SubPlot(ax1, SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, 'Time (min)', 0.,
            (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., " ", 0, 0, "upper right", "A",
            Legendfontsize="small", Legendfancybox=True, Black = "Yes")

    ax2 = plt.subplot(122)

    SubPlot(ax2, SimulationFeedingControlsDataRMStimes, ListForPlottingHSP, 'Time (min)', 0.,
            (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., " ", 0, 0, "upper right", "B",
            Legendfontsize="small", Legendfancybox=True, Black = "Yes")
    
    ############ and Data Points   

    ListOfFeedingKeys= ["stau", "radi", "ChloCyc", "canav", "Gelda", "CaChel"]
    DictionaryForLegend = {"stau": "Staurosporine", 
                           "radi": "Radicicol", 
                           "ChloCyc": "Chlor. / Cyclo.", 
                           "canav": "Canavanine", 
                           "Gelda": "Geldanamicil", 
                           "CaChel": "Calcium Chelator"}

    i = 0
    for key in ListOfFeedingKeys:

        ListOfDataHSFNORM = []
        k = 0
        for val in ListOfDataHSF[key]:
             ListOfDataHSFNORM.append(ListOfDataHSF[key][k]/max(ListOfDataHSF[key]))
             k = k + 1 

        ListOfDataHSP90aNORM = []
        k = 0
        for val in ListOfDataHSP90a[key]:
            ListOfDataHSP90aNORM.append(ListOfDataHSP90a[key][k]/max(ListOfDataHSP90a[key]))
            k = k + 1 

        DataLegend = [r"Data Control " + DictionaryForLegend[key] + " Exp."]
        DataSubPlotMOD(ax1, ListOfDataTimes[key], [ListOfDataTimes[key], np.asarray(ListOfDataHSFNORM)], 'Time (min)', 0.,
                    (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., r"Concentration of mRNA$_{HSF}$ (normalized to max)", 0., 1.,  "upper right", DataLegend, "",
                    Legendfontsize="small", Legendfancybox=True, ColorNumber = i)

        DataLegend = [r"Data Control " + DictionaryForLegend[key] + " Exp."]
        DataSubPlotMOD(ax2, ListOfDataTimes[key], [ListOfDataTimes[key], np.asarray(ListOfDataHSP90aNORM)], 'Time (min)', 0.,
                    (timesetDataRMS.CurrentParams["t_stop"] - vorl) / 60., r"Concentration of mRNA$_{HSP}$ (normalized to max)", 0., 1.,  "upper right", DataLegend, "",
                    Legendfontsize="small", Legendfancybox=True, ColorNumber = i)
    
        i = i+1

    PlotAndSave(fig, FolderContainingDataVsSimuCalibration + "FittingToDataPAPERversion" + str(FigureName) + FigureExtension, "PS", 1, 0)



def ExtractDataControlsFeedingExperimentsFromFilesIntoListOfDictionaries():
    """ EXTRACT EXPERIMENTAL DATA (FEEDING EXPERIMENTS) FROM FILES ALEXANDER SKUPIN """

    ListOfFeedingKeys= ["stau", "radi", "ChloCyc", "canav", "Gelda", "CaChel"]

    DictionaryOfHSFcontrolFiles, DictionaryOfHSP90acontrolFiles = {}, {}
    for key in ListOfFeedingKeys:
        DictionaryOfHSFcontrolFiles.update({key : "DataFilesOriginals/" + "hsfcontrol_" + str(key) + ".csv"}) 
        DictionaryOfHSP90acontrolFiles.update({key : "DataFilesOriginals/" + "hsp90acontrol_" + str(key) + ".csv"})

    DictionaryOfListsOfDataHSF, DictionaryOfListsOfDataTimes, DictionaryOfListsOfDataHSP90a = {}, {}, {}
    for key in ListOfFeedingKeys:
        ListOfListsOfExtractedData = ExtractDataControlFeedExpHSFandHSP90aFromFiles(DictionaryOfHSFcontrolFiles[key], DictionaryOfHSP90acontrolFiles[key])
        DictionaryOfListsOfDataHSF.update({key : ListOfListsOfExtractedData[0]})
        DictionaryOfListsOfDataTimes.update({key : ListOfListsOfExtractedData[1]})
        DictionaryOfListsOfDataHSP90a.update({key : ListOfListsOfExtractedData[2]})

    # Create Temperature settings and Time settings reproducing the experimental setup of controls for feeding experiments
    TsetSchmol2013dataRMS = ParametersSet({"Ttype": 1, "Tin": 25., "Tup": 40., "tau": 5., "ta": 0. * 60. + vorl})
    timeset240minsDataRMS = ParametersSet({"t_start": 0., "t_stop": 240. * 60. + vorl, "delta_t": 5.0})

    AllDataControlsFeeding = (DictionaryOfListsOfDataTimes, DictionaryOfListsOfDataHSF, DictionaryOfListsOfDataHSP90a, TsetSchmol2013dataRMS, timeset240minsDataRMS, ListOfFeedingKeys)

    return AllDataControlsFeeding



def ComputeRMSfeedingForGivenParameterSet(ParamSetRates, ParamSetForREACTIONS, ParamSetInitCond, OutputHSFandHSPtoPlot, AllDataControlsFeeding):
    """ Function to compute RMS w.r.t. data of the controls of the feeding experiments from Schmollinger et al. 2013"""
    # OutputHSFandHSPtoPlot = "Yes" or "No", for creating also output for time course plots or not, respectively.

    ############ 1] ############ NEED TO GENERATE A MODEL FOR EVERY PARAMS SET

    # Temeprature and Time parameters sets reproducing the experiments
    TsetSchmol2013dataRMS = AllDataControlsFeeding[3]
    timeset240minsDataRMS = AllDataControlsFeeding[4]

    # Set default parameter values
    TestParamSetIC = ParametersSet(ParamSetInitCond)
    TestParamSetRATES = ParametersSet(ParamSetRates)

    ParamSetForREACTIONS["piRFconst"] = ParamSetRates["ketaRF"]*ParamSetRates["ketaF"]/ParamSetRates["kFpi0"]*0.17/8.
    ParamSetForREACTIONS["piRHPconst"] = ParamSetRates["ketaRHP"]*ParamSetRates["ketaHP"]/ParamSetRates["kpiHP"]*17.5
    ParamSetForREACTIONSobject = ParametersSet(ParamSetForREACTIONS)

    # Create an object of the class "Heat shock models" with these parameters
    TestHSM = HeatShockModel(TestParamSetIC, TestParamSetRATES, ParamSetForREACTIONSobject)

    ############ 2] ############ NEXT, FOR EACH MODEL COMPUTE pieces necessary to compute the RMS W.R.T. DATA 

    ##### 2-A: SIMULATE CONTROLS OF ALL FEEDING EXPERIMENTS (one simulation only for all the datasets of feeding experiments!!!)
    SimulationFeedingControlsDataRMS = Simulate(TestHSM, timeset240minsDataRMS, TsetSchmol2013dataRMS, "xyzUSELESSxyz") # "testFeedingControls" + str(i) + FigureExtension
    SimulationFeedingControlsDataRMS.TimeRun(AvoidPlots="Yes")

    ##### 2-B : Extract from input the dictionaries containing the data
    DictionaryOfListsOfDataTimes = AllDataControlsFeeding[0]
    DictionaryOfListsOfDataHSF = AllDataControlsFeeding[1]
    DictionaryOfListsOfDataHSP90a = AllDataControlsFeeding[2]
    ListOfFeedingKeys = AllDataControlsFeeding[5]

    ##### 2-C: Compute Part Of RMS Simulation Vs Data FOR EACH DIFFERENT FEEDING EXPERIMENT
    SumOverDatapointsFeeding = {}
    for key in ListOfFeedingKeys:
        SumOverDatapointsFeeding.update( {key : ComputePartOfRMSSimulationVsData(DictionaryOfListsOfDataTimes[key], timeset240minsDataRMS, SimulationFeedingControlsDataRMS, DictionaryOfListsOfDataHSF[key], DictionaryOfListsOfDataHSP90a[key])} )

    ############ 3] ############ Put together the pieces for each different dataset into 1 single RMS value!!!
    NumberOfDataPoints, SumOverDatapoints = 0., 0.
    for key in ListOfFeedingKeys:
        NumberOfDataPoints = NumberOfDataPoints + len(DictionaryOfListsOfDataHSF[key]) + len(DictionaryOfListsOfDataHSP90a[key])
        SumOverDatapoints = SumOverDatapoints + SumOverDatapointsFeeding[key]

    RootMeanSquareDeviation = math.sqrt( SumOverDatapoints / NumberOfDataPoints )

    print("\n" + str(RootMeanSquareDeviation) + "\n")

    if OutputHSFandHSPtoPlot == "No":
        output = (RootMeanSquareDeviation)
    elif OutputHSFandHSPtoPlot == "Yes":
        output = (RootMeanSquareDeviation, SimulationFeedingControlsDataRMS.RF, SimulationFeedingControlsDataRMS.RHP, SimulationFeedingControlsDataRMS.t)
    else:
        print("\nError in RMS feeding function!!!\n")

    return output



def Convert_ORIGINAL_to_RESCALED_ParameterSet(ORIGINAL_ParameterSetDictionary, RescalingFactorsDictionary):

    RESCALED_ParameterSetDictionary = {}
    for key in RescalingFactorsDictionary:
        RESCALED_ParameterSetDictionary.update({ key : deepcopy(ORIGINAL_ParameterSetDictionary[key])/deepcopy(RescalingFactorsDictionary[key]) })
    return deepcopy(RESCALED_ParameterSetDictionary)



def Convert_RESCALED_to_ORIGINAL_ParameterSet(RESCALED_ParameterSetDictionary, RescalingFactorsDictionary):

    ORIGINAL_ParameterSetDictionary = {}
    for key in RescalingFactorsDictionary:
        ORIGINAL_ParameterSetDictionary.update({ key : deepcopy(RESCALED_ParameterSetDictionary[key]) * deepcopy(RescalingFactorsDictionary[key]) })
    return deepcopy(ORIGINAL_ParameterSetDictionary)



def ComputeRMSfeedingForGivenParameterSet_RESCALED_PARAMETERS(ParamSetRates_RESCALED, ParamSetForREACTIONS, ParamSetInitCond, OutputHSFandHSPtoPlot, AllDataControlsFeeding, RescalingFactorsDictionary):
    """ Function that does exactly what ComputeRMSfeedingForGivenParameterSet does, but parameters are rescaled to their fiducial value."""
    # This is a PRECONDITIONING, it serves to have a function which is easier to treat numerically with optimization algorithms as the gradient search
    ParamSetRates_ORIGINAL = deepcopy(Convert_RESCALED_to_ORIGINAL_ParameterSet(ParamSetRates_RESCALED, RescalingFactorsDictionary))
    Output = ComputeRMSfeedingForGivenParameterSet(ParamSetRates_ORIGINAL, ParamSetForREACTIONS, ParamSetInitCond, OutputHSFandHSPtoPlot, AllDataControlsFeeding)
    return Output



def FindMinimumOfFunctionUsingGoldenRatioBisectionMethod(FunctionToMinimize, LowerBound, UpperBound, Tolerance):

    GoldenRatio = 2./(math.sqrt(5.) + 1)

    ### Use the golden ratio to set the initial test points
    x1 = UpperBound - GoldenRatio*(UpperBound - LowerBound)
    x2 = LowerBound + GoldenRatio*(UpperBound - LowerBound)

    ### Evaluate the function at the test points
    f1 = FunctionToMinimize(x1)
    f2 = FunctionToMinimize(x2)

    i = 0

    while ( (abs(UpperBound - LowerBound) > Tolerance) and i <= 15):

        i = i + 1

        if f2 > f1:
        # then the minimum is to the left of x2
        # let x2 be the new upper bound
        # let x1 be the new upper test point

            ### Set the new upper bound
            UpperBound = deepcopy(x2)
            ### Set the new upper test point
            ### Use the special result of the golden ratio
            x2 = deepcopy(x1)
            f2 = deepcopy(f1)

            ### Set the new lower test point
            x1 = UpperBound - GoldenRatio*(UpperBound - LowerBound)
            f1 = FunctionToMinimize(x1) 

        elif f2 < f1: 
            # the minimum is to the right of x1
            # let x1 be the new lower bound
            # let x2 be the new lower test point

            ### Set the new lower bound
            LowerBound = deepcopy(x1)

            ### Set the new lower test point
            x1 = deepcopy(x2)
            f1 = deepcopy(f2)

            ### Set the new upper test point
            x2 = LowerBound + GoldenRatio*(UpperBound - LowerBound)
            f2 = FunctionToMinimize(x2)

        else:
            print("Error in Golden Rule minimization algorithm!")

        print(str(i) + "     " + str(x1) + "  " + str(x2) + "     " + str(f1) + "  " + str(f2))

    ### Use the mid-point of the final interval as the estimate of the optimzer
    EstimatedMinimizer = (LowerBound + UpperBound)/2.
    return EstimatedMinimizer



def NormalizeRnaCurvesFromSimulationsToMaxForPlot(RMSvalue, mRNA_HSF_simulation, mRNA_HSP_simulation, ListForPlottingHSF, ListForPlottingHSP, CurveNameInLegend, timesetDataRMS):

    IndexOfT0seconds = int(vorl/timesetDataRMS.CurrentParams["delta_t"]) # ( seconds/seconds = adimensional )
    RangeOfIndexesForPositiveTimes = range(IndexOfT0seconds,len(mRNA_HSF_simulation),1)

    Max_HSF_FeedingControls_simulation = max(np.max(mRNA_HSF_simulation[kkk]) for kkk in RangeOfIndexesForPositiveTimes)
    Max_HSP_FeedingControls_simulation = max(np.max(mRNA_HSP_simulation[kkk]) for kkk in RangeOfIndexesForPositiveTimes)

    Y_HSF_FeedingControls_simulationNORM = (np.asarray(mRNA_HSF_simulation) / Max_HSF_FeedingControls_simulation) # * 100. for %
    Y_HSP_FeedingControls_simulationNORM = (np.asarray(mRNA_HSP_simulation) / Max_HSP_FeedingControls_simulation) # * 100. for %

    ListForPlottingHSF.append([CurveNameInLegend, Y_HSF_FeedingControls_simulationNORM])
    ListForPlottingHSP.append([CurveNameInLegend, Y_HSP_FeedingControls_simulationNORM])


##############################################################################################
##################################### RMS for DOUBLE HS ######################################
##############################################################################################

def ExtractDataDoubleHSFromFiles(DataFileName, ListOfDoubleHSKeys):

    ColumnNumber = 6
    ListOfDataArray = []
    FromDataFileToArrays(DataFileName, ColumnNumber, ListOfDataArray)  # Read data file, put in list of arrays

    TimePointsDoubleHS = np.array(ListOfDataArray[0])

    DictionaryOfResults2HS = {}
    index = 1
    for key in ListOfDoubleHSKeys:
        DictionaryOfResults2HS.update( {key : ListOfDataArray[index]} )
        index = index + 1

    ListToBeFilledWithResults = []
    ListToBeFilledWithResults.append(TimePointsDoubleHS)
    ListToBeFilledWithResults.append(DictionaryOfResults2HS)

    return ListToBeFilledWithResults

#############################

def ExtractDataControlsDoubleHSExperimentFromFilesIntoListOfDictionaries():
    """ EXTRACT EXPERIMENTAL DATA (DOUBLE HEAT SHOCK EXPERIMENTS) FROM FILES """

    ListOfDoubleHSKeys= ["singleHS", "2hdoubleHS", "3hdoubleHS", "4hdoubleHS", "5hdoubleHS"]
    PathOfFileContainingAll2HSData = "DataFiles/DataShroda2000ARSFig7b.dat" 
    ListToBeFilledWithResults = ExtractDataDoubleHSFromFiles(PathOfFileContainingAll2HSData, ListOfDoubleHSKeys)

    #ArrayTimeDataPointsDoubleHS = ListToBeFilledWithResults[0] # (in seconds)
    DictionaryOfArraysOfData2HS = ListToBeFilledWithResults[1]

    DictionaryTimeDataPointsDoubleHS = {
    "singleHS" : np.array([0.00, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0]), 
    "2hdoubleHS" : np.array([0.00, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0]), 
    "3hdoubleHS" : np.array([0.00, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]), 
    "4hdoubleHS" : np.array([0.00, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0, 390.0]), 
    "5hdoubleHS" : np.array([0.00, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0, 390.0, 420.0]),
    }

    # Create Temperature settings and Time settings reproducing the starting experimental setup (they will be modified when solving the ODE for different 2HSs)
    HSduration = 30.  # (min)
    TsetDoubleHSdataRMS = ParametersSet({"Ttype": 2, "Tin": 23., "Tup": 40., "tau": 5., "ta": 0. * 60. + vorl, "tb": HSduration * 60. + vorl})
    timesetDoubleHSDataRMS = ParametersSet({"t_start": 0., "t_stop": (2. * HSduration + 5 * 60. + 60.) * 60 + vorl, "delta_t": 5.0})

    USELESS = "This should not appear. If you see it anywere, it means something is wrong. It is needed to keep the number of elements of AllDataControlsDoubleHS"

    AllDataControlsDoubleHS = (DictionaryTimeDataPointsDoubleHS, DictionaryOfArraysOfData2HS, USELESS, TsetDoubleHSdataRMS, timesetDoubleHSDataRMS, ListOfDoubleHSKeys)

    return AllDataControlsDoubleHS

#############################

def ComputePartOfRMSSimulationVsDataDoubleHS(ListOfDataTimes, timesetDoubleHSDataRMS, TimePointsForAllSimulations, ARSconcentrationForOneSetup, ListOfDataARSactivity, AbsoluteMaxData, AbsoluteMaxSimulation):
    """ Compute the Sum over all datapoints of ( Xth - Xdata )^2, for 1 feeding experiment, taking into caccount HSF + HSP """

    # OutputOfDoubleHSsimulation = [self.t, [ [Legend_singleHS, Yval_singleHS=[]], ...,  [Legend_2HS5h,Yval_2HS5h=[]] ] ]

    ##### 2-C1.1: EXTRACT SIMULATION VALUES AT PROPER TIME POINTS FOR COMPARISON WITH DATA
    ListOfARSsimulation = []
    ListOfTimesForDatapoints = deepcopy(ListOfDataTimes) #= [0., 30., 60., 90., 120., etc.]

    for val in ListOfTimesForDatapoints:
        j = ( val * 60. + vorl ) / timesetDoubleHSDataRMS.CurrentParams["delta_t"] # ( seconds/seconds = adimensional )
        ListOfARSsimulation.append(ARSconcentrationForOneSetup[j])
    ArrayOfARSsimulation = np.array([val for sublist in ListOfARSsimulation for val in sublist])

    #print("We now want to compare these...")
    #print(ListOfDataARSactivity)
    #print()
    #print(ArrayOfARSsimulation)
    #print()

    ##### 2-C1.2: COMPUTE pieces of LS TH VS DATA - STAUROSPORINE

    k = 0
    SumOverDataPointsARS = 0.
    for val in ArrayOfARSsimulation:
         DeviationHSF = ArrayOfARSsimulation[k]/AbsoluteMaxSimulation - ListOfDataARSactivity[k]/AbsoluteMaxData
         SumOverDataPointsARS = SumOverDataPointsARS + pow(DeviationHSF, 2)
         #print(str(k) + "   " + str(SumOverDataPointsARS))
         k = k + 1

    return SumOverDataPointsARS



def ComputeRMSdoubleHSforGivenParameterSet(ParamSetRates, ParamSetForREACTIONS, ParamSetInitCond, OutputARSactiviryToPlot, AllDataControlsDoubleHS):
    """ Function to compute RMS w.r.t. data of the controls of the feeding experiments from Schmollinger et al. 2013"""
    # OutputARSactiviryToPlot = "Yes" or "No", for creating also output for time course plots or not, respectively.

    ############ 1] ############ NEED TO GENERATE A MODEL FOR EVERY PARAMS SET

    # Temeprature and Time parameters sets reproducing the experiments
    TsetDoubleHSdataRMS = AllDataControlsDoubleHS[3]
    timesetDoubleHSDataRMS = AllDataControlsDoubleHS[4]

    # Set default parameter values
    TestParamSetIC = ParametersSet(ParamSetInitCond)
    TestParamSetRATES = ParametersSet(ParamSetRates)

    ParamSetForREACTIONS["piRFconst"] = ParamSetRates["ketaRF"]*ParamSetRates["ketaF"]/ParamSetRates["kFpi0"]*0.17/8.
    ParamSetForREACTIONS["piRHPconst"] = ParamSetRates["ketaRHP"]*ParamSetRates["ketaHP"]/ParamSetRates["kpiHP"]*17.5
    ParamSetForREACTIONSobject = ParametersSet(ParamSetForREACTIONS)

    # Create an object of the class "Heat shock models" with these parameters
    TestHSM = HeatShockModel(TestParamSetIC, TestParamSetRATES, ParamSetForREACTIONSobject)

    ############ 2] ############ NEXT, COMPUTE pieces necessary to compute the RMS W.R.T. DATA 

    ##### 2-A: SIMULATE the single HS + the 4 2HS
    SimulationARSdoubleHSdataRMS = Simulate(TestHSM, timesetDoubleHSDataRMS, TsetDoubleHSdataRMS, "Useless")
    EmptyListToExtractOutput = []
    SimulationARSdoubleHSdataRMS.TimeRunPlusARSdoubleHS(EmptyListToExtractOutput, AvoidPlots="Yes")

    ##### 2-B : Extract from input the dictionaries containing the data
    DictionaryOfListsOfDataTimes = AllDataControlsDoubleHS[0]
    DictionaryOfListsOfDataARSactivity = AllDataControlsDoubleHS[1]
    ListOfDoubleHSKeys = AllDataControlsDoubleHS[5]

    OutputOfDoubleHSsimulation = EmptyListToExtractOutput[0]
    #print("I am the one " + str(OutputOfDoubleHSsimulation))

    TimePointsForAllSimulations = OutputOfDoubleHSsimulation[0]

    ARSconcentrationForEachHSsetupDictionary = {
    "singleHS" : OutputOfDoubleHSsimulation[1][0][1],
    "2hdoubleHS" : OutputOfDoubleHSsimulation[1][1][1],
    "3hdoubleHS" : OutputOfDoubleHSsimulation[1][2][1],
    "4hdoubleHS" : OutputOfDoubleHSsimulation[1][3][1],
    "5hdoubleHS" : OutputOfDoubleHSsimulation[1][4][1]
    }

    EmptyARSdataMaximaList = []
    for key in ListOfDoubleHSKeys:
        massimo = deepcopy(max(DictionaryOfListsOfDataARSactivity[key]))
        EmptyARSdataMaximaList.append(massimo)
    AbsoluteMaxData = max(EmptyARSdataMaximaList)

    EmptyARSSimulationMaximaList = []
    for key in ListOfDoubleHSKeys:
        massimo = deepcopy(max(ARSconcentrationForEachHSsetupDictionary[key]))
        EmptyARSSimulationMaximaList.append(massimo)
    AbsoluteMaxSimulationList = max(EmptyARSSimulationMaximaList)
    AbsoluteMaxSimulation = AbsoluteMaxSimulationList[0]

    #print()
    #print("MAX")
    #print()
    #print(AbsoluteMaxData)
    #print()
    #print(AbsoluteMaxSimulation)
    #print()

    ##### 2-C: Compute Part Of RMS Simulation Vs Data FOR EACH DIFFERENT FEEDING EXPERIMENT
    SumOverDatapointsDoubleHS = {}
    for key in ListOfDoubleHSKeys:
        #print("I am into this loop!!!")
        SumOverDatapointsDoubleHS.update( {key : ComputePartOfRMSSimulationVsDataDoubleHS(DictionaryOfListsOfDataTimes[key], timesetDoubleHSDataRMS, TimePointsForAllSimulations, ARSconcentrationForEachHSsetupDictionary[key], DictionaryOfListsOfDataARSactivity[key], AbsoluteMaxData, AbsoluteMaxSimulation)} )

    ############ 3] ############ Put together the pieces for each different dataset into 1 single RMS value!!!
    NumberOfDataPoints, SumOverDatapoints = 0., 0.
    for key in ListOfDoubleHSKeys:
        NumberOfDataPoints = NumberOfDataPoints + len(DictionaryOfListsOfDataTimes[key])
        SumOverDatapoints = SumOverDatapoints + SumOverDatapointsDoubleHS[key]
        #print(key)
        #print(NumberOfDataPoints)
        #print(SumOverDatapoints)

    RootMeanSquareDeviationDoubleHS = math.sqrt( SumOverDatapoints / NumberOfDataPoints )

    print("\n" + str(RootMeanSquareDeviationDoubleHS) + "\n")

    #if OutputHSFandHSPtoPlot == "No":
    output = (RootMeanSquareDeviationDoubleHS)
    #elif OutputHSFandHSPtoPlot == "Yes":
    #    output = (RootMeanSquareDeviation, SimulationFeedingControlsDataRMS.RF, SimulationFeedingControlsDataRMS.RHP, SimulationFeedingControlsDataRMS.t)
    #else:
    #    print("\nError in RMS feeding function!!!\n")

    return output



def ComputeRMStotalForGivenParameterSet_RESCALED_PARAMETERS(ParamSetRates_RESCALED, ParamSetForREACTIONS, ParamSetInitCond, AllDataControlsFeeding, RescalingFactorsDictionary, AllDataControlsDoubleHS):
    """ Function that does exactly what ComputeRMSfeedingForGivenParameterSet does, but parameters are rescaled to their fiducial value."""

    # This is a PRECONDITIONING, it serves to have a function which is easier to treat numerically with optimization algorithms as the gradient search
    ParamSetRates_ORIGINAL = deepcopy(Convert_RESCALED_to_ORIGINAL_ParameterSet(ParamSetRates_RESCALED, RescalingFactorsDictionary))

    RMS_Feeding=ComputeRMSfeedingForGivenParameterSet(ParamSetRates_ORIGINAL,ParamSetForREACTIONS,ParamSetInitCond,"No",AllDataControlsFeeding)
    RMS_DoubleHS=ComputeRMSdoubleHSforGivenParameterSet(ParamSetRates_ORIGINAL,ParamSetForREACTIONS,ParamSetInitCond,"No",AllDataControlsDoubleHS)
    Output = deepcopy(RMS_Feeding + RMS_DoubleHS)

    return Output






