


from HSM_SimulateClass import *
from HSM_StudyEquilibrium import *
from HSM_StudyHPproduction import *
from HSM_VaryParamsRMSvsData import *



def GenerateMCRandomOrNotParSetsAndComputeRMSFeeding(SwitchRandomSetsOrParametersK1by1Sets, MyNumberOfRandomSets, MyNumberOfValuesForEachParameterk, FactorOfRandom, FactorOfK1by1, FolderContainingCsvFiles, FolderContainingDataVsSimuCalibration, NameOfOutputFileRMSmanyParamsSets, NameOfOutputFileKeys, StartingParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, FigureExtension):

    """BIG function 1 of calibration: Generates a random (MC) params sets or change params 1by1, and computes RMS Feeding"""

    ############ 1) ############ GENERATE SETS OF PARAMETERS FOR RATES

    ### 1-a) GENERATE SETS OF RANDOM PARAMETERS FROM A FLAT DISTRIBUTION CENTERED AROUND FIDUCIAL VALUES AND WITHIN A FACTOR OF 2
    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
        NumberOfRandomSets = MyNumberOfRandomSets # 3 # we used up to 100000 
        FactorOf = FactorOfRandom # 0.5 # 0.5 means 50% variation of the parameter
        ListOfManyDictionariesOfParameters = GenerateRandomParametersSets(NumberOfRandomSets, FactorOf, StartingParamSetRATES)
    ### 1-b) GENERATE SETS OF PARAMETERS By Changing Only 1 PARAMETER AT A TIME
    elif SwitchRandomSetsOrParametersK1by1Sets == "ParametersK1by1":
        NumberOfValuesForEachParameterk = MyNumberOfValuesForEachParameterk # 1 # we used 100 
        FactorOf = FactorOfK1by1 # 0.5 # 0.5 means 50% variation of the parameter
        ListOfManyDictionariesOfParameters = GenerateParametersSetsChangingOneParameter(NumberOfValuesForEachParameterk, FactorOf, StartingParamSetRATES)#, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding)
    else:
        print("Error in the value of the switch SwitchRandomSetsOrParametersK1by1Sets!!!")

    ############ 2) ############ NOW NEED TO LOOP OVER EVERY PARAMETER SET

    ### Prepare list to be filled with RMS values
    ListOfRootMeanSquares = []

    ### Prepare lists necessary for plotting the HSP(t) and HSF(t) curves afterward
    ListForPlottingHSF, ListForPlottingHSP = [], []

    i = 0
    for ParamSetRates in ListOfManyDictionariesOfParameters:

        ##### A: Call the function which computes RMS w.r.t. feeding controls data (and rna(t) for HSP and HSF for plotting)
        ResultOfComputingRMSfromFeeding = ComputeRMSfeedingForGivenParameterSet(ParamSetRates, TestParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsFeeding)

        ##### B: Split its output
        RMS_simulation = ResultOfComputingRMSfromFeeding[0]
        mRNA_HSF_simulation = ResultOfComputingRMSfromFeeding[1]
        mRNA_HSP_simulation = ResultOfComputingRMSfromFeeding[2]
        SimulationFeedingControlsDataRMStimes = ResultOfComputingRMSfromFeeding[3]

        ##### C: Append RMS value to the list of all RMS
        ListOfRootMeanSquares.append(RMS_simulation)

        ##### D: EXTRACT AND NORMALIZE EACH SIMULATION CURVE TO ITS MAXIMUM (FOR THE PLOTTING ONLY!)

        CurveNameInLegend = "Param Set " + str(i) + ", RMS = " + str(round(ListOfRootMeanSquares[i], 4))
        if len(ListForPlottingHSF) <= 10:
            NormalizeRnaCurvesFromSimulationsToMaxForPlot(ListOfRootMeanSquares[i], mRNA_HSF_simulation, mRNA_HSP_simulation, 
                                                          ListForPlottingHSF, ListForPlottingHSP, CurveNameInLegend, AllDataControlsFeeding[4])#timeset240minsDataRMS)
        ##### E: ITERATE
        i = i + 1

    ############ 3) ############ ...AND FIND THE MINIMUM OF THE RMS (and the associated parameter set)
    
    MinimumRMSdistance = min(ListOfRootMeanSquares)
    IndexOfMinRMSdistance = ListOfRootMeanSquares.index(MinimumRMSdistance)
    ParameterSetWithMinRMSdistance = ListOfManyDictionariesOfParameters[IndexOfMinRMSdistance]
    
    print()
    print("The list of RMS values determined is:")
    print(ListOfRootMeanSquares)
    print()
    print("The parameter set which minimizes the RMS distance from the data is number " + str(IndexOfMinRMSdistance) + ", corresponding to a RMS of " + str(MinimumRMSdistance) + ", and it has the parameters' values:")
    print(ParameterSetWithMinRMSdistance)
    print()
    
    ############ 4) ############ Verify which is the RMS of the fiducial parameter set, and add corresponding curves to the plots

    FiducialParameterSetRates = deepcopy(StartingParamSetRATES)

    ResultOfComputingRMSfromFeeding_Fiducial = ComputeRMSfeedingForGivenParameterSet(FiducialParameterSetRates, TestParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsFeeding)

    RMS_simulation_Fiducial = ResultOfComputingRMSfromFeeding_Fiducial[0]
    mRNA_HSF_simulation_Fiducial = ResultOfComputingRMSfromFeeding_Fiducial[1]
    mRNA_HSP_simulation_Fiducial = ResultOfComputingRMSfromFeeding_Fiducial[2]
    SimulationFeedingControlsDataRMStimes_Fiducial = ResultOfComputingRMSfromFeeding_Fiducial[3]

    print()
    print("The FIDUCIAL parameter set, i.e. the one we start with, is:\n")
    print(FiducialParameterSetRates)
    print("\nand corresponds to a RMS of\n" + str(RMS_simulation_Fiducial))
    print()

    ############ 5) ############ plot mRNAs(t) to see if it makes sense!
    """
    ### But before, add the curves corresponding to the fiducail parameters set
    CurveNameInLegend = "Fiducial Params Set" + ", RMS = " + str(round(RMS_simulation_Fiducial, 4))
    NormalizeRnaCurvesFromSimulationsToMaxForPlot(RMS_simulation_Fiducial, mRNA_HSF_simulation_Fiducial, mRNA_HSP_simulation_Fiducial,
                                                  ListForPlottingHSF, ListForPlottingHSP, CurveNameInLegend, AllDataControlsFeeding[4])#timeset240minsDataRMS)

    ### And then plot everything
    if len(ListForPlottingHSF) <= 10:
        for key in AllDataControlsFeeding[5]:#ListOfFeedingKeys:
            #PlotSimulationVsDataFeeding(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, AllDataControlsFeeding[4], DictionaryOfListsOfDataHSF[key], DictionaryOfListsOfDataHSP90a[key], DictionaryOfListsOfDataTimes[key], str(key)+"Conrol", FigureExtension, FolderContainingDataVsSimuCalibration)
            PlotSimulationVsDataFeeding(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, AllDataControlsFeeding[4], AllDataControlsFeeding[1][key], AllDataControlsFeeding[2][key], AllDataControlsFeeding[0][key], str(key)+"Conrol", FigureExtension, FolderContainingDataVsSimuCalibration)
    """
    ############ 5bis) ############ plot mRNAs(t) HSP and HSF data vs Model Fit, for paper (compact version)
    ### But before, add the curves corresponding to the fiducail parameters set
    """
    CurveNameInLegend = "Fiducial Params Set" + ", RMS = " + str(round(RMS_simulation_Fiducial, 4))
    NormalizeRnaCurvesFromSimulationsToMaxForPlot(RMS_simulation_Fiducial, mRNA_HSF_simulation_Fiducial, mRNA_HSP_simulation_Fiducial,
                                                  ListForPlottingHSF, ListForPlottingHSP, CurveNameInLegend, AllDataControlsFeeding[4])#timeset240minsDataRMS)

    ### And then plot everything
    #if len(ListForPlottingHSF) <= 10:
    #for key in AllDataControlsFeeding[5]:#ListOfFeedingKeys:
            #PlotSimulationVsDataFeeding(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, AllDataControlsFeeding[4], DictionaryOfListsOfDataHSF[key], DictionaryOfListsOfDataHSP90a[key], DictionaryOfListsOfDataTimes[key], str(key)+"Conrol", FigureExtension, FolderContainingDataVsSimuCalibration)
    PlotSimulationVsDataFeedingModelVSFittedData(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, AllDataControlsFeeding[4], AllDataControlsFeeding[1], AllDataControlsFeeding[2], AllDataControlsFeeding[0], " Conrol", FigureExtension, FolderContainingDataVsSimuCalibration)
    """

    ############ 6) ############ NOW PUT RMS VALUES AND PARAMETER VALUES IN A FILE, FOR SUBSEQUENT USE.

    NumberOfParameters = len(StartingParamSetRATES) # + len(OtherParamSet)

    OutputFile = open(FolderContainingCsvFiles + NameOfOutputFileRMSmanyParamsSets, 'w')

    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
        TotalNumberOfParameterSets = NumberOfRandomSets
    elif SwitchRandomSetsOrParametersK1by1Sets == "ParametersK1by1":
        TotalNumberOfParameterSets = (NumberOfValuesForEachParameterk+1)*NumberOfParameters
    else:
        print("Error N 2 in the value of the switch SwitchRandomSetsOrParametersK1by1Sets!!!")

    for Index in range(TotalNumberOfParameterSets):
        OutputFile.write(str(ListOfRootMeanSquares[Index]) + " " + str(Index))
        for key in ListOfManyDictionariesOfParameters[Index]:
            OutputFile.write(" " + str(ListOfManyDictionariesOfParameters[Index][key]))
        OutputFile.write("\n")

    OutputFile.close()

    # Save in a separate file the order of the keys in the file containing all dictionaries
    OutputFile2 = open(FolderContainingCsvFiles + NameOfOutputFileKeys, 'w')
    for key in ListOfManyDictionariesOfParameters[Index]:
        OutputFile2.write(str(key) + " ")
    OutputFile2.close()

















def PlotResultOfBestFitToData(FolderContainingDataVsSimuCalibration, FINALParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, FigureExtension):

    ############ 4) ############ Verify which is the RMS of the fiducial parameter set, and add corresponding curves to the plots

    FINALParameterSetRates = deepcopy(FINALParamSetRATES)

    ResultOfComputingRMSfromFeeding_FINAL = ComputeRMSfeedingForGivenParameterSet(FINALParameterSetRates, TestParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsFeeding)

    RMS_simulation_FINAL = ResultOfComputingRMSfromFeeding_FINAL[0]
    mRNA_HSF_simulation_FINAL = ResultOfComputingRMSfromFeeding_FINAL[1]
    mRNA_HSP_simulation_FINAL = ResultOfComputingRMSfromFeeding_FINAL[2]
    SimulationFeedingControlsDataRMStimes_FINAL = ResultOfComputingRMSfromFeeding_FINAL[3]

    print()
    print("The FIDUCIAL parameter set, i.e. the one we start with, is:\n")
    print(FINALParameterSetRates)
    print("\nand corresponds to a RMS of\n" + str(RMS_simulation_FINAL))
    print()

    ListForPlottingHSF = []
    ListForPlottingHSP = []
    ############ 5bis) ############ plot mRNAs(t) HSP and HSF data vs Model Fit, for paper (compact version)
    ### But before, add the curves corresponding to the fiducail parameters set
    CurveNameInLegend = "Final Params Set" + ", RMS = " + str(round(RMS_simulation_FINAL, 4))
    NormalizeRnaCurvesFromSimulationsToMaxForPlot(RMS_simulation_FINAL, mRNA_HSF_simulation_FINAL, mRNA_HSP_simulation_FINAL,
                                                  ListForPlottingHSF, ListForPlottingHSP, CurveNameInLegend, AllDataControlsFeeding[4])#timeset240minsDataRMS)

    ### And then plot everything
    #if len(ListForPlottingHSF) <= 10:
    #for key in AllDataControlsFeeding[5]:#ListOfFeedingKeys:
            #PlotSimulationVsDataFeeding(SimulationFeedingControlsDataRMStimes, ListForPlottingHSF, ListForPlottingHSP, AllDataControlsFeeding[4], DictionaryOfListsOfDataHSF[key], DictionaryOfListsOfDataHSP90a[key], DictionaryOfListsOfDataTimes[key], str(key)+"Conrol", FigureExtension, FolderContainingDataVsSimuCalibration)
    PlotSimulationVsDataFeedingModelVSFittedData(SimulationFeedingControlsDataRMStimes_FINAL, ListForPlottingHSF, ListForPlottingHSP, AllDataControlsFeeding[4], AllDataControlsFeeding[1], AllDataControlsFeeding[2], AllDataControlsFeeding[0], "Conrol", FigureExtension, FolderContainingDataVsSimuCalibration)
    

















def PlotRMSvaluesAsFunctionOfParametersFromFile(FolderContainingCsvFiles, FolderContaining1ParametrsRMSplots, FolderContaining2ParametrsRMSplots, FileNameManyParamsSetsRMS, FileNameKeysNamesParamsSets, MyNumberOfBestRMSparamsSetsPlotted, StartingParamSetRATES, SwitchRandomSetsOrParametersK1by1Sets, FigureExtension, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding):

    """BIG function 2 of calibration: Plots RMS values as function of parameters, from file"""

    ############ 1) ############ plot RMS as function of param value, for every param

    NumberOfParameters = len(StartingParamSetRATES) # + len(OtherParamSet)

    ListOfOutputArrays = []
    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
        NumberOfBestRMSparamsSetsPlotted = MyNumberOfBestRMSparamsSetsPlotted # 5000
        SortDataFileByValuesInColonN(FileNameManyParamsSetsRMS, 2+NumberOfParameters, 0, FolderContainingCsvFiles + 'ORDEREDOutputFileRMSmanyParamsSets.csv')
        ExtractFirstNLinesOfFileInputIntoFileOutput(FolderContainingCsvFiles + 'ORDEREDOutputFileRMSmanyParamsSets.csv', FolderContainingCsvFiles + 'CutORDEREDOutputFileRMSmanyParamsSets.csv', NumberOfBestRMSparamsSetsPlotted)
        InvertLinesOrderInFile(FolderContainingCsvFiles + 'CutORDEREDOutputFileRMSmanyParamsSets.csv', FolderContainingCsvFiles + 'InvertedCutORDEREDOutputFileRMSmanyParamsSets.csv')
        FromDataFileToArrays(FolderContainingCsvFiles + 'InvertedCutORDEREDOutputFileRMSmanyParamsSets.csv', 2+NumberOfParameters, ListOfOutputArrays)
    elif SwitchRandomSetsOrParametersK1by1Sets == "ParametersK1by1":
        FromDataFileToArrays(FileNameManyParamsSetsRMS, 2+NumberOfParameters, ListOfOutputArrays)
    else:
        print("Error here!")

    # Prepare Lebels for the plots with the model's parameters on the axes
    Datafile = open(FileNameKeysNamesParamsSets, 'r')
    DataLines = Datafile.readlines()
    Datafile.close()

    MyList = []
    for line in DataLines:
        SplittedLine = line.split()
        for j in range(0, NumberOfParameters):
            MyList.append(SplittedLine[j])

    ListOfParametersNames = deepcopy(MyList)
    ListOfParametersNamesForLegend = []
    DictionaryParamNamesLabels = {"kP0":r"$k_P$ ($(\mu M$ $s)^{-1}$)", "kP0p":r"$k'_P$ ($s^{-1}$)", "kS":r"$k_S$ ($s^{-1}$)", "kSp0":r"$k'_S$ ($s^{-1}$)", "kFp0":r"$k'_F$ ($s^{-1}$)", "kF0":r"$k_F$ ($(\mu M$ $s)^{-1}$)", "kFpi0":r"$k_{\pi_{F}}$ ($s^{-1}$)", "kFGp":r"$k'_{FG}$ ($s^{-1}$)", "kFG":r"$k_{FG}$ ($(\mu M$ $s)^{-1}$)", "ketaF":r"$d_F$ ($s^{-1}$)", "kFsG":r"$k_{F^*G}$ ($(\mu M$ $s)^{-1}$)", "kFsGp":r"$k'_{F^*G}$ ($s^{-1}$)", "kFsp":r"$k'_{F^*}$ ($s^{-1}$)", "kFs":r"$k_{F^*}$ ($s^{-1}$)", "kpiRF":r"$k_{\pi_{RF}}$ ($s^{-1}$)", "kpiRH":r"$k_{\pi_{RH}}$ ($s^{-1}$)", "kpiHP":r"$k_{\pi_{HP}}$ ($s^{-1}$)", "ketaHP":r"$d_{HP}$ ($s^{-1}$)", "ketaRF":r"$d_{RF}$ ($s^{-1}$)", "ketaRHP":r"$d_{RP}$ ($s^{-1}$)"}
    for i in range(len(ListOfParametersNames)):
        for keyParamName in DictionaryParamNamesLabels.keys():
            if ListOfParametersNames[i] == keyParamName:
                Label = DictionaryParamNamesLabels[keyParamName]
        ListOfParametersNamesForLegend.append(Label)

    ArraysToPlot = [["RMS", ListOfOutputArrays[0]]]
    """
    for i in range(NumberOfParameters):
        SingleScatterPlot(figure(), ListOfOutputArrays[2+i], ArraysToPlot, ListOfParametersNamesForLegend[i], 0, 0, "Root Mean Square (adim.)", 0, 0, 'upper right', FolderContaining1ParametrsRMSplots + "RMSparam"+str(i)+FigureExtension)

    ############ 2) ############ do plots with RMS as function of 2 params to study correlation!

    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
        for i in range(NumberOfParameters):
            for j in range(NumberOfParameters):
                fig = figure()
                plt.scatter(ListOfOutputArrays[2+i], ListOfOutputArrays[2+j], c=ListOfOutputArrays[0], cmap=plt.cm.rainbow_r)
                cbar = plt.colorbar()
                cbar.set_label('Root Mean Square (adim.)')
                plt.xlabel(ListOfParametersNamesForLegend[i], fontsize = 18)
                plt.ylabel(ListOfParametersNamesForLegend[j], fontsize = 18)
                PlotAndSave(fig, FolderContaining2ParametrsRMSplots + "TwoParamsRMS_" + str(i) + "_" + str(j), "S", 0, 1)
                plt.close()
    """
    ############ 3) ############ do 1 big cumulative plot with all the single parameter plots
    
    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets" or SwitchRandomSetsOrParametersK1by1Sets == "ParametersK1by1" :

        Nlines = 2
        Ncols = 10
        fig, axes = plt.subplots(nrows=Nlines, ncols=Ncols)
        fig.subplots_adjust(hspace=0.10, wspace=0.10)

        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            # Set up ticks only on one side for the "edge" subplots...
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
            if ax.is_last_col():
                ax.yaxis.set_ticks_position('right')
                ax.yaxis.set_label_position('right')
            if ax.is_first_row():
                ax.xaxis.set_ticks_position('top')
                ax.xaxis.set_label_position('top')
            if ax.is_last_row():
                ax.xaxis.set_ticks_position('bottom')

        ListOfParamsValuesToPlot = []

        for i in range(Nlines):
            for j in range(Ncols):
                print(" ")
                print(i)
                print(" ")
                if i == 0:
                    k = j
                elif i == 1:
                    k = j + 10

                ###############################################

                print("\nSTARTING TO CHARGE PARAMTER SET FROM FILE...JUST TO ADD TO SCATTER PLOT..\n")
                # Read the whole file into a variable which is a list of every row of the file.
                #if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
                #    Datafile = open('OutputFileBestParametersSet.csv', 'r')
                #elif SwitchRandomSetsOrParametersK1by1Sets == "ParametersK1by1":
                #    Datafile = open('./CsvFilesWithRMSandParamSets/OutputFileRMSmanyParamsSets.csv', 'r')
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

                print("The parameter set charged from the txt file is:")
                print(BestParameterSetFromGradientSearchFromFile)
                
                #if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets": 
                ThisParametrSet = deepcopy(BestParameterSetFromGradientSearchFromFile)
                RMSFeedingList = ComputeRMSfeedingForGivenParameterSet(ThisParametrSet, DefaultParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsFeeding)
                RMSFeedingBestParSetGradSearch = RMSFeedingList[0]
                print("look here!!!")
                print(RMSFeedingBestParSetGradSearch)

                CurrentParameterNameLaTeX = ListOfParametersNamesForLegend[k]
                print(CurrentParameterNameLaTeX)

                DictionaryConvertsinParNamesLaTeXToNormal = {"$k'_{F^*}$ ($s^{-1}$)":'kFsp', '$d_F$ ($s^{-1}$)':'ketaF', '$k_F$ ($(\\mu M$ $s)^{-1}$)':'kF0', '$k_P$ ($(\\mu M$ $s)^{-1}$)':'kP0', '$d_{HP}$ ($s^{-1}$)':'ketaHP', '$k_{\\pi_{HP}}$ ($s^{-1}$)':'kpiHP', "$k'_F$ ($s^{-1}$)":'kFp0', '$k_{\\pi_{F}}$ ($s^{-1}$)':'kFpi0', '$k_{FG}$ ($(\\mu M$ $s)^{-1}$)':'kFG', "$k'_S$ ($s^{-1}$)":'kSp0', '$k_{\\pi_{RF}}$ ($s^{-1}$)':'kpiRF', '$k_{F^*}$ ($s^{-1}$)':'kFs', '$d_{RF}$ ($s^{-1}$)':'ketaRF', '$k_S$ ($s^{-1}$)':'kS', '$d_{RP}$ ($s^{-1}$)':'ketaRHP', "$k'_{FG}$ ($s^{-1}$)":'kFGp', '$k_{\\pi_{RH}}$ ($s^{-1}$)':'kpiRH', "$k'_P$ ($s^{-1}$)":'kP0p', '$k_{F^*G}$ ($(\\mu M$ $s)^{-1}$)':'kFsG', "$k'_{F^*G}$ ($s^{-1}$)":'kFsGp'}

                CurrentParameterNameNormal = DictionaryConvertsinParNamesLaTeXToNormal[CurrentParameterNameLaTeX]
                print(CurrentParameterNameNormal)

                #if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets": 
                CurrentParamValueFromBestParSetGradSearch = BestParameterSetFromGradientSearchFromFile[CurrentParameterNameNormal]
                print(CurrentParamValueFromBestParSetGradSearch)
                ListOfParamsValuesToPlot.append(CurrentParamValueFromBestParSetGradSearch)
                #print(ListOfParamsValuesToPlot)

                #RMSlist = [RMSFeedingBestParSetGradSearch]*len(ListOfParamsValuesToPlot)
                #print(RMSlist)

                CurrentParamFiducialValue = StartingParamSetRATES[CurrentParameterNameNormal]

                ##############################################

                #if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
                #    MyColor = 'blue'
                #elif SwitchRandomSetsOrParametersK1by1Sets == "ParametersK1by1":
                #    ThisArray = ListOfOutputArrays[2+k]
                #    CentralValueOfThisArray = ThisArray[(len(ThisArray)-1)/2.0]
                #    if CurrentParamFiducialValue == CentralValueOfThisArray :
                #        MyColor = 'green'
                #    else:
                #        MyColor = 'blue'

                MyXarray = ListOfOutputArrays[2+k]
                MyYarray = ListOfOutputArrays[0]
                #MyColor = 'blue'

                ll = 0
                TempArrayX = []
                TempArrayY = []
                for Datum in MyXarray:
                    if Datum != CurrentParamFiducialValue:
                        TempArrayX.append(Datum)
                        TempArrayY.append(MyYarray[ll])
                    ll = ll+1
                MyXarrayParamPERTURBED = TempArrayX
                MyYarrayParamPERTURBED = TempArrayY
                MyColorParamPERTURBED = 'blue'

                #axes[i,j].scatter(MyXarray, MyYarray, s=18, edgecolor = '', color = MyColor) # Here!!!!!!! 
                axes[i,j].scatter(MyXarrayParamPERTURBED, MyYarrayParamPERTURBED, s=18, edgecolor = '', color = MyColorParamPERTURBED) # Here!!!!!!! 
                axes[i,j].axvline(x=CurrentParamFiducialValue, color = 'red', linewidth = 1.5)    
                if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
                    axes[i,j].scatter(CurrentParamValueFromBestParSetGradSearch, RMSFeedingBestParSetGradSearch, s=400, color = 'yellow', marker = "*", edgecolor = 'black', linewidths = 2) 

                axes[i,j].set_xlim(min(ListOfOutputArrays[2+k]),max(ListOfOutputArrays[2+k]))
                axes[i,j].set_ylim(min(ListOfOutputArrays[0]),max(ListOfOutputArrays[0]))
                #plt.setp(axes[i,j].get_xticklabels(), rotation=45, horizontalalignment='left')

                axes[0,j].xaxis.set_visible(True)
                axes[0,j].locator_params(axis='x',nbins=5)
                plt.setp(axes[0,j].get_xticklabels(), rotation=90, horizontalalignment='left')
                if i == 0:
                    axes[0,j].set_xlabel(ListOfParametersNamesForLegend[k], fontsize='x-large')
                axes[0,j].get_xaxis().set_tick_params(direction='out')

                axes[Nlines-1,j].xaxis.set_visible(True)
                axes[Nlines-1,j].locator_params(axis='x',nbins=5)
                plt.setp(axes[Nlines-1,j].get_xticklabels(), rotation=90, horizontalalignment='right')
                if i == 1:
                    axes[Nlines-1,j].set_xlabel(ListOfParametersNamesForLegend[k], fontsize='x-large')
                axes[Nlines-1,j].get_xaxis().set_tick_params(direction='out')

                axes[i,0].yaxis.set_visible(True)
                axes[i,0].locator_params(axis='y',nbins=8)
                axes[i,0].set_ylabel("Root Mean Square (adim.)", fontsize='x-large')
                axes[i,0].get_yaxis().set_tick_params(direction='out')

                axes[i,Ncols-1].yaxis.set_visible(True)
                axes[i,Ncols-1].locator_params(axis='y',nbins=8)
                axes[i,Ncols-1].set_ylabel("Root Mean Square (adim.)", fontsize='x-large')
                axes[i,Ncols-1].get_yaxis().set_tick_params(direction='out')
                print("k is ")
                print(k)
        print(ListOfParametersNamesForLegend)

        plt.show()




    ############ 4) ############ do 1 big cumulative plot with all the plots above

    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":

        fig, axes = plt.subplots(nrows=NumberOfParameters, ncols=NumberOfParameters)
        fig.subplots_adjust(hspace=0.00, wspace=0.00)

        for ax in axes.flat:
            # Hide all ticks and labels
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            # Set up ticks only on one side for the "edge" subplots...
            if ax.is_first_col():
                ax.yaxis.set_ticks_position('left')
            if ax.is_last_col():
                ax.yaxis.set_ticks_position('right')
                ax.yaxis.set_label_position('right')
            if ax.is_first_row():
                ax.xaxis.set_ticks_position('top')
                ax.xaxis.set_label_position('top')
            if ax.is_last_row():
                ax.xaxis.set_ticks_position('bottom')

        for i in range(NumberOfParameters):
            for j in range(NumberOfParameters):
                DotSize = 14
                if i>j:
                    if (i % 2) != 0:
                        if (j % 2) != 0:
                            im = axes[i,j].scatter(ListOfOutputArrays[2+j], ListOfOutputArrays[2+i], s=DotSize, c=ListOfOutputArrays[0],
         edgecolor = '', cmap=plt.cm.rainbow_r)
                    elif (i % 2) == 0:
                        if  (j % 2) == 0:
                            im = axes[i,j].scatter(ListOfOutputArrays[2+j], ListOfOutputArrays[2+i], s=DotSize, c=ListOfOutputArrays[0],  edgecolor = '', cmap=plt.cm.rainbow_r)
                elif i<j:
                    if (i % 2) != 0:
                        if (j % 2) == 0:
                            im = axes[i,j].scatter(ListOfOutputArrays[2+j], ListOfOutputArrays[2+i], s=DotSize, c=ListOfOutputArrays[0], edgecolor = '', cmap=plt.cm.rainbow_r)
                    elif (i % 2) == 0:
                        if (j % 2) != 0:
                            im = axes[i,j].scatter(ListOfOutputArrays[2+j], ListOfOutputArrays[2+i], s=DotSize, c=ListOfOutputArrays[0], edgecolor = '', cmap=plt.cm.rainbow_r)
                elif i==j:
                    # Label the diagonal subplots...
                    axes[i,j].annotate(ListOfParametersNamesForLegend[i].split(" ", 1)[0], (0.5, 0.5), xycoords='axes fraction', ha='center', va='center',)
                    #axes[i,j].scatter(ListOfOutputArrays[2+i], ListOfOutputArrays[0], s=8)               

                axes[i,j].set_xlim(min(ListOfOutputArrays[2+j]),max(ListOfOutputArrays[2+j]))
                axes[i,j].set_ylim(min(ListOfOutputArrays[2+i]),max(ListOfOutputArrays[2+i]))
                plt.setp(axes[i,j].get_xticklabels(), rotation=45, horizontalalignment='left')

        cbar = fig.colorbar(im, ax=axes.ravel().tolist())
        cbar.set_label('Root Mean Square (adim.)')

        for k in range(round(NumberOfParameters/2)):

            axes[0,2*k+1].xaxis.set_visible(True)
            axes[0,2*k+1].locator_params(axis='x',nbins=3)
            plt.setp(axes[0,2*k+1].get_xticklabels(), rotation=45, horizontalalignment='left')
            axes[0,2*k+1].set_xlabel(ListOfParametersNamesForLegend[2*k+1])
            axes[0,2*k+1].get_xaxis().set_tick_params(direction='out')

            axes[NumberOfParameters-1,2*k].xaxis.set_visible(True)
            axes[NumberOfParameters-1,2*k].locator_params(axis='x',nbins=3)
            plt.setp(axes[NumberOfParameters-1,2*k].get_xticklabels(), rotation=45, horizontalalignment='right')
            axes[NumberOfParameters-1,2*k].set_xlabel(ListOfParametersNamesForLegend[2*k])
            axes[NumberOfParameters-1,2*k].get_xaxis().set_tick_params(direction='out')

            axes[2*k,0].yaxis.set_visible(True)
            axes[2*k,0].locator_params(axis='y',nbins=3)
            axes[2*k,0].set_ylabel(ListOfParametersNamesForLegend[2*k])
            axes[2*k,0].get_yaxis().set_tick_params(direction='out')

            axes[2*k+1,NumberOfParameters-1].yaxis.set_visible(True)
            axes[2*k+1,NumberOfParameters-1].locator_params(axis='y',nbins=3)
            axes[2*k+1,NumberOfParameters-1].set_ylabel(ListOfParametersNamesForLegend[2*k+1])
            axes[2*k+1,NumberOfParameters-1].get_yaxis().set_tick_params(direction='out')
        
        plt.show()
        #PlotAndSave(fig, "MCexplorationALLinONE" + FigureExtension, "PS", 0, 0)





    ############ 5) ############ Identify the parameters set in OutputFile wich minimizes the RMS
    if SwitchRandomSetsOrParametersK1by1Sets == "RandomSets":
        Datafile = open(FolderContainingCsvFiles + "CutORDEREDOutputFileRMSmanyParamsSets.csv", 'r')
        DataLines = Datafile.readlines()
        Datafile.close()
        i = 0

        for line in DataLines:
            if not line.startswith("#"):
                SplittedLine = line.split()
                break
            i = i + 1 

        BestRMSDictionaryParameters = {}
        for i in range(len(SplittedLine)-2):
                BestRMSDictionaryParameters.update({deepcopy(ListOfParametersNames[i]) : float(SplittedLine[i+2])})

        print()
        print("RMS best = " + str(SplittedLine[0]))
        print()
        print(BestRMSDictionaryParameters)
        print()






def ForBestRMSfeedingPointspointsComputeRMSDoubleHS(FolderRMS1vs2, FolderContainingCsvFiles, NameFigureRMS1vs2, FileNameManyParamsSetsRMSDoubleHS, StartingParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsDoubleHS, FigureExtension):

    """BIG function 3 of calibration: For the best 5000 points wrt RSM Feeding, compute the corresponding RMS w.r.t. Double HS"""

    ############ 1) ############ Select datafile with parameters sets and RMS values, extract it and put in arrays

    ListOfOutputArraysFromDatatfile = []
    FromDataFileToArrays(FileNameManyParamsSetsRMSDoubleHS, 2+20, ListOfOutputArraysFromDatatfile)

    FileNameKeysNamesParamsSetsDoubleHS = FolderContainingCsvFiles + 'OutputFileKeys.csv'
    ListOfKeys = []
    FromDataFileToArrays(FileNameKeysNamesParamsSetsDoubleHS, 20, ListOfKeys, Strings="Yes")

    DictionaryOfNewParameterSetsDictionariesFromTheNumberOfTheBest = {}
    for i in range(len(ListOfOutputArraysFromDatatfile[0])):
        DictionaryForParameters = {}
        for j in range(20):
            NewKeyList = deepcopy(ListOfKeys[j])
            NewKey = NewKeyList[0]
            NewParameterValue = deepcopy(ListOfOutputArraysFromDatatfile[j+2][i])
            DictionaryForParameters.update({NewKey:NewParameterValue})
        DictionaryOfNewParameterSetsDictionariesFromTheNumberOfTheBest.update({ str(ListOfOutputArraysFromDatatfile[0][i]) : DictionaryForParameters })

    ############ 2) ############ NOW NEED TO LOOP OVER EVERY PARAMETER SET

    #### Prepare list to be filled with RMS values
    ListOfRootMeanSquaresFeeding = []
    ListOfRootMeanSquaresDoubleHS = []

    i = 0
    for key in DictionaryOfNewParameterSetsDictionariesFromTheNumberOfTheBest:
        ParamSetRates = DictionaryOfNewParameterSetsDictionariesFromTheNumberOfTheBest[key]
        RMSfromDoubleHS = ComputeRMSdoubleHSforGivenParameterSet(ParamSetRates, TestParamSetForREACTIONS, DefaultParamSetInitCond, "Yes", AllDataControlsDoubleHS)
        print()
        print(ParamSetRates)
        print()
        print(RMSfromDoubleHS)
        print()
        print(float(key))
        print()

        ListOfRootMeanSquaresFeeding.append(float(key))
        ListOfRootMeanSquaresDoubleHS.append(RMSfromDoubleHS)
    
    ############ 3) ############ NOW PLOT THE RESULTS
    ArraysToPlot = [["Best prameters sets", ListOfRootMeanSquaresDoubleHS]]
    SingleScatterPlot(figure(), ListOfRootMeanSquaresFeeding, ArraysToPlot, "RMS Feeding (adim.)", 0, 0, "RMS double HS (adim.)", 0, 0, 'upper right', FolderRMS1vs2 + NameFigureRMS1vs2 + FigureExtension)






def ExecuteGradientSearch(FolderContainingGradientSearchPlots, UseRMSForFeedingOrTotal, MaxNumberOfIterations, MyNumberOfIterationsForAverage, ThresholdAverageRMSdecrease, IncrementInComputingDerivative, NameOutputFileBestParametersSet, StartingParamSetRATES, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, GammaMin, GammaMax, GammaBisectionStep, FigureExtension):

    """BIG function 4 of calibration: Implement the gradient search"""

    ############ 0) ############ Check if ok preconditioning of RMS, i.e. rescale all parameters by their fiducial value for better applicability of gradient search

    RescalingFactorsDictionary = deepcopy(StartingParamSetRATES)
    ORIGINAL_ParameterSetDictionary1 = deepcopy(StartingParamSetRATES)

    ############ 1) ############ Take the fiducial parameter set as the starting point in the parameter space

    ListOfIndexesOfGammaValuesChosenByLineSearch = []

    ListOfParametersSetsKs = []
    ListOfRootMeanSquares = []

    StartingParameterSetKs = Convert_ORIGINAL_to_RESCALED_ParameterSet(ORIGINAL_ParameterSetDictionary1, RescalingFactorsDictionary)
    ListOfParametersSetsKs.append(StartingParameterSetKs)

    DictionaryOfListsOfParametersValuesForPlotting = {}
    for key in StartingParameterSetKs:
        DictionaryOfListsOfParametersValuesForPlotting.update({key:[deepcopy(StartingParameterSetKs[key])]})

    ############ 2) ############ For every point in the parameter space do the following
    i = 0
    ListOfAverageRMSdecreaseOverLastTenIterations = []
    for i in range(MaxNumberOfIterations):  # 150

        ##### A: Compute the RMS corresponding to the current point in the parameter space
        CurrentParamsSet = deepcopy(ListOfParametersSetsKs[i])
        if UseRMSForFeedingOrTotal == "Feeding":
            RMS = ComputeRMSfeedingForGivenParameterSet_RESCALED_PARAMETERS(CurrentParamsSet, TestParamSetForREACTIONS, DefaultParamSetInitCond, "No", AllDataControlsFeeding, RescalingFactorsDictionary)
        elif UseRMSForFeedingOrTotal == "FeedingPlusDouble":
            RMS = ComputeRMStotalForGivenParameterSet_RESCALED_PARAMETERS(CurrentParamsSet, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, RescalingFactorsDictionary, AllDataControlsDoubleHS)
        else :
            print("Error1 in gradient search")

        ListOfRootMeanSquares.append(RMS)

        ##### F: Verify if the stop condition is fullfilled or not

        NumberOfIterationsForAverage = MyNumberOfIterationsForAverage # 10

        if i >= (NumberOfIterationsForAverage-1):
            SumRMSdecreases = 0.
            ListOfAverages = []
            print("HERE " + str(ListOfRootMeanSquares))
            for n in range(NumberOfIterationsForAverage-1):
                RMSdecrease = abs(deepcopy(ListOfRootMeanSquares[i-n]) - deepcopy(ListOfRootMeanSquares[i-n-1]))
                SumRMSdecreases = SumRMSdecreases + deepcopy(RMSdecrease)
                print(str(i) + "   " + str(n) + "   " + str(RMSdecrease) + "   " + str(SumRMSdecreases))
            AverageRMSdecrease = deepcopy(SumRMSdecreases)/deepcopy(NumberOfIterationsForAverage-1)
            ListOfAverages.append(AverageRMSdecrease)
            print(AverageRMSdecrease)
            if AverageRMSdecrease < ThresholdAverageRMSdecrease: # 0.00001
                OutputFileAverages = open('OutputBREAKINGgradient.txt', 'w')
                OutputFileAverages.write(str(ListOfRootMeanSquares) + "\n")
                OutputFileAverages.write(str(ListOfAverages) + "\n")
                OutputFileAverages.write(str(AverageRMSdecrease) + "\n")
                OutputFileAverages.close()
                break
        elif i < (NumberOfIterationsForAverage-1):
            pass
        else:
            print("Error again in gradient search")

        ##### B: Initialize the gradient as an empty dictionary
        GradientOfRMS = {}

        ##### C: Fill in the gradient by iterating over all the parameters corresponding to the current point
        for key in CurrentParamsSet:

            ### C.1: Compute Dk
            OldValueKey = deepcopy(CurrentParamsSet[key])
            Dkey = IncrementInComputingDerivative # 1.e-6

            ### C.2a: Compute RMS+
            CurrentParamsSet.update({key : OldValueKey + Dkey})
            CurrentParamsSetPlus = deepcopy(CurrentParamsSet)
            if UseRMSForFeedingOrTotal == "Feeding":
                RMSplus = ComputeRMSfeedingForGivenParameterSet_RESCALED_PARAMETERS(CurrentParamsSetPlus, TestParamSetForREACTIONS, DefaultParamSetInitCond, "No", AllDataControlsFeeding, RescalingFactorsDictionary)
            elif UseRMSForFeedingOrTotal == "FeedingPlusDouble":
                RMSplus = ComputeRMStotalForGivenParameterSet_RESCALED_PARAMETERS(CurrentParamsSetPlus, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, RescalingFactorsDictionary, AllDataControlsDoubleHS)
            else :
                print("Error2 in gradient search")

            ### C.2b: Compute RMS-
            CurrentParamsSet.update({key : OldValueKey - Dkey})
            CurrentParamsSetMinus = deepcopy(CurrentParamsSet)
            if UseRMSForFeedingOrTotal == "Feeding":
                RMSminus = ComputeRMSfeedingForGivenParameterSet_RESCALED_PARAMETERS(CurrentParamsSetMinus, TestParamSetForREACTIONS, DefaultParamSetInitCond, "No", AllDataControlsFeeding, RescalingFactorsDictionary)
            elif UseRMSForFeedingOrTotal == "FeedingPlusDouble":
                RMSminus = ComputeRMStotalForGivenParameterSet_RESCALED_PARAMETERS(CurrentParamsSetMinus, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, RescalingFactorsDictionary, AllDataControlsDoubleHS)
            else :
                print("Error2 in gradient search")

            ### C.2c: Put back to original the Param set
            CurrentParamsSet.update({key : OldValueKey})

            ### C.3: Compute the approximate partial derivative of RMS w.r.t. the current parameter
            PartialDerivRMSwrtKeyApprox = (RMSplus - RMSminus) / (2*Dkey)
            print(str(key)+"  "+str(OldValueKey)+"  "+str(Dkey)+"  "+str(RMS)+"  "+str(RMSplus)+"  "+str(RMSminus)+"  "+str(PartialDerivRMSwrtKeyApprox))

            ### C.4: Append this partial derivative to the gradient
            GradientOfRMS.update({key : PartialDerivRMSwrtKeyApprox})


        ##### D: Choose the value of the multiplier which determines the lenght of the next step

        def MoveInParameterSpaceInDirectionOfGradientByStepOfLenghtGamma(OldParamsSet, GradientRMS, gamma):
            #Starting from OldParamsSet it provides a new ParamsSet moving in the direction of the gradient for a step of lenght gamma
            NewParamsSet = deepcopy(OldParamsSet)
            for key in OldParamsSet:
                ### G.1: Extract approximated partial derivative from gradient
                PartialRMSkey = deepcopy(GradientRMS[key])
                ### G.2: Change the value of one parameter using the partial derivative
                NewValueKey = deepcopy(OldParamsSet[key]) - gamma * PartialRMSkey
                NewParamsSet.update({key:NewValueKey})
            return deepcopy(NewParamsSet)

        def RMSFunctionToBeMinimizedWRTgamma(gamma):
            #RMS function computed along the direction of the gradient, as a function of the step lenght gamma, for line search minimization
            TestGammaParamsSet = MoveInParameterSpaceInDirectionOfGradientByStepOfLenghtGamma(CurrentParamsSet, GradientOfRMS, gamma)
            if UseRMSForFeedingOrTotal == "Feeding":
                RMS = ComputeRMSfeedingForGivenParameterSet_RESCALED_PARAMETERS(TestGammaParamsSet, TestParamSetForREACTIONS, DefaultParamSetInitCond, "No", AllDataControlsFeeding, RescalingFactorsDictionary)
            elif UseRMSForFeedingOrTotal == "FeedingPlusDouble":
                RMS = ComputeRMStotalForGivenParameterSet_RESCALED_PARAMETERS(TestGammaParamsSet, TestParamSetForREACTIONS, DefaultParamSetInitCond, AllDataControlsFeeding, RescalingFactorsDictionary, AllDataControlsDoubleHS)
            else :
                print("Error2 in gradient search")
            return deepcopy(RMS)

        GammaApproxMinimizingRMS = FindMinimumOfFunctionUsingGoldenRatioBisectionMethod(RMSFunctionToBeMinimizedWRTgamma, GammaMin, GammaMax, GammaBisectionStep)
        ListOfIndexesOfGammaValuesChosenByLineSearch.append(GammaApproxMinimizingRMS)

        print("\nGamma = "+str(GammaApproxMinimizingRMS)+"\n")

        ##### G: Compute the next point of the parameter space, i.e. iterate over all the parameters to:

        ### G.3: Update the dictionary representing the point in the parameter space
        NextParameterSetKs = MoveInParameterSpaceInDirectionOfGradientByStepOfLenghtGamma(CurrentParamsSet, GradientOfRMS, GammaApproxMinimizingRMS)

        for key in CurrentParamsSet:
            DictionaryOfListsOfParametersValuesForPlotting[key].append(NextParameterSetKs[key])

        ListOfParametersSetsKs.append(NextParameterSetKs)

        ##### H: Increase loop index by 1.
        i = i + 1

    FinalParamSetMinimizingRMS = deepcopy(ListOfParametersSetsKs[i-1])

    print("\n"+str(ListOfParametersSetsKs))
    print("\n"+str(ListOfRootMeanSquares))
    print("\n"+str(GradientOfRMS))
    print("\n"+str(DictionaryOfListsOfParametersValuesForPlotting))
    print()
    print(FinalParamSetMinimizingRMS)

    ############ 3) ############ Plot RMS as a function of i ( iterations number)
    fig = figure()
    plt.plot( range(len(ListOfRootMeanSquares)), ListOfRootMeanSquares, marker="o", linewidth=1) 
    plt.xlabel("Number of iterations")
    plt.ylabel("RMS")
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    PlotAndSave(fig, FolderContainingGradientSearchPlots + "GradientSearch_RMS" + FigureExtension, "PS", 0, 0)

    ############ 4) ############ Plot Gamma Value (from line search minimization of RMS) as a function of i (iterations number)
    fig = figure()
    plt.plot( range(len(ListOfIndexesOfGammaValuesChosenByLineSearch)), ListOfIndexesOfGammaValuesChosenByLineSearch, marker="o", linewidth=1) 
    plt.xlabel("Number of iterations")
    plt.ylabel("Gamma")
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)
    PlotAndSave(fig, FolderContainingGradientSearchPlots + "GradientSearch_GammaLineSearch" + FigureExtension, "PS", 0, 0)

    ############ 5) ############ Plot 20 plots with k_i as a function of i (parameter value as a function of iterations number)
    for key in StartingParameterSetKs:
        fig = figure()
        plt.plot( range(len(DictionaryOfListsOfParametersValuesForPlotting[key])), DictionaryOfListsOfParametersValuesForPlotting[key], marker="o", linewidth=1) 
        plt.xlabel("Number of iterations")
        plt.ylabel(str(key) + "/" + str(key) + "_FIDUCIAL")
        ax = plt.gca()
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_scientific(False)
        PlotAndSave(fig, FolderContainingGradientSearchPlots + "GradientSearch_" + str(key) + FigureExtension, "PS", 0, 0)

    ############ 6) ############ Rescale back (invert the preconditioning on the RMS, to go back to the original values of the parameters)
    RESCALEDFinalParamSetMinimizingRMS = FinalParamSetMinimizingRMS
    ORIGINALFinalParamSetMinimizingRMS = Convert_RESCALED_to_ORIGINAL_ParameterSet(RESCALEDFinalParamSetMinimizingRMS, RescalingFactorsDictionary)
    print("\n"+str(RESCALEDFinalParamSetMinimizingRMS)+"\n"+str(ORIGINALFinalParamSetMinimizingRMS)+"\n")

    ############ 7) ############ Write the best parameter set so found into a file, for subsequent use
    
    OutputFile = open(NameOutputFileBestParametersSet, 'w')

    for key in ORIGINALFinalParamSetMinimizingRMS:
        OutputFile.write(str(key))
        OutputFile.write(" ")
        OutputFile.write(str(ORIGINALFinalParamSetMinimizingRMS[key]))
        OutputFile.write("\n")

    OutputFile.close()

