from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def AddLetterToSubplot(ax, SubPlotLabel, Xmultiplyer, Ymultiplyer):
    xlimMIN = ax.get_xlim()[0]
    xlimMAX = ax.get_xlim()[1]
    ylimMIN = ax.get_ylim()[0]
    ylimMAX = ax.get_ylim()[1]
    ax.text(xlimMIN + Xmultiplyer * (xlimMAX - xlimMIN), ylimMIN + Ymultiplyer * (ylimMAX - ylimMIN), SubPlotLabel,
            fontdict={'fontsize': "large"})


def SubPlot(ax, time, ArraysToPlot, xLabel, xMin, xMax, yLabel, yMin, yMax, LegendPosition, SubPlotLabel,
            Legendfontsize="x-small", Legendfancybox=True, UsePointMarker="No", LineStyleListInverted="No", Black="No"):
    ColorsList = ['red', 'blue', 'green', 'black', 'orange', 'gray', 'violet', 'yellow', 'cyan', 'chartreuse']

    if LineStyleListInverted == "No":
        LinestylesList = ['-', '--', '-.', ':', '-', '--', '-.', ':']
        LinesizesList = [1., 1.3, 1.6, 2., 1., 1.3, 1.6, 2.]
    elif LineStyleListInverted == "Yes":
        LinestylesList = ['--', '-', '-.', ':', '--', '-', '-.', ':']
        LinesizesList = [1.3, 1., 1.6, 2., 1.3, 1., 1.6, 2.]

    if Black == "Yes":
        ColorsList = ['black']
        LinesizesList = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

    i = 0
    for i in range(len(ArraysToPlot)):
        if UsePointMarker == "No":
            ax.plot(time, ArraysToPlot[i][1], color=ColorsList[i], linewidth=LinesizesList[i], linestyle=LinestylesList[i], label=ArraysToPlot[i][0])
        elif UsePointMarker == "Yes":
            ax.plot(time, ArraysToPlot[i][1], color=ColorsList[i], marker="o", label=ArraysToPlot[i][0])
        else:
            print("Error in SubPlot")
        i = i + 1
    ax.set_xlim(xMin, xMax)
    if (yMin - yMax) != 0:
        ax.set_ylim(yMin, yMax)
    ax.set_xlabel(xLabel, fontsize="x-small")
    ax.set_ylabel(yLabel, fontsize="x-small")
    AddLetterToSubplot(ax, SubPlotLabel, -0.2, 1.065)
    ax.legend(loc=LegendPosition, numpoints=1, fontsize=Legendfontsize, fancybox=Legendfancybox)


def DataSubPlot(ax, time, ListOfOutputArrays, xLabel, xMin, xMax, yLabel, yMin, yMax, LegendPosition, DataLegendList,
                SubPlotLabel, Legendfontsize="x-small", Legendfancybox=True, Black = "No", MyLinestyle='-'):
    MarkersList = [ ".", "o", "s", "^", "v", ">", "<", ",",]
    ColorsList = ['red', 'blue', 'green', 'black', 'orange', 'gray', 'violet', 'yellow', 'cyan']
    if Black == "Yes":
        ColorsList = ['black']
    for k in range(1, len(ListOfOutputArrays)):
        ax.plot(time, ListOfOutputArrays[k], marker=MarkersList[k], color=ColorsList[k - 1],
                label=DataLegendList[k - 1], markersize = 10, linestyle = MyLinestyle)
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    ax.set_xlabel(xLabel, fontsize="medium")
    ax.set_ylabel(yLabel, fontsize="medium")
    AddLetterToSubplot(ax, SubPlotLabel, -0.25, 1.055)
    ax.legend(loc=LegendPosition, numpoints=1, fontsize=Legendfontsize, fancybox=Legendfancybox)


def DataSubPlotMOD(ax, time, ListOfOutputArrays, xLabel, xMin, xMax, yLabel, yMin, yMax, LegendPosition, DataLegendList,
                SubPlotLabel, Legendfontsize="x-small", Legendfancybox=True, ColorNumber=0):
    MarkersList = [ "o", "s", "^", "v", ">", "<", ",", "."]
    ColorsList = ['red', 'blue', 'green', 'orange', 'gray', 'violet', 'yellow', 'cyan']
    for k in range(1, len(ListOfOutputArrays)):
        ax.plot(time, ListOfOutputArrays[k], marker=MarkersList[ColorNumber], color=ColorsList[ColorNumber],
                label=DataLegendList[k - 1], markersize = 12)
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    ax.set_xlabel(xLabel, fontsize="medium")
    ax.set_ylabel(yLabel, fontsize="medium")
    AddLetterToSubplot(ax, SubPlotLabel, -0.25, 1.055)
    ax.legend(loc=LegendPosition, numpoints=1, fontsize=Legendfontsize, fancybox=Legendfancybox)


def PlotAndSave(figure, figureName, PlotSave, TightLayout, Wait):
    """ Plot all and save """

    if TightLayout == 1:
        figure.tight_layout()
    elif TightLayout == 0:
        pass
    else:
        print("Error with TightLayout")

    if PlotSave == "P" or PlotSave == "PS" or PlotSave == "SP":
        if Wait == 1:
            figure.show()
        elif Wait == 0:
            plt.show()
        else:
            print("Error in PlotAndSave")
    if PlotSave == "S" or PlotSave == "PS" or PlotSave == "SP":
        FileNameHSR = figureName
        figure.savefig(FileNameHSR)
    if PlotSave != "P" and PlotSave != "PS" and PlotSave != "SP" and PlotSave != "S":
        print("Error in PlotAndSave")
    
    #plt.clf()
    #plt.cla()
    #plt.close()


def FromDataFileToArrays(DataFileName, NumOfColumns, ListOfOutputArrays, Strings="No"):
    """ Read data from a file and put each column in an array, thus filling a list of arrays """

    # Read the whole file into a variable which is a list of every row of the file.
    Datafile = open(DataFileName, 'r')
    DataLines = Datafile.readlines()
    Datafile.close()

    # Initialize the lists which will contain the data:
    ListOfOutputLists = []
    for i in range(NumOfColumns):
        ListOfOutputLists.append([])

        # Scan the rows of the file stored in lines, and put the values into some variables:
    for line in DataLines:
        if not line.startswith("#"):
            SplittedLine = line.split()
            for j in range(0, NumOfColumns):
                if Strings == "No":
                    ListOfOutputLists[j].append(float(SplittedLine[j]))
                elif Strings == "Yes":
                    ListOfOutputLists[j].append(SplittedLine[j])
                else:
                    print("Problem in FromDataFileToArrays!")

    # Make arrays out of them
    for k in range(0, NumOfColumns):
        ListOfOutputArrays.append(np.array(ListOfOutputLists[k]))


def SortDataFileByValuesInColonN(DataFileName, NumOfColumns, N, OutputFileName):

    # Read the whole file into a variable which is a list of every row of the file.
    Datafile = open(DataFileName, 'r')
    DataLines = Datafile.readlines()
    Datafile.close()

    # Initialize the lists which will contain the data:
    ListOfOutputLists = []
    
    for i in range(len(DataLines)):
        ListOfOutputLists.append([])

        # Scan the rows of the file stored in lines, and put the values into some variables:
    
    LineIndex = 0
    for line in DataLines:
        if not line.startswith("#"):
            SplittedLine = line.split()
            #print(SplittedLine)
            for j in range(0, NumOfColumns):
                ListOfOutputLists[LineIndex].append(float(SplittedLine[j]))
            LineIndex = LineIndex + 1
            #print(LineIndex)
        #if LineIndex >= MaxNumberOfOutputBestLines:
        #    print("Brake at line " + str(LineIndex) + "\n")
        #    break
    
    def getRMSForOrdering(List):
        return List[0]

    ListOfOutputListsORDERED = sorted( deepcopy(ListOfOutputLists), key = getRMSForOrdering)

    #print(ListOfOutputLists)
    #print(ListOfOutputListsORDERED)

    OutputFile = open(OutputFileName, 'w')

    NumberOfLinesNewFile = len(DataLines)

    for Index in range(NumberOfLinesNewFile):
        for i in range(len(ListOfOutputListsORDERED[Index])):
            OutputFile.write(str(ListOfOutputListsORDERED[Index][i]) + " ")
        OutputFile.write("\n")

    OutputFile.close()


def ExtractFirstNLinesOfFileInputIntoFileOutput(InputFileName, OutputFileName, NumberOfLinesExtracted):

    # Read the whole file into a variable which is a list of every row of the file.
    Datafile = open(InputFileName, 'r')
    DataLines = Datafile.readlines()
    Datafile.close()

    # Initialize the lists which will contain the data:
    ListOfOutputLists = []
    
    for i in range(len(DataLines)):
        ListOfOutputLists.append([])
    
    OutputFile = open(OutputFileName, 'w')

    LineIndex = 0
    for line in DataLines:
        if LineIndex < NumberOfLinesExtracted:
            OutputFile.write(str(line)) # + "\n"
            #print(LineIndex)            
            LineIndex = LineIndex + 1

    OutputFile.close()


def InvertLinesOrderInFile(InputFileName, OutputFileName):

    # Read the whole file into a variable which is a list of every row of the file.
    Datafile = open(InputFileName, 'r')
    DataLines = Datafile.readlines()
    Datafile.close()

    # Initialize the lists which will contain the data:
    ListOfOutputLists = []
    
    for i in range(len(DataLines)):
        ListOfOutputLists.append([])
    
    OutputFile = open(OutputFileName, 'w')

    LineIndex = 0
    for line in reversed(DataLines):
            OutputFile.write(str(line)) # + "\n"
            #print(LineIndex)            
            LineIndex = LineIndex + 1

    OutputFile.close()


def SinglePlot(fig, x, ArraysToPlot, xLabel, xMin, xMax, yLabel, yMin, yMax, LegendPosition, figurename, labelssize):

    ColorsList = ['Magenta', 'Green', 'blue', 'red', 'cyan', 'orange', 'gray', 'black']
    i = 0
    for i in range(len(ArraysToPlot)):
        plt.plot(x, ArraysToPlot[i][1], color=ColorsList[i], linewidth=1., label=ArraysToPlot[i][0])
        i = i + 1
    plt.xlim(xMin, xMax)
    if (yMin - yMax) != 0:
        plt.ylim(yMin, yMax)
    plt.xlabel(xLabel, fontsize=labelssize)
    plt.ylabel(yLabel, fontsize=labelssize)

    plt.legend(loc=LegendPosition)
    PlotAndSave(fig, figurename, "PS", 0, 0)


def SingleScatterPlot(fig, x, ArraysToPlot, xLabel, xMin, xMax, yLabel, yMin, yMax, LegendPosition, figurename, labelssize='medium'):

    ColorsList = ['Magenta', 'Green', 'blue', 'red', 'cyan', 'orange', 'gray', 'black']
    i = 0
    for i in range(len(ArraysToPlot)):
        plt.plot(x, ArraysToPlot[i][1], color=ColorsList[i], marker='o', label=ArraysToPlot[i][0], linestyle="None")
        i = i + 1
    if (xMin - xMax) != 0:
        plt.xlim(xMin, xMax)
    if (yMin - yMax) != 0:
        plt.ylim(yMin, yMax)
    plt.xlabel(xLabel, fontsize=labelssize)
    plt.ylabel(yLabel, fontsize=labelssize)

    #plt.legend(loc=LegendPosition)
    PlotAndSave(fig, figurename, "PS", 0, 0)



SetOfReferenceValuesForPlotting = {
"Ptot" : 100001./ 1000.,  # (mM) total amount of P# + P
"SK" : 0.105* 1000.,      # (nM) total amount of SK* + SK
"HSF" : 300.,             # (muM) reference value
"Genes" : 0.0022* 1000.,  # (nM) total amount of G + GHSF + G*HSF
"mRNA" : 15.,             # (muM) reference value
"HSP" : 10000./ 1000.,    # (mM) reference value                                     
                                   }



