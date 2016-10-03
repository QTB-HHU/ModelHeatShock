from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

from HSM_HSModelClass import *
from HSM_PlottingFunctions import *
from HSM_FunctsOfSimulateClass import *


class Simulate:
    """ Creates a simulation (of which one can make time run, plot the results, make simulation of experiments """

    def __init__(self, Model, timeset, TparamSet, NameOfSimulation):

        self.name = NameOfSimulation
        self.model = Model

        self.timeset = timeset
        self.TparamSet = TparamSet

        # Number of time steps, +1 for initial condit.
        self.StepsNum = np.floor(
            (self.timeset.CurrentParams["t_stop"] - self.timeset.CurrentParams["t_start"]) / self.timeset.CurrentParams[
                "delta_t"]) + 1

        # Create vectors to store and plot trajectories
        for attr in ["t", "Tplot", "P", "Ph", "S", "Ss", "F", "Fs", "G", "FsG", "FG", "RF", "RHP", "HP"]:
            setattr(self, attr, np.zeros((self.StepsNum, 1)))
        # Add two vectors for the additional variables used to fit exp of Schroda et al 2000 (not part of original model)
        for attr in ["RHP_ARS", "HP_ARS"]:
            setattr(self, attr, np.zeros((self.StepsNum, 1)))

        # Fill in the initial conditions
        self.t[0] = -vorl / 60.
        self.Tplot[0] = self.TparamSet.CurrentParams["Tin"]
        self.P[0] = self.model.ParamsSetIC.CurrentParams["Pin"] / 1000.
        self.Ph[0] = self.model.ParamsSetIC.CurrentParams["Phin"] / 1000.
        self.S[0] = self.model.ParamsSetIC.CurrentParams["Sin"] * 1000.
        self.Ss[0] = self.model.ParamsSetIC.CurrentParams["Ssin"] * 1000.
        self.F[0] = self.model.ParamsSetIC.CurrentParams["Fin"]
        self.Fs[0] = self.model.ParamsSetIC.CurrentParams["Fsin"]
        self.G[0] = self.model.ParamsSetIC.CurrentParams["Gin"] * 1000.
        self.FsG[0] = self.model.ParamsSetIC.CurrentParams["FsGin"] * 1000.
        self.FG[0] = self.model.ParamsSetIC.CurrentParams["FGin"] * 1000.
        self.RF[0] = self.model.ParamsSetIC.CurrentParams["RFin"]
        self.RHP[0] = self.model.ParamsSetIC.CurrentParams["RHPin"]
        self.HP[0] = self.model.ParamsSetIC.CurrentParams["HPin"] / 1000.

    def PlotTemperature(self, tminMANUAL=("No",123456789)):
        """ Plot the temperature T(t) """
        FuncPlotTemperature(self, tminMANUAL=tminMANUAL)

    def TimeRun(self, ModifyKset="No", ParamsToBeModif={}, TimeOfMod=vorl, DirectControlnuPp=["No",123456789], AvoidPlots="No", ZoomInPanelA="No",tminMANUAL=("No",123456789)):
        """Solve ODEs system for given parameter values, plot time behaviour of T and concentrations"""
        #print("!!! In TimeRun   " + str(DirectControlnuPp))
        FuncTimeRun(self, DirectControlnuPp, AvoidPlots, ModifyKset=ModifyKset, ParamsToBeModif=ParamsToBeModif, TimeOfMod=TimeOfMod, ZoomInPanelA=ZoomInPanelA, tminMANUAL=tminMANUAL)

    def TimeRun10eqs(self, DirectControlnuPp=["No",123456789], AvoidPlots="No"):
        """Solve ODEs system for given parameter values, plot time behaviour of T and concentrations"""
        FuncTimeRun10eqs(self, DirectControlnuPp, AvoidPlots)

    def TimeRun9eqs(self, DirectControlnuPp=["No",123456789], AvoidPlots="No"):
        """Solve ODEs system for given parameter values, plot time behaviour of T and concentrations"""
        FuncTimeRun9eqs(self, DirectControlnuPp, AvoidPlots)

    def FeedingExperimentStaur2Dplots(self, ParamName, ListOf3Values):
        """ Plots mRNAf and Fs for different values of ParamName, reproduces staurosporine plots """
        FuncFeedingExperimentStaur2Dplots(self, ParamName, ListOf3Values)

    def FeedingExperiment3Dplots(self, switchHPmRNAF, ParamName, ParamLabel, ParMin, ParMAX, ParStepsNumber, Camera):
        """ Make a 3D plot of some concentrations as functions of time while variating ParamName between ParMin and ParMAX """
        FuncFeedingExperiment3Dplots(self, switchHPmRNAF, ParamName, ParamLabel, ParMin, ParMAX, ParStepsNumber, Camera)

    def FeedingExperimentPlotsVsData(self, ParamName, ListOfValues,
                                     ModelLegend, YVariabToPlot, DataFileName, DataLegend, Ylabel, NameOfFigure,
                                     LegendPosition, ColumnNumber, MultipleT="No", TParamsToUpdate={}):
        """ Simulate feeding experiment by variating a k between values, and compare with literature data """
        FuncFeedingExperimentPlotsVsData(self, ParamName, ListOfValues,
                                         ModelLegend, YVariabToPlot, DataFileName, DataLegend, Ylabel, NameOfFigure,
                                         LegendPosition, ColumnNumber, MultipleT, TParamsToUpdate)

    def TimeCourseVsDataPlot(self):
        """ Compare time course with literature data """
        FuncTimeCourseVsDataPlot(self)

    def TimeRunPlusARS(self):
        """ ...WORK IN PROGRESS... """
        FuncTimeRunPlusARS(self)

    def TimeRunPlusARSdoubleHS(self, EmptyListToBeFilled, AvoidPlots="No"):
        """ ...WORK IN PROGRESS... """
        FuncTimeRunPlusARSdoubleHS(self, EmptyListToBeFilled, AvoidPlots)





