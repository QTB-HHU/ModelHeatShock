import scipy.optimize as sciopt
import numdifftools as nd
import numpy.linalg as lg
import numpy as np

from HSM_HSModelClass import *
from HSM_PlottingFunctions import *
from HSM_SimulateClass import *


def FinalValuesOf(Simulation):
    """ return final values of the 12 y concentrations of a simulation run up to a certain time """

    LAST = (len(Simulation.t) - 1)

    # (microM) P, Ph, S, Ss, F, Fs, G, FsG, FG, RF, RHP, HP
    yEquilibriumGuessUNFLATTENED = (Simulation.P[LAST] * 1000., Simulation.Ph[LAST] * 1000., Simulation.S[LAST] / 1000.,
                                    Simulation.Ss[LAST] / 1000., Simulation.F[LAST], Simulation.Fs[LAST],
                                    Simulation.G[LAST] / 1000., Simulation.FsG[LAST] / 1000.,
                                    Simulation.FG[LAST] / 1000.,
                                    Simulation.RF[LAST], Simulation.RHP[LAST], Simulation.HP[LAST] * 1000.)

    yEquilibriumGuess = [val for sublist in yEquilibriumGuessUNFLATTENED for val in
                         sublist]  # Remove unnecessary parentheses
    print("Final values (after %r minutes) of the concentrations y(t) of species in the system:" % (
    Simulation.timeset.CurrentParams["t_stop"] / 60.))
    print("P, Ph, S, Ss, F, Fs, G, FsG, FG, RF, RHP, HP (microM)")
    print(yEquilibriumGuess)
    print("This can be taken as an initial guess for an equilibrium point of the system.")
    print("")
    return yEquilibriumGuess


def FinalValuesOf_10eqs(Simulation):
    """ return final values of the 10 y concentrations of a simulation run up to a certain time """

    LAST = (len(Simulation.t) - 1)

    # (microM) Ph, Ss, F, Fs, G, FsG, FG, RF, RHP, HP
    yEquilibriumGuessUNFLATTENED = (Simulation.Ph[LAST] * 1000.,
                                    Simulation.Ss[LAST] / 1000., Simulation.F[LAST], Simulation.Fs[LAST],
                                    Simulation.G[LAST] / 1000., Simulation.FsG[LAST] / 1000.,
                                    Simulation.FG[LAST] / 1000.,
                                    Simulation.RF[LAST], Simulation.RHP[LAST], Simulation.HP[LAST] * 1000.)

    yEquilibriumGuess = [val for sublist in yEquilibriumGuessUNFLATTENED for val in
                         sublist]  # Remove unnecessary parentheses
    print("Final values (after %r minutes) of the concentrations y(t) of species in the system:" % (
    Simulation.timeset.CurrentParams["t_stop"] / 60.))
    print("Ph, Ss, F, Fs, G, FsG, FG, RF, RHP, HP (microM)")
    print(yEquilibriumGuess)
    print("This can be taken as an initial guess for an equilibrium point of the system.")
    print("")
    return yEquilibriumGuess


def FinalValuesOf_9eqs(Simulation):
    """ return final values of the 10 y concentrations of a simulation run up to a certain time """

    LAST = (len(Simulation.t) - 1)

    # (microM) Ph, Ss, F, Fs, FsG, FG, RF, RHP, HP
    yEquilibriumGuessUNFLATTENED = (Simulation.Ph[LAST] * 1000.,
                                    Simulation.Ss[LAST] / 1000., Simulation.F[LAST], Simulation.Fs[LAST], 
                                    Simulation.FsG[LAST] / 1000.,
                                    Simulation.FG[LAST] / 1000.,
                                    Simulation.RF[LAST], Simulation.RHP[LAST], Simulation.HP[LAST] * 1000.)

    yEquilibriumGuess = [val for sublist in yEquilibriumGuessUNFLATTENED for val in
                         sublist]  # Remove unnecessary parentheses
    print("Final values (after %r minutes) of the concentrations y(t) of species in the system:" % (
    Simulation.timeset.CurrentParams["t_stop"] / 60.))
    print("Ph, Ss, F, Fs, FsG, FG, RF, RHP, HP (microM)")
    print(yEquilibriumGuess)
    print("This can be taken as an initial guess for an equilibrium point of the system.")
    print("")
    return yEquilibriumGuess


def FindRootOfFunc(func, yGuess, Arguments=()):
    """ find roots of f """
    # return sciopt.fsolve(func, yGuess, Arguments, full_output=0)
    return sciopt.root(func, yGuess, Arguments)


def FindRootOfFuncAndPrint(func, yGuess, Arguments=()):
    """ find roots of f """

    # An equilibrium corresponds to no motion, i.e. when the right-hand sides of all the ODEs are zero-valued. 
    # So for a system dy/dt = f(y), where y is a vector, you just need to solve f(y) = 0,
    # which means finding the zeros of the vector field f(y). This may have multiple solutions, of course.
    # It's best to have a good initial guess as a start condition for searching, and integration can help find that.

    Radici = FindRootOfFunc(func, yGuess, Arguments)

    print("Let's look for roots of the equations, to find equilibrium points.")
    print("Have I found a root? " + str(Radici.success) + ". " + str(Radici.message))
    print("The root found is:")
    print(Radici.x)
    print("")

    return Radici.x


def ODEsSysthAsFunction(y, ksetDict, REACparamSet, FixedTemperature, DirectControlnuPp, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG):
    TparamSet = ParametersSet({"Ttype": 0, "Tin": FixedTemperature})
    if NumbEquations9or10or12 == 12:
        return f(0., y, ksetDict, TparamSet, REACparamSet, DirectControlnuPp)
    elif NumbEquations9or10or12 == 10:
        return f10eqs(0., y, ksetDict, TparamSet, REACparamSet, DirectControlnuPp, IC_PplusPp, IC_SplusSs)
    elif NumbEquations9or10or12 == 9:
        return f9eqs(0., y, ksetDict, TparamSet, REACparamSet, DirectControlnuPp, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)
    # t enters only via the time dependent temperature, which we set here to be a constant, so we can put t to any value


def ODEsSysthResultInColumn(y, ksetDict, REACparamSet, FixedTemperature, DirectControlnuPp, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG):
    return np.reshape(ODEsSysthAsFunction(y, ksetDict, REACparamSet, FixedTemperature, DirectControlnuPp, 
                                          NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG), (-1, 1))


def EigenvaluesOfJacobianAtEquilibrium(yEquilibriumPoint, ParamSetRATES, ParamSetForREACTIONS, FixedTemp, DirectControlnuPp, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG):
    """ Compute the Jacobian of the vector field associated to our ODEs system and its eigenvalues """

    # Compute the Jacobian of the vector field f(y) associated to our ODEs system 
    # (i.e. the square matrix df(y)/dy, where f and y are vectors of the same dimension)
    JacobOfODEsSysth = nd.Jacobian(ODEsSysthResultInColumn)

    # Evaluate this Jacobian (for given parameters' values) at the y point which should correspond to the equilibrium
    JacobAtEquilibriumPoint = JacobOfODEsSysth(yEquilibriumPoint, ParamSetRATES, ParamSetForREACTIONS, FixedTemp, DirectControlnuPp, 
                                               NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)

    # Compute the eigenvalues of the Jacobian above, to study stability of the equilibrium point
    # All eigenvalues are negative --> STABLE equilibrium
    # Some eigenvalues are positive --> UNSTABLE equilibrium
    # Some eigenvalues are 0 --> requires to be studied with Jacobi method
    # Some eigenvalues are complex --> ...
    eigenvalues, eigenvectors = lg.eig(JacobAtEquilibriumPoint)
    print(
        "The eigenvalues of the Jacobian of the vecotr field associated to our ODEs system, \ncomputed at the equilibrium point provided, are: ")
    print(eigenvalues)
    print("")
    print("")
    return eigenvalues


def UnfoldingRateFuncOfT(kP0p, n1, T0const, T):
    """Computing the rate of unfolding proteins as a function of T"""
    TparamSet = ParametersSet({"Ttype": 0, "Tin": T})
    result = nuPp(1, 0., kP0p, n1, T0const, TparamSet, ["No",123456789.]) # Computed with [P]=1 (i.e. divided by [P])
    return result


### The following set of functions are built by using the ones above and used to do a systematc study of how 
# the equilibrium (steady state) point and the corresponding eigenvalues of the Jacobian (computed at that equilibrium point) evolve. ###

def StudyUnfoldingRateFuncOfT(kP0p, Tmin, Tmax, LegendPosition, figurename, n1List=[1, 5, 10], T0constList=[36., 39.]):
    """Studying how the rate of unfolding proteins as function of T changes by using different parameters"""

    ArraysToPlot=[]
    for n1 in n1List:
        for T0 in T0constList:
            TArray=[]
            Array=[]
            for T in range(Tmin,Tmax+1,1):
                rate = UnfoldingRateFuncOfT(kP0p, n1, T0, T)
                #print(T, n1, T0, rate)
                TArray.append(T)
                Array.append(rate)
            ArraysToPlot.append([r'T$_0$=%r$^\circ$C, n=%r' % (round(T0), n1), Array])
    fig = figure()
    SinglePlot(fig, TArray, ArraysToPlot, r"Temperature ($^\circ$C)", Tmin, Tmax, 
                r"$\frac{\nu'_P}{\left[P\right]}\times\frac{%r s^{-1}}{k'_P}$ $(s^{-1})$" % round(kP0p), 0., 0., LegendPosition, figurename, "large")


def PlotTrajectoriesOfEquilibriumPoint(nupPoverP, ListOfConcentrationsAtEquilibrium, FigureNameEquilibrium, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG):
    """ Plot the evolution of the equilibrium point in 6 separated subplots """

    nupPoverPmin = min(nupPoverP)
    nupPoverPmax = max(nupPoverP)
    ylabel = r"$\nu'_P$/[P] (s$^{-1}$)"

    if NumbEquations9or10or12 == 12:
        P_equilibrium = ListOfConcentrationsAtEquilibrium[0]
        Ph_equilibrium = ListOfConcentrationsAtEquilibrium[1]
        S_equilibrium = ListOfConcentrationsAtEquilibrium[2]
        Ss_equilibrium = ListOfConcentrationsAtEquilibrium[3]
        F_equilibrium = ListOfConcentrationsAtEquilibrium[4]
        Fs_equilibrium = ListOfConcentrationsAtEquilibrium[5]
        G_equilibrium = ListOfConcentrationsAtEquilibrium[6]
        FsG_equilibrium = ListOfConcentrationsAtEquilibrium[7]
        FG_equilibrium = ListOfConcentrationsAtEquilibrium[8]
        RF_equilibrium = ListOfConcentrationsAtEquilibrium[9]
        RHP_equilibrium = ListOfConcentrationsAtEquilibrium[10]
        HP_equilibrium = ListOfConcentrationsAtEquilibrium[11]

    elif NumbEquations9or10or12 == 10:
        Ph_equilibrium = ListOfConcentrationsAtEquilibrium[0]
        Ss_equilibrium = ListOfConcentrationsAtEquilibrium[1]
        F_equilibrium = ListOfConcentrationsAtEquilibrium[2]
        Fs_equilibrium = ListOfConcentrationsAtEquilibrium[3]
        G_equilibrium = ListOfConcentrationsAtEquilibrium[4]
        FsG_equilibrium = ListOfConcentrationsAtEquilibrium[5]
        FG_equilibrium = ListOfConcentrationsAtEquilibrium[6]
        RF_equilibrium = ListOfConcentrationsAtEquilibrium[7]
        RHP_equilibrium = ListOfConcentrationsAtEquilibrium[8]
        HP_equilibrium = ListOfConcentrationsAtEquilibrium[9]

        P_equilibrium = []
        S_equilibrium = []
        for i in range(len(Ph_equilibrium)):
            P_equilibrium.append(IC_PplusPp - Ph_equilibrium[i])
            S_equilibrium.append(IC_SplusSs - Ss_equilibrium[i])

    elif NumbEquations9or10or12 == 9:
        Ph_equilibrium = ListOfConcentrationsAtEquilibrium[0]
        Ss_equilibrium = ListOfConcentrationsAtEquilibrium[1]
        F_equilibrium = ListOfConcentrationsAtEquilibrium[2]
        Fs_equilibrium = ListOfConcentrationsAtEquilibrium[3]
        FsG_equilibrium = ListOfConcentrationsAtEquilibrium[4]
        FG_equilibrium = ListOfConcentrationsAtEquilibrium[5]
        RF_equilibrium = ListOfConcentrationsAtEquilibrium[6]
        RHP_equilibrium = ListOfConcentrationsAtEquilibrium[7]
        HP_equilibrium = ListOfConcentrationsAtEquilibrium[8]

        P_equilibrium = []
        S_equilibrium = []
        G_equilibrium = []
        for i in range(len(Ph_equilibrium)):
            P_equilibrium.append(IC_PplusPp - Ph_equilibrium[i])
            S_equilibrium.append(IC_SplusSs - Ss_equilibrium[i])
            G_equilibrium.append(IC_GplusFsGplusFG - FsG_equilibrium[i] - FG_equilibrium[i])

    else:
        print("Error in PlotTrajectoriesOfEquilibriumPoint")

    fig1 = figure()

    ax1 = plt.subplot(321)
    SubPlot(ax1, nupPoverP, [["P", P_equilibrium], [r"P$^\#$", Ph_equilibrium]],
            ylabel, nupPoverPmin, nupPoverPmax, 'Concentration\nat steady state\n(fraction of total)', 0., 0., "center right", "A", UsePointMarker="Yes")

    ax2 = plt.subplot(322)
    SubPlot(ax2, nupPoverP, [[r"SK", S_equilibrium], [r"SK$^*$", Ss_equilibrium]],
            ylabel, nupPoverPmin, nupPoverPmax, 'Concentration\nat steady state\n(fraction of total)', 0., 0., "center right", "B", UsePointMarker="Yes")

    ax3 = plt.subplot(323)
    SubPlot(ax3, nupPoverP, [[r"HSF", F_equilibrium], [r"HSF$^*$", Fs_equilibrium]],
            ylabel, nupPoverPmin, nupPoverPmax, 'Concentration\nat steady state\n(a.u.)', 0, 0, "center right", "C", UsePointMarker="Yes")

    ax4 = plt.subplot(324)
    SubPlot(ax4, nupPoverP, [[r"G", G_equilibrium], [r"HSF$^*$G", FsG_equilibrium], [r"HSFG", FG_equilibrium]],
            ylabel, nupPoverPmin, nupPoverPmax, 'Concentration\nat steady state\n(fraction of total)', 0, 0, "center right", "D", UsePointMarker="Yes")

    ax5 = plt.subplot(325)
    SubPlot(ax5, nupPoverP, [[r"mR$_{F}$", RF_equilibrium], [r"mR$_{HP}$", RHP_equilibrium]],
            ylabel, nupPoverPmin, nupPoverPmax, 'Concentration\nat steady state\n(a.u.)', 0, 0, "center right", "E", UsePointMarker="Yes")

    ax6 = plt.subplot(326)
    SubPlot(ax6, nupPoverP, [[r"HP", HP_equilibrium]],
            ylabel, nupPoverPmin, nupPoverPmax, 'Concentration\nat steady state\n(a.u.)', 0, 0, "center right", "F", UsePointMarker="Yes")

    PlotAndSave(fig1, FigureNameEquilibrium, "PS", 1, 1)


    
def PlotEvolutionOfEigenvalues(nupPoverP, ListOfEigenvalues, FigureNameEigenvalues, NumbEquations9or10or12):
    """ Plot the evolution of the 12 eigenvalues in 6 separated subplots """

    nupPoverPmin = min(nupPoverP)
    nupPoverPmax = max(nupPoverP)
    ylabel = r"$\nu'_P$/[P] (s$^{-1}$)"

    fig1 = figure()

    if NumbEquations9or10or12 == 12:
        NumberHorizontalPanelsInPlot = 3
        NumberVerticalPanelsInPlot = 4
    elif NumbEquations9or10or12 == 10:
        NumberHorizontalPanelsInPlot = 2
        NumberVerticalPanelsInPlot = 5
    elif NumbEquations9or10or12 == 9:
        NumberHorizontalPanelsInPlot = 3
        NumberVerticalPanelsInPlot = 3
    else:
        print("Error in PlotEvolutionOfEigenvalues")

    ListLegend = [r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$", r"$\lambda_4$", 
                  r"$\lambda_5$", r"$\lambda_6$", r"$\lambda_7$", r"$\lambda_8$", 
                  r"$\lambda_9$", r"$\lambda_10$", r"$\lambda_11$", r"$\lambda_12$"]

    Alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]

    ax = []
    for k in range(NumbEquations9or10or12):
        ax.append(plt.subplot(NumberHorizontalPanelsInPlot, NumberVerticalPanelsInPlot, k+1))
        SubPlot(ax[k], nupPoverP, [[ListLegend[k], ListOfEigenvalues[k]]],
                ylabel, nupPoverPmin, nupPoverPmax, r'Eigenvalue ($\mu$M ???)', 0., 0., "center right", Alphabet[k], UsePointMarker="Yes")

    PlotAndSave(fig1, FigureNameEigenvalues, "PS", 0, 0)


def StudySteadyStateAndEquilibriumEvolution(MyHSM, FigureExtension, DefaultParamSetRATES, DefaultParamSetForREACTIONS, nupPoverPmax, Nsteps, FigureNameEquilibrium, FigureNameEigenvalues, NumbEquations9or10or12 = 12, IC_PplusPp = 0., IC_SplusSs = 0., IC_GplusFsGplusFG = 0.):
    """Study the evolution of the concentrations at equilibrium and the stability of such equilibrium as function of P degeneration rate"""

    # 0) Prepare parameter sets and empty arrays that will be used in the for loop

    Tis0 = 0. # Set T=0°C so that you have nu'P=0 NO MATTER what k'P, [P], T0 and n1 are!!!
    Tset = ParametersSet({"Ttype": 0, "Tin": Tis0})
    # Note that this is used just because an argument is required, but since then DirectControlnuPp=["Yes",nupPoverP]
    # will be specified, the value of Temperature will not be used (instead the value of nupP/P will be changed directly)

    TimeAfterWhichYouTest = vorl
    timeset = ParametersSet({"t_start": 0., "t_stop": TimeAfterWhichYouTest, "delta_t": 5.0})

    nupPoverPvector = []
    yEquilibriumInitialGuess = []
    yEquilibriumPoint = []
    eigenvalues = []       
       
    i = 0
    for nupPoverP in np.linspace(0., nupPoverPmax, Nsteps):  # the loop increases (rate nupP)/(Concentration of P) from 0 to nupPmax,

        # 1) register the value of nupPoverP into a vector
        print()
        print(str(i) + "   The value of nupPoverP is " + str(nupPoverP))
        nupPoverPvector.append(nupPoverP)

        # 2) simulate the system for a very long time specified by timeset
        # and use as a guess for the equilibrium (steady state) point the variables' values at the final time of the simulation above
        Simulation = Simulate(MyHSM, timeset, Tset, "SimulationEvolutionOfSteadyState" + FigureExtension)
        if NumbEquations9or10or12 == 12:
            Simulation.TimeRun(DirectControlnuPp=["Yes",nupPoverP], AvoidPlots="Yes")
            yEquilibriumInitialGuess.append(FinalValuesOf(Simulation))
        elif NumbEquations9or10or12 == 10:
            Simulation.TimeRun10eqs(DirectControlnuPp=["Yes",nupPoverP], AvoidPlots="Yes")
            yEquilibriumInitialGuess.append(FinalValuesOf_10eqs(Simulation))
        elif NumbEquations9or10or12 == 9:
            Simulation.TimeRun9eqs(DirectControlnuPp=["Yes",nupPoverP], AvoidPlots="Yes")
            yEquilibriumInitialGuess.append(FinalValuesOf_9eqs(Simulation))
        else: 
            print("Error in StudySteadyStateAndEquilibriumEvolution")

        print(yEquilibriumInitialGuess[i])

        # 3) Find a root of the system f(y), i.e. a point of equilibrium, starting from the guess for the array of concentrations found above
        yEquilibriumPoint.append(FindRootOfFuncAndPrint(ODEsSysthAsFunction, yEquilibriumInitialGuess[i],
                                                        (DefaultParamSetRATES, DefaultParamSetForREACTIONS, Tis0, ["Yes",nupPoverP], 
                                                         NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG))) 
       
        # 4) Compute the Jacobian of f(y) and find its eigenvalues to determine the stability of the equilibrium point found above
        eigenvalues.append(EigenvaluesOfJacobianAtEquilibrium(yEquilibriumPoint[i], DefaultParamSetRATES, 
                           DefaultParamSetForREACTIONS, Tis0, ["Yes",nupPoverP], NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG))

        i = i + 1  # increase the index used to put the results of each loop into arrays

    # 5) Split the arrays "yEquilibriumPoint" and "eigenvalues" into arrays containing only one concentration or one eigenvalue, for plotting
    ListOfConcentrationsAtEquilibrium = []
    ListOfEigenvalues = []
    for k in range(NumbEquations9or10or12):
        ListOfConcentrationsAtEquilibrium.append([])
        ListOfEigenvalues.append([])

    eigenvalues_flattened = [val for sublist in eigenvalues for val in sublist]  # Flatten the lists of eigenvalues to remove unnecessary parentheses

    for j in range(len(yEquilibriumPoint)):
        for k in range(NumbEquations9or10or12):
            ListOfConcentrationsAtEquilibrium[k].append(yEquilibriumPoint[j][k])
            ListOfEigenvalues[k].append(eigenvalues_flattened[j][k])
  
    # Print everything on the terminal just to check that it looks OK
    print()
    print(nupPoverPvector)
    print()
    print(yEquilibriumInitialGuess)
    print()
    print(yEquilibriumPoint)
    print()
    print(eigenvalues)
    print()
    print(eigenvalues_flattened)
    print()
    for k in range(NumbEquations9or10or12):
        print(k)
        print(ListOfEigenvalues[k])
    print()

    # 6) Plot the evolution of the EQUILIBRIUM point as a function of nupP/[P]  
    PlotTrajectoriesOfEquilibriumPoint(nupPoverPvector, ListOfConcentrationsAtEquilibrium, FigureNameEquilibrium, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)
    
    # 7) Plot the evolution of the EIGENVALUES as a function of nupP/[P]
    PlotEvolutionOfEigenvalues(nupPoverPvector, ListOfEigenvalues, FigureNameEigenvalues, NumbEquations9or10or12)



def PlotTrajectoriesOfEquilibriumPointFuncOfTEMPERATURE(VectorOfTEMPERATURES, ListOfConcentrationsAtEquilibrium, FigureNameEquilibrium, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG):
    """ Plot the evolution of the equilibrium point in 6 separated subplots """

    TEMPERAUREmin = min(VectorOfTEMPERATURES)
    TEMPERAUREmax = max(VectorOfTEMPERATURES)
    xlabel = r"T (°C)"

    if NumbEquations9or10or12 == 12:
        P_equilibrium = ListOfConcentrationsAtEquilibrium[0]
        Ph_equilibrium = ListOfConcentrationsAtEquilibrium[1]
        S_equilibrium = ListOfConcentrationsAtEquilibrium[2]
        Ss_equilibrium = ListOfConcentrationsAtEquilibrium[3]
        F_equilibrium = ListOfConcentrationsAtEquilibrium[4]
        Fs_equilibrium = ListOfConcentrationsAtEquilibrium[5]
        G_equilibrium = ListOfConcentrationsAtEquilibrium[6]
        FsG_equilibrium = ListOfConcentrationsAtEquilibrium[7]
        FG_equilibrium = ListOfConcentrationsAtEquilibrium[8]
        RF_equilibrium = ListOfConcentrationsAtEquilibrium[9]
        RHP_equilibrium = ListOfConcentrationsAtEquilibrium[10]
        HP_equilibrium = ListOfConcentrationsAtEquilibrium[11]

    elif NumbEquations9or10or12 == 10:
        Ph_equilibrium = ListOfConcentrationsAtEquilibrium[0]
        Ss_equilibrium = ListOfConcentrationsAtEquilibrium[1]
        F_equilibrium = ListOfConcentrationsAtEquilibrium[2]
        Fs_equilibrium = ListOfConcentrationsAtEquilibrium[3]
        G_equilibrium = ListOfConcentrationsAtEquilibrium[4]
        FsG_equilibrium = ListOfConcentrationsAtEquilibrium[5]
        FG_equilibrium = ListOfConcentrationsAtEquilibrium[6]
        RF_equilibrium = ListOfConcentrationsAtEquilibrium[7]
        RHP_equilibrium = ListOfConcentrationsAtEquilibrium[8]
        HP_equilibrium = ListOfConcentrationsAtEquilibrium[9]

        P_equilibrium = []
        S_equilibrium = []
        for i in range(len(Ph_equilibrium)):
            P_equilibrium.append(IC_PplusPp - Ph_equilibrium[i])
            S_equilibrium.append(IC_SplusSs - Ss_equilibrium[i])

    elif NumbEquations9or10or12 == 9:
        Ph_equilibrium = ListOfConcentrationsAtEquilibrium[0]
        Ss_equilibrium = ListOfConcentrationsAtEquilibrium[1]
        F_equilibrium = ListOfConcentrationsAtEquilibrium[2]
        Fs_equilibrium = ListOfConcentrationsAtEquilibrium[3]
        FsG_equilibrium = ListOfConcentrationsAtEquilibrium[4]
        FG_equilibrium = ListOfConcentrationsAtEquilibrium[5]
        RF_equilibrium = ListOfConcentrationsAtEquilibrium[6]
        RHP_equilibrium = ListOfConcentrationsAtEquilibrium[7]
        HP_equilibrium = ListOfConcentrationsAtEquilibrium[8]

        P_equilibrium = []
        S_equilibrium = []
        G_equilibrium = []
        for i in range(len(Ph_equilibrium)):
            P_equilibrium.append(IC_PplusPp - Ph_equilibrium[i])
            S_equilibrium.append(IC_SplusSs - Ss_equilibrium[i])
            G_equilibrium.append(IC_GplusFsGplusFG - FsG_equilibrium[i] - FG_equilibrium[i])

    else:
        print("Error in PlotTrajectoriesOfEquilibriumPoint")

    # RESCALING FOR PLOTTING
    P_equilibrium_rescaled = []
    Ph_equilibrium_rescaled = []
    S_equilibrium_rescaled = []
    Ss_equilibrium_rescaled = []
    F_equilibrium_rescaled = []
    Fs_equilibrium_rescaled = []
    G_equilibrium_rescaled = []
    FsG_equilibrium_rescaled = []
    FG_equilibrium_rescaled = []
    RF_equilibrium_rescaled = []
    RHP_equilibrium_rescaled = []
    HP_equilibrium_rescaled = []
    i = 0
    for element in P_equilibrium:
        P_equilibrium_rescaled.append(P_equilibrium[i]/(1000.*SetOfReferenceValuesForPlotting["Ptot"]))
        Ph_equilibrium_rescaled.append(Ph_equilibrium[i]/(1000.*SetOfReferenceValuesForPlotting["Ptot"]))
        S_equilibrium_rescaled.append(1000.*S_equilibrium[i]/SetOfReferenceValuesForPlotting["SK"])
        Ss_equilibrium_rescaled.append(1000.*Ss_equilibrium[i]/SetOfReferenceValuesForPlotting["SK"])
        F_equilibrium_rescaled.append(F_equilibrium[i]/SetOfReferenceValuesForPlotting["HSF"])
        Fs_equilibrium_rescaled.append(Fs_equilibrium[i]/SetOfReferenceValuesForPlotting["HSF"])
        G_equilibrium_rescaled.append(1000.*G_equilibrium[i]/SetOfReferenceValuesForPlotting["Genes"])
        FsG_equilibrium_rescaled.append(1000.*FsG_equilibrium[i]/SetOfReferenceValuesForPlotting["Genes"])
        FG_equilibrium_rescaled.append(1000.*FG_equilibrium[i]/SetOfReferenceValuesForPlotting["Genes"])
        RF_equilibrium_rescaled.append(RF_equilibrium[i]/SetOfReferenceValuesForPlotting["mRNA"])
        RHP_equilibrium_rescaled.append(RHP_equilibrium[i]/SetOfReferenceValuesForPlotting["mRNA"])
        HP_equilibrium_rescaled.append(HP_equilibrium[i]/(1000.*SetOfReferenceValuesForPlotting["HSP"]))
        i=i+1

    fig1 = figure()

    ax1 = plt.subplot(321)
    #SubPlot(ax1, VectorOfTEMPERATURES, [["P", P_equilibrium], [r"P$^\#$", Ph_equilibrium]],
    #        xlabel, TEMPERAUREmin, TEMPERAUREmax, r'Protein ($\mu$M)', 0., 0., "center left", "A", UsePointMarker="Yes")
    SubPlot(ax1, VectorOfTEMPERATURES, [["P", P_equilibrium_rescaled], [r"P$^\#$", Ph_equilibrium_rescaled]],
            xlabel, TEMPERAUREmin, TEMPERAUREmax, 'Concentration\nat steady state\n(fraction of total)', 0., 0., "center left", "A", UsePointMarker="Yes")

    ax2 = plt.subplot(322)
    SubPlot(ax2, VectorOfTEMPERATURES, [[r"SK", S_equilibrium_rescaled], [r"SK$^*$", Ss_equilibrium_rescaled]],
            xlabel, TEMPERAUREmin, TEMPERAUREmax, 'Concentration\nat steady state\n(fraction of total)', 0., 0., "center right", "B", UsePointMarker="Yes")

    ax3 = plt.subplot(323)
    SubPlot(ax3, VectorOfTEMPERATURES, [[r"HSF", F_equilibrium_rescaled], [r"HSF$^*$", Fs_equilibrium_rescaled]],
            xlabel, TEMPERAUREmin, TEMPERAUREmax, 'Concentration\nat steady state\n(a.u.)', 0, 0, "center left", "C", UsePointMarker="Yes")

    ax4 = plt.subplot(324)
    SubPlot(ax4, VectorOfTEMPERATURES, [[r"G", G_equilibrium_rescaled], [r"HSF$^*$G", FsG_equilibrium_rescaled], [r"HSFG", FG_equilibrium_rescaled]],
            xlabel, TEMPERAUREmin, TEMPERAUREmax, 'Concentration\nat steady state\n(fraction of total)', 0, 0, "center right", "D", UsePointMarker="Yes")

    ax5 = plt.subplot(325)
    SubPlot(ax5, VectorOfTEMPERATURES, [[r"mR$_{F}$", RF_equilibrium_rescaled], [r"mR$_{HP}$", RHP_equilibrium_rescaled]],
            xlabel, TEMPERAUREmin, TEMPERAUREmax, 'Concentration\nat steady state\n(a.u.)', 0, 0, "center left", "E", UsePointMarker="Yes")

    ax6 = plt.subplot(326)
    SubPlot(ax6, VectorOfTEMPERATURES, [[r"HP", HP_equilibrium_rescaled]],
            xlabel, TEMPERAUREmin, TEMPERAUREmax, 'Concentration\nat steady state\n(a.u.)', 0, 0, "center right", "F", UsePointMarker="Yes")

    PlotAndSave(fig1, FigureNameEquilibrium, "PS", 1, 1)


    
def PlotEvolutionOfEigenvaluesFuncOfTEMPERATURE(VectorOfTEMPERAUTRES, ListOfEigenvalues, FigureNameEigenvalues, NumbEquations9or10or12):
    """ Plot the evolution of the 12 eigenvalues in 6 separated subplots """

    TEMPERATUREmin = min(VectorOfTEMPERAUTRES)
    TEMPERATUREmax = max(VectorOfTEMPERAUTRES)
    xlabel = r"T (°C)"

    fig1 = figure()

    if NumbEquations9or10or12 == 12:
        NumberHorizontalPanelsInPlot = 3
        NumberVerticalPanelsInPlot = 4
    elif NumbEquations9or10or12 == 10:
        NumberHorizontalPanelsInPlot = 2
        NumberVerticalPanelsInPlot = 5
    elif NumbEquations9or10or12 == 9:
        NumberHorizontalPanelsInPlot = 3
        NumberVerticalPanelsInPlot = 3
    else:
        print("Error in PlotEvolutionOfEigenvalues")

    ListLegend = [r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$", r"$\lambda_4$", 
                  r"$\lambda_5$", r"$\lambda_6$", r"$\lambda_7$", r"$\lambda_8$", 
                  r"$\lambda_9$", r"$\lambda_10$", r"$\lambda_11$", r"$\lambda_12$"]

    Alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]

    ax = []
    for k in range(NumbEquations9or10or12):
        ax.append(plt.subplot(NumberHorizontalPanelsInPlot, NumberVerticalPanelsInPlot, k+1))
        SubPlot(ax[k], VectorOfTEMPERAUTRES, [[ListLegend[k], ListOfEigenvalues[k]]],
                xlabel, TEMPERATUREmin, TEMPERATUREmax, r'Eigenvalue ($\mu$M ???)', 0., 0., "center right", Alphabet[k], UsePointMarker="Yes")

    PlotAndSave(fig1, FigureNameEigenvalues, "PS", 0, 0)



def StudySteadyStateAndEquilibriumEvolutionFuncOfTEMPERATURE(MyHSM, FigureExtension, DefaultParamSetRATES, DefaultParamSetForREACTIONS, TEMPERATUREmin, TEMPERATUREmax, Nsteps, FigureNameEquilibrium1, FigureNameEquilibrium2, FigureNameEigenvalues1, FigureNameEigenvalues2, NumbEquations9or10or12 = 12, IC_PplusPp = 0., IC_SplusSs = 0., IC_GplusFsGplusFG = 0.):
    """Study the evolution of the concentrations at equilibrium and the stability of such equilibrium as function of P degeneration rate"""

    # 0) Prepare parameter sets and empty arrays that will be used in the for loop

    Tis0 = 0. # Set T=0°C so that you have nu'P=0 NO MATTER what k'P, [P], T0 and n1 are!!!
    Tset = ParametersSet({"Ttype": 0, "Tin": Tis0})
    # Note that this is used just because an argument is required, but since then DirectControlnuPp=["Yes",nupPoverP]
    # will be specified, the value of Temperature will not be used (instead the value of nupP/P will be changed directly)

    TimeAfterWhichYouTest = vorl
    timeset = ParametersSet({"t_start": 0., "t_stop": TimeAfterWhichYouTest, "delta_t": 5.0})

    # XXX-1) GENERATE T INTERVAL
    Tmin = TEMPERATUREmin
    Tmax = TEMPERATUREmax
    VectorOfTEMPERATURES = []
    VectorOfNuPoverP = []
    kP0pHERE = DefaultParamSetRATES["kP0p"]
    # XXX-2) CONVER IT INTO NUP
    for TEMPERATURE in np.linspace(Tmin,Tmax,Nsteps):
        Tset2 = ParametersSet({"Ttype": 0, "Tin": TEMPERATURE})
        NuPp = nuPp(1., 0., kP0pHERE, 123456789, 123456789, Tset2, ["No",123456789]) # [P] = 1. means we are using NuPp/[P] and not NuPp
        VectorOfTEMPERATURES.append(TEMPERATURE)
        VectorOfNuPoverP.append(NuPp)
    #print("")
    #print("LOOK HERE!!!!!!!!!!!!!")
    #print("")
    #print(VectorOfTEMPERATURES)
    #print("")
    #print(VectorOfNuPoverP)
    #print("")

    # XXX-3) COMPUTE STST
    nupPoverPvector = []
    yEquilibriumInitialGuess = []
    yEquilibriumPoint = []
    eigenvalues = []       
       
    i = 0
    for nupPoverP in VectorOfNuPoverP:  # the loop increases (rate nupP)/(Concentration of P) from 0 to nupPmax,

        # 1) register the value of nupPoverP into a vector
        print()
        print(str(i) + "   The value of nupPoverP is " + str(nupPoverP))
        #nupPoverPvector.append(nupPoverP)

        # 2) simulate the system for a very long time specified by timeset
        # and use as a guess for the equilibrium (steady state) point the variables' values at the final time of the simulation above
        Simulation = Simulate(MyHSM, timeset, Tset, "SimulationEvolutionOfSteadyState" + FigureExtension)
        if NumbEquations9or10or12 == 12:
            Simulation.TimeRun(DirectControlnuPp=["Yes",nupPoverP], AvoidPlots="Yes")
            yEquilibriumInitialGuess.append(FinalValuesOf(Simulation))
        elif NumbEquations9or10or12 == 10:
            Simulation.TimeRun10eqs(DirectControlnuPp=["Yes",nupPoverP], AvoidPlots="Yes")
            yEquilibriumInitialGuess.append(FinalValuesOf_10eqs(Simulation))
        elif NumbEquations9or10or12 == 9:
            Simulation.TimeRun9eqs(DirectControlnuPp=["Yes",nupPoverP], AvoidPlots="Yes")
            yEquilibriumInitialGuess.append(FinalValuesOf_9eqs(Simulation))
        else: 
            print("Error in StudySteadyStateAndEquilibriumEvolution")

        print(yEquilibriumInitialGuess[i])

        # 3) Find a root of the system f(y), i.e. a point of equilibrium, starting from the guess for the array of concentrations found above
        yEquilibriumPoint.append(FindRootOfFuncAndPrint(ODEsSysthAsFunction, yEquilibriumInitialGuess[i],
                                                        (DefaultParamSetRATES, DefaultParamSetForREACTIONS, Tis0, ["Yes",nupPoverP], 
                                                         NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG))) 
       
        # 4) Compute the Jacobian of f(y) and find its eigenvalues to determine the stability of the equilibrium point found above
        eigenvalues.append(EigenvaluesOfJacobianAtEquilibrium(yEquilibriumPoint[i], DefaultParamSetRATES, 
                           DefaultParamSetForREACTIONS, Tis0, ["Yes",nupPoverP], NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG))

        i = i + 1  # increase the index used to put the results of each loop into arrays

    # 5) Split the arrays "yEquilibriumPoint" and "eigenvalues" into arrays containing only one concentration or one eigenvalue, for plotting
    ListOfConcentrationsAtEquilibrium = []
    ListOfEigenvalues = []
    for k in range(NumbEquations9or10or12):
        ListOfConcentrationsAtEquilibrium.append([])
        ListOfEigenvalues.append([])

    eigenvalues_flattened = [val for sublist in eigenvalues for val in sublist]  # Flatten the lists of eigenvalues to remove unnecessary parentheses

    for j in range(len(yEquilibriumPoint)):
        for k in range(NumbEquations9or10or12):
            ListOfConcentrationsAtEquilibrium[k].append(yEquilibriumPoint[j][k])
            ListOfEigenvalues[k].append(eigenvalues_flattened[j][k])
  
    # Print everything on the terminal just to check that it looks OK
    print()
    print(VectorOfNuPoverP)
    print()
    print(yEquilibriumInitialGuess)
    print()
    print(yEquilibriumPoint)
    print()
    print(eigenvalues)
    print()
    print(eigenvalues_flattened)
    print()
    for k in range(NumbEquations9or10or12):
        print(k)
        print(ListOfEigenvalues[k])
    print()

    # 6) Plot the evolution of the EQUILIBRIUM point as a function of nupP/[P]  
    PlotTrajectoriesOfEquilibriumPoint(VectorOfNuPoverP, ListOfConcentrationsAtEquilibrium, FigureNameEquilibrium1, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)
    
    # 7) Plot the evolution of the EIGENVALUES as a function of nupP/[P]
    PlotEvolutionOfEigenvalues(VectorOfNuPoverP, ListOfEigenvalues, FigureNameEigenvalues1, NumbEquations9or10or12)

    # XXX-4) PLOT VECTOR T AND VECTOR STST
    # 8) Plot the evolution of the EQUILIBRIUM point as a function of TEMPERATURE  
    PlotTrajectoriesOfEquilibriumPointFuncOfTEMPERATURE(VectorOfTEMPERATURES, ListOfConcentrationsAtEquilibrium, FigureNameEquilibrium2, NumbEquations9or10or12, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)
    
    # 9) Plot the evolution of the EIGENVALUES as a function of TEMPERATURE
    PlotEvolutionOfEigenvaluesFuncOfTEMPERATURE(VectorOfTEMPERATURES, ListOfEigenvalues, FigureNameEigenvalues2, NumbEquations9or10or12)

