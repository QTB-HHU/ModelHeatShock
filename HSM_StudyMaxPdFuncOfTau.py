
import math

from HSM_SimulateClass import *
from HSM_StudyEquilibrium import *
from HSM_StudyHPproduction import *
from HSM_VaryParamsRMSvsData import *
from HSM_StudyMaxPdFuncOfTau import *
from HSM_calibrationRMSmainFunctions import *

def ComputeMaxUnfoldedProteinsAsFunctionOfTimeToIncreaseTemperature(MyHSM, TauMin, TauMax, NumberOfSteps, SimulationName, FigureName, Tup, Tdown):

    Dtau = (math.log10(TauMax)-math.log10(TauMin)) / (NumberOfSteps)
    ListPd = []
    ListTau = []
    for i in range(NumberOfSteps+1):

        tau = TauMin * math.pow(10,i*Dtau)
        ListTau.append(tau)
        #ListTau.append(TauMin+i*(TauMax-TauMin)/NumberOfSteps)

        TsetPdTau = ParametersSet({"Ttype": 8, "Tin": Tdown, "Tup": Tup, "tau": ListTau[i], "ta": 2. * 60. + vorl})
        timesetPdTau = ParametersSet({"t_start": 0., "t_stop": (vorl + ListTau[i] + 62. * 60.), "delta_t": 5.0})
        SimulationPdTau = Simulate(MyHSM, timesetPdTau, TsetPdTau, SimulationName)
        SimulationPdTau.TimeRun(AvoidPlots="Yes")

        ScalingFactorToAUforPd = SetOfReferenceValuesForPlotting["Ptot"] 

        FirstIndex = round(vorl / timesetPdTau.DefaultParams["delta_t"])

        NewList = []
        for i in range(len(SimulationPdTau.Ph)-FirstIndex):
            NewList.append(SimulationPdTau.Ph[i+FirstIndex][0]/ScalingFactorToAUforPd)

        ListPd.append(max(NewList))
        print(max(NewList))

    ListTauLog = []
    for j in range(len(ListTau)):
        ListTauLog.append(math.log10(ListTau[j]/60.))

    fig = figure()
    plt.plot(ListTauLog, ListPd, color="blue", marker="o")
    plt.xlim(math.log10(TauMin/60.), math.log10(TauMax/60.))
    plt.xlabel(r"$\tau$ (min)", fontsize="large")
    plt.ylabel(r"$P^\#_{MAX}$ (a.u.)", fontsize="large")
    plt.legend(loc="upper right")

    # Set axis' ticks at proper places (convert the axis from linear in the exponents to log)
    plt.tick_params(axis='x', labelsize=14)
    ax=plt.gca()
    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels(["1","10","100","1000"])
    ax.xaxis.set_ticks([math.log10(0.2),math.log10(0.30),math.log10(0.40),math.log10(0.50),math.log10(0.60),math.log10(0.70),math.log10(0.80),
    math.log10(0.90),math.log10(2.0),math.log10(3.0),math.log10(4.0),math.log10(5.0),math.log10(6.0),math.log10(7.0),math.log10(8.0),math.log10(9.0),
    math.log10(20),math.log10(30),math.log10(40),math.log10(50),math.log10(60),math.log10(70),math.log10(80),math.log10(90),
    math.log10(200),math.log10(300),math.log10(400),math.log10(500),math.log10(600),math.log10(700),math.log10(800),math.log10(900)], minor = True)

    plt.show()

    PlotAndSave(fig, FigureName, "PS", 0, 0)


def PlotTemperatureManyTau(Tup, Tdown, FigureExtension):

    """ Plot the temperature T(t) """

    TsetTau1 = ParametersSet({"Ttype": 8, "Tin": Tdown, "Tup": Tup, "tau": 1*60, "ta": 2. * 60.})
    TsetTau10 = ParametersSet({"Ttype": 8, "Tin": Tdown, "Tup": Tup, "tau": 10*60, "ta": 2. * 60.})       
    TsetTau50 = ParametersSet({"Ttype": 8, "Tin": Tdown, "Tup": Tup, "tau": 50*60, "ta": 2. * 60.})         
    TsetTau100 = ParametersSet({"Ttype": 8, "Tin": Tdown, "Tup": Tup, "tau": 100*60, "ta": 2. * 60.})

    tlist = []
    Tlist1 = []
    Tlist10 = []  
    Tlist50 = []      
    Tlist100 = []

    tstart = -50
    tstop = 150
    tcurrent = tstart
    for i in range(100*(tstop-tstart)):
        
        tcurrent = tcurrent + 0.01
        #print(" ")
        #print(tcurrent)
        T1 = T(tcurrent*60., TsetTau1)
        T10 = T(tcurrent*60., TsetTau10)    
        T50 = T(tcurrent*60., TsetTau50)    
        T100 = T(tcurrent*60., TsetTau100)

        tlist.append(tcurrent)
        Tlist1.append(T1)
        Tlist10.append(T10)
        Tlist50.append(T50)
        Tlist100.append(T100)

    fig0 = figure()

    plt.xlabel('Time (min)', fontsize=18)
    plt.ylabel(r'Temperature ($^\circ$C)', fontsize=18)
    #plt.xlim(0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    plt.ylim(24., 45.)
    plt.xlim(-50., 150.)

    plt.plot(tlist, Tlist1, color='b', linestyle=":", linewidth=3., label=r"Temperature, $\tau=1$ min")
    plt.plot(tlist, Tlist10, color='orange', linestyle="-.", linewidth=2.5, label=r"Temperature, $\tau=10$ min")
    plt.plot(tlist, Tlist50, color='g', linestyle="--", linewidth=2., label=r"Temperature, $\tau=50$ min")
    plt.plot(tlist, Tlist100, color='r', linestyle="-", linewidth=1.5, label=r"Temperature, $\tau=100$ min")

    plt.legend(loc="lower right", fontsize="medium", fancybox=True)

    PlotAndSave(fig0, "TemperatureManyTaus"+FigureExtension, "PS", 0, 1)


