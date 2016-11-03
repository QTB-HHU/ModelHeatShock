from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

from HSM_ODEsSystem10or9eqs import *
from HSM_HSModelClass import *
from HSM_PlottingFunctions import *

import time


def SolveODEsystem(self, DirectControlnuPp, ModifyKset="No", ParamsToBeModif={}, TimeOfMod=vorl, Add3rdmRNA_ARS="No"):

    """ Set the system of Ordinary Differential Equations and solve it """

    #start_time = time.time()

    #print("!!! In SolveODEsystem   " + str(DirectControlnuPp))
    # Put the initial conditions into an array (keep the right ordering!)
    y0 = np.zeros((12, 1))
    if Add3rdmRNA_ARS == "Yes":
        y0 = np.zeros((14, 1))

    i = 0
    for key in ["Pin", "Phin", "Sin", "Ssin", "Fin", "Fsin", "Gin", "FsGin", "FGin", "RFin", "RHPin", "HPin"]:
        y0[i] = self.model.ParamsSetIC.CurrentParams[key]
        i = i + 1
    if Add3rdmRNA_ARS == "Yes":
        y0[12] = 0.  # 0.00360 # (microM) mRNA_HP_ARS, mRNA for ARS enzime of which gene has been attached to HSP promot
        y0[13] = 0.  # 1.      # (microM) ARS enzime of which the gene has been attached to the HSP promoter

    # Copy the set of k (rates) parameters and that off params for reactions
    kset = self.model.ParamsSetRATES.CurrentParams.copy()
    ksetMOD = self.model.ParamsSetRATES.CurrentParams.copy()
    ParSetREAC = self.model.ParamsSetForREAC.CurrentParams.copy()

    # Specify integrator ('lsoda' choses automatically between 'bdf' and 'Adams'):
    ODEs = ode(f).set_integrator('lsoda', nsteps=1000, max_step = 60.0) # max_step is in seconds.
    # IMPORTANT: "nsteps" tells the maximum number of steps allowed. if lower the integrator might crash.
    # IMPORTANT: "max_step" tells the maximal lenght of a step in the integration. If not specified or too big, 
    # the integrator can "miss" an heat shock if it is too short (like few minutes) compared to vorl.
    # Thus, max_step should be of same order (or smaller, which is even better!) than the minimal duration of an HS.

    # Set initial conditions for the integrating variable and time
    ODEs.set_initial_value(y0, self.timeset.CurrentParams["t_start"]).set_f_params(kset, self.TparamSet, ParSetREAC,
                                                                                   DirectControlnuPp, Add3rdmRNA_ARS)

    # Integrate the ODEs across each delta_t timestep (i.e. compute y(t) for every t = n * delta_t)
    i = 1
    while ODEs.successful() and i < self.StepsNum:

        if ModifyKset == "Yes" and ODEs.t + self.timeset.CurrentParams["delta_t"] >= TimeOfMod:
            for key in ParamsToBeModif.keys():
                ksetMOD[key] = ParamsToBeModif[key]
            ODEs.set_initial_value(ODEs.y, ODEs.t).set_f_params(ksetMOD, self.TparamSet, ParSetREAC, DirectControlnuPp, Add3rdmRNA_ARS)
            ModifyKset = "No"

        ODEs.integrate(ODEs.t + self.timeset.CurrentParams["delta_t"]) # , step=True

        # Store the results to plot later
        self.t[i] = (ODEs.t - vorl) / 60.  # (min)
        self.Tplot[i] = T(i * self.timeset.CurrentParams["delta_t"], self.TparamSet)  # (°C)
        self.P[i] = ODEs.y[0] / 1000.      # (mM)
        self.Ph[i] = ODEs.y[1] / 1000.     # (mM)
        self.S[i] = ODEs.y[2] * 1000.      # (nM)
        self.Ss[i] = ODEs.y[3] * 1000.     # (nM)
        self.F[i] = ODEs.y[4]              # (microM)
        self.Fs[i] = ODEs.y[5]             # (microM)
        self.G[i] = ODEs.y[6] * 1000.      # (nM)
        self.FsG[i] = ODEs.y[7] * 1000.    # (nM)
        self.FG[i] = ODEs.y[8] * 1000.     # (nM)
        self.RF[i] = ODEs.y[9]             # (microM)
        self.RHP[i] = ODEs.y[10]           # (microM)
        self.HP[i] = ODEs.y[11] / 1000.    # (mM)
        if Add3rdmRNA_ARS == "Yes":
            self.RHP_ARS[i] = ODEs.y[12]   # (microM)
            self.HP_ARS[i] = ODEs.y[13]    # (microM)

        i += 1

    # end_time = time.time()
    # print("Elapsed time was %g seconds" % (end_time - start_time))


def SolveODEsystem10eqs(self, DirectControlnuPp):

    """ Set the system of Ordinary Differential Equations, with only 10 ODEs, and solve it """

    # Put the initial conditions into an array (keep the right ordering!)
    y0 = np.zeros((10, 1))

    i = 0
    for key in ["Phin", "Ssin", "Fin", "Fsin", "Gin", "FsGin", "FGin", "RFin", "RHPin", "HPin"]:
        y0[i] = self.model.ParamsSetIC.CurrentParams[key]
        i = i + 1

    # Copy the set of k (rates) parameters and that off params for reactions
    kset = self.model.ParamsSetRATES.CurrentParams.copy()
    ksetMOD = self.model.ParamsSetRATES.CurrentParams.copy()
    ParSetREAC = self.model.ParamsSetForREAC.CurrentParams.copy()

    # Specify integrator ('lsoda' choses automatically between 'bdf' and 'Adams'):
    ODEs = ode(f10eqs).set_integrator('lsoda', nsteps=1000, max_step = 60.0) # max_step is in seconds.

    IC_PplusPp = self.model.ParamsSetIC.CurrentParams["Pin"] + self.model.ParamsSetIC.CurrentParams["Phin"]
    IC_SplusSs = self.model.ParamsSetIC.CurrentParams["Sin"] + self.model.ParamsSetIC.CurrentParams["Ssin"]
    # Set initial conditions for the integrating variable and time
    ODEs.set_initial_value(y0, self.timeset.CurrentParams["t_start"]).set_f_params(kset, self.TparamSet, ParSetREAC,
                                                                                   DirectControlnuPp, IC_PplusPp, IC_SplusSs)

    # Integrate the ODEs across each delta_t timestep (i.e. compute y(t) for every t = n * delta_t)
    i = 1
    while ODEs.successful() and i < self.StepsNum:

        ODEs.integrate(ODEs.t + self.timeset.CurrentParams["delta_t"]) 

        # Store the results to plot later
        self.t[i] = (ODEs.t - vorl) / 60.  # (min)
        self.Tplot[i] = T(i * self.timeset.CurrentParams["delta_t"], self.TparamSet)  # (°C)

        self.Ph[i] = ODEs.y[0] / 1000.     # (mM)
        self.Ss[i] = ODEs.y[1] * 1000.     # (nM)
        self.F[i] = ODEs.y[2]              # (microM)
        self.Fs[i] = ODEs.y[3]             # (microM)
        self.G[i] = ODEs.y[4] * 1000.      # (nM)
        self.FsG[i] = ODEs.y[5] * 1000.    # (nM)
        self.FG[i] = ODEs.y[6] * 1000.     # (nM)
        self.RF[i] = ODEs.y[7]             # (microM)
        self.RHP[i] = ODEs.y[8]           # (microM)
        self.HP[i] = ODEs.y[9] / 1000.    # (mM)

        self.P[i] = ( self.P[0] + self.Ph[0] ) - self.Ph[i]  # (mM)
        self.S[i] = ( self.S[0] + self.Ss[0] ) - self.Ss[i]  # (nM)

        i += 1


def SolveODEsystem9eqs(self, DirectControlnuPp):

    """ Set the system of Ordinary Differential Equations, with only 10 ODEs, and solve it """

    # Put the initial conditions into an array (keep the right ordering!)
    y0 = np.zeros((9, 1))

    i = 0
    for key in ["Phin", "Ssin", "Fin", "Fsin", "FsGin", "FGin", "RFin", "RHPin", "HPin"]:
        y0[i] = self.model.ParamsSetIC.CurrentParams[key]
        i = i + 1

    # Copy the set of k (rates) parameters and that off params for reactions
    kset = self.model.ParamsSetRATES.CurrentParams.copy()
    ksetMOD = self.model.ParamsSetRATES.CurrentParams.copy()
    ParSetREAC = self.model.ParamsSetForREAC.CurrentParams.copy()

    # Specify integrator ('lsoda' choses automatically between 'bdf' and 'Adams'):
    ODEs = ode(f9eqs).set_integrator('lsoda', nsteps=1000, max_step = 60.0) # max_step is in seconds.

    IC_PplusPp = self.model.ParamsSetIC.CurrentParams["Pin"] + self.model.ParamsSetIC.CurrentParams["Phin"]
    IC_SplusSs = self.model.ParamsSetIC.CurrentParams["Sin"] + self.model.ParamsSetIC.CurrentParams["Ssin"]
    IC_GplusFsGplusFG = self.model.ParamsSetIC.CurrentParams["Gin"] + self.model.ParamsSetIC.CurrentParams["FsGin"] + self.model.ParamsSetIC.CurrentParams["FGin"] 
    # Set initial conditions for the integrating variable and time
    ODEs.set_initial_value(y0, self.timeset.CurrentParams["t_start"]).set_f_params(kset, self.TparamSet, ParSetREAC,
                                                                                   DirectControlnuPp, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG)

    # Integrate the ODEs across each delta_t timestep (i.e. compute y(t) for every t = n * delta_t)
    i = 1
    while ODEs.successful() and i < self.StepsNum:

        ODEs.integrate(ODEs.t + self.timeset.CurrentParams["delta_t"]) 

        # Store the results to plot later
        self.t[i] = (ODEs.t - vorl) / 60.  # (min)
        self.Tplot[i] = T(i * self.timeset.CurrentParams["delta_t"], self.TparamSet)  # (°C)

        self.Ph[i] = ODEs.y[0] / 1000.     # (mM)
        self.Ss[i] = ODEs.y[1] * 1000.     # (nM)
        self.F[i] = ODEs.y[2]              # (microM)
        self.Fs[i] = ODEs.y[3]             # (microM)
        self.FsG[i] = ODEs.y[4] * 1000.    # (nM)
        self.FG[i] = ODEs.y[5] * 1000.     # (nM)
        self.RF[i] = ODEs.y[6]             # (microM)
        self.RHP[i] = ODEs.y[7]           # (microM)
        self.HP[i] = ODEs.y[8] / 1000.    # (mM)

        self.P[i] = ( self.P[0] + self.Ph[0] ) - self.Ph[i]  # (mM)
        self.S[i] = ( self.S[0] + self.Ss[0] ) - self.Ss[i]  # (nM)
        self.G[i] = ( self.G[0] + self.FsG[0] + self.FG[0] ) - self.FsG[i] - self.FG[i] # (nM)

        i += 1


def FuncPlotTemperature(self, tminMANUAL=("No",0.)):
    """ Plot the temperature T(t) """

    fig0 = figure()

    plt.xlabel('Time (min)', fontsize=18)
    if self.TparamSet.CurrentParams["Ttype"] == 6 or (
            self.TparamSet.CurrentParams["Ttype"] == 2 and self.TparamSet.CurrentParams["tb"] >= 12 * 60. * 60.):
        plt.xticks([0., 6. * 60., 12. * 60., 18. * 60., 24 * 60., (24. + 8.) * 60.],
                   ["0", "6", "12", "18", "24", "24+8"])
        plt.xlabel('Time (h)', fontsize="small")
    #plt.xlim(0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    #if tminMANUAL[0] == "Yes":
        #print("Sbleught")
    #print(str(tminMANUAL[1]/60.))
    #print(str((self.timeset.CurrentParams["t_stop"] - vorl) / 60.))
    plt.xlim(0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    plt.ylim(16., 44.)
    plt.ylabel(r'Temperature ($^\circ$C)', fontsize=18)
    plt.plot(self.t, self.Tplot, 'g', linewidth=1.)

    if self.TparamSet.CurrentParams["Ttype"] == 7:
        plt.xticks([0., 3. * 60., 6. * 60., 9. * 60., 12. * 60., 15. * 60., 18. * 60., 21 * 60., 24 * 60.],
                   ["3:00", "6:00", "9:00", "12:00", "15:00", "18:00", "21:00", "24:00", "3:00"])
        plt.xlabel('Time (h)', fontsize="small")
    if self.TparamSet.CurrentParams["Ttype"] == 7 and tminMANUAL[0] == "Yes":
        plt.xticks([0.+tminMANUAL[1]/60., 3. * 60.+tminMANUAL[1]/60., 6. * 60.+tminMANUAL[1]/60., 9. * 60.+tminMANUAL[1]/60., 12. * 60.+tminMANUAL[1]/60., 15. * 60.+tminMANUAL[1]/60., 18. * 60.+tminMANUAL[1]/60., 21 * 60.+tminMANUAL[1]/60., 24 * 60.+tminMANUAL[1]/60.],
                   ["3:00", "6:00", "9:00", "12:00", "15:00", "18:00", "21:00", "24:00", "3:00"])
        plt.xlim(tminMANUAL[1]/60., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
        plt.xlabel('Time (h)', fontsize="small")
    #plt.xlim(0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    plt.ylim(16., 44.)
    plt.ylabel(r'Temperature ($^\circ$C)', fontsize="small")
    plt.plot(self.t, self.Tplot, 'g', linewidth=1.)

    PlotAndSave(fig0, "Temperature_" + self.name, "PS", 0, 1)


def PlotTrajectories(self, ModifyKset="No", ParamsToBeModif={}, TimeOfMod=vorl, ZoomInPanelA="No", tminMANUAL=("No",123456789)):
    """ Plot the trajectories in 6 separated subplots """

    tmax = (self.timeset.CurrentParams["t_stop"] - vorl) / 60.
    tmin = 0.
    if tminMANUAL[0] == "Yes":
        tmin = 0. + tminMANUAL[1]
    #print("\ntminMANUAL[0] = " + str(tminMANUAL[0]) + "   tmin = " + str(tmin) + "\n")
    indexoft0 = int(vorl/self.timeset.CurrentParams["delta_t"])

    fig1 = figure()

    # RESCALING FOR PLOTTING
    Prescaled = []
    Phrescaled = []
    Srescaled = []
    Ssrescaled = []
    Frescaled = []
    Fsrescaled = []
    Grescaled = []
    FsGrescaled = []
    FGrescaled = []
    RFrescaled = []
    RHPrescaled = []
    HPrescaled = []
    i = 0
    for element in self.t:
        Prescaled.append(self.P[i]/SetOfReferenceValuesForPlotting["Ptot"])
        Phrescaled.append(self.Ph[i]/SetOfReferenceValuesForPlotting["Ptot"])
        Srescaled.append(self.S[i]/SetOfReferenceValuesForPlotting["SK"])
        Ssrescaled.append(self.Ss[i]/SetOfReferenceValuesForPlotting["SK"])
        Frescaled.append(self.F[i]/SetOfReferenceValuesForPlotting["HSF"])
        Fsrescaled.append(self.Fs[i]/SetOfReferenceValuesForPlotting["HSF"])
        Grescaled.append(self.G[i]/SetOfReferenceValuesForPlotting["Genes"])
        FsGrescaled.append(self.FsG[i]/SetOfReferenceValuesForPlotting["Genes"])
        FGrescaled.append(self.FG[i]/SetOfReferenceValuesForPlotting["Genes"])
        RFrescaled.append(self.RF[i]/SetOfReferenceValuesForPlotting["mRNA"])
        RHPrescaled.append(self.RHP[i]/SetOfReferenceValuesForPlotting["mRNA"])
        HPrescaled.append(self.HP[i]/SetOfReferenceValuesForPlotting["HSP"])
        i=i+1

    if ZoomInPanelA=="No":
        ax1 = plt.subplot(321)
        ax1.set_title("Protein", fontsize = "small")
        #SubPlot(ax1, self.t, [["P", self.P], [r"P$^\#$", self.Ph]],
        #        'Time (min)', tmin, tmax, 'Concentration (mM)', 0., 100.5, "center right", "A", LineStyleListInverted = "Yes")
        SubPlot(ax1, self.t, [["P", Prescaled], [r"P$^\#$", Phrescaled]],
                'Time (min)', tmin, tmax, 'Concentration (a.u.)', 0., 100.5/SetOfReferenceValuesForPlotting["Ptot"], "center right", "A", LineStyleListInverted = "Yes")

    elif ZoomInPanelA=="Yes" or ZoomInPanelA=="Yes2" or ZoomInPanelA=="Yes3" :

        CutY = []
        #for i in range(size(self.Ph)-indexoft0):
        #    CutY.append(self.Ph[i+indexoft0])
        for i in range(size(Phrescaled)-indexoft0):
            CutY.append(Phrescaled[i+indexoft0])
        MaxY = 1.4 * max(CutY)

        ax1 = plt.subplot(321)
        ax1.set_title("Protein", fontsize = "small")
        #SubPlot(ax1, self.t, [[r"P$^\#$", self.Ph]],
        #        'Time (min)', tmin, tmax, 'Concentration (mM)', 0., MaxY, "center right", "A", LineStyleListInverted = "No")
        SubPlot(ax1, self.t, [[r"P$^\#$", Phrescaled]],
                'Time (min)', tmin, tmax, 'Concentration (a.u.)', 0., MaxY, "center right", "A", LineStyleListInverted = "No")

    #elif ZoomInPanelA=="Yes":
    #    ax1 = plt.subplot(321)
    #    ax1.set_title("Protein", fontsize = "small")
    #    SubPlot(ax1, self.t, [[r"P$^\#$", self.Ph]],
    #            'Time (min)', tmin, tmax, 'Concentration (mM)', 0., 1.2, "center right", "A", LineStyleListInverted = "Yes")
    #    #AddLetterToSubplot(ax1, "A", -0.2, 1.065)
    #
    #elif ZoomInPanelA=="Yes2":
    #    ax1 = plt.subplot(321)
    #    ax1.set_title("Protein", fontsize = "small")
    #    SubPlot(ax1, self.t, [[r"P$^\#$", self.Ph]],
    #            'Time (min)', tmin, tmax, 'Concentration (mM)', 0., 4.5, "center right", "A", LineStyleListInverted = "Yes")
    #    #AddLetterToSubplot(ax1, "A", -0.2, 1.065)

    ax2 = plt.subplot(322)
    ax2.set_title("Stress Kinases", fontsize = "small")
    #SubPlot(ax2, self.t, [[r"SK", self.S], [r"SK$^*$", self.Ss]],
    #        'Time (min)', tmin, tmax, 'Concentration (nM)', 0., 120., "center right", "B", LineStyleListInverted = "Yes")
    SubPlot(ax2, self.t, [[r"SK", Srescaled], [r"SK$^*$", Ssrescaled]],
            'Time (min)', tmin, tmax, 'Concentration (a.u.)', 0., 120./SetOfReferenceValuesForPlotting["SK"], "center right", "B", LineStyleListInverted = "Yes")

    ax3 = plt.subplot(323)
    ax3.set_title("Heat Shock Factor", fontsize = "small")
    SubPlot(ax3, self.t, [[r"HSF", Frescaled], [r"HSF$^*$", Fsrescaled]],
            'Time (min)', tmin, tmax, r'Concentration (a.u.)', 0, 0, "center right", "C", LineStyleListInverted = "Yes")

    ax4 = plt.subplot(324)
    ax4.set_title("Gene", fontsize = "small")
    SubPlot(ax4, self.t, [[r"G", Grescaled], [r"HSF$^*$G", FsGrescaled], [r"HSFG", FGrescaled]],
            'Time (min)', tmin, tmax, 'Concentration (a.u.)', 0, 0, "center right", "D", LineStyleListInverted = "Yes")

    ax5 = plt.subplot(325)
    ax5.set_title("mRNA", fontsize = "small")
    SubPlot(ax5, self.t, [[r"mR$_{F}$", RFrescaled], [r"mR$_{HP}$", RHPrescaled]],
            'Time (min)', tmin, tmax, r'Concentration (a.u.)', 0, 0, "center right", "E", LineStyleListInverted = "Yes")

    ax6 = plt.subplot(326)
    ax6.set_title("Heat Shock Protein", fontsize = "small")
    SubPlot(ax6, self.t, [[r"HP", HPrescaled]],
            'Time (min)', tmin, tmax, r'Concentration (a.u.)', 0, 0, "center right", "F", LineStyleListInverted = "No")

    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        MakeGrayBackgroundTemperature(ax, self.timeset, self.TparamSet)

    # Make the time measured in hours for certain plots
    if self.TparamSet.CurrentParams["Ttype"] == 2 and self.TparamSet.CurrentParams["tb"] >= 12 * 60. * 60. + vorl:
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.set_xticks([0., 6. * 60., 12. * 60., 18. * 60., 24 * 60., (24. + 8.) * 60.])
            ax.set_xticklabels(["0", "6", "12", "18", "24", "24+8"])
            ax.set_xlabel('Time (h)', fontsize="small")
            ax.set_xlim(0. - 60., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)

    if self.TparamSet.CurrentParams["Ttype"] == 7:
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.set_xticks([0., 3. * 60., 6. * 60., 9. * 60., 12. * 60., 15. * 60., 18. * 60., 21 * 60., 24 * 60.])
            if tminMANUAL[0] == "Yes":            
                ax.set_xticks([0.+tminMANUAL[1]/60., 3. * 60.+tminMANUAL[1]/60., 6. * 60.+tminMANUAL[1]/60., 9. * 60.+tminMANUAL[1]/60., 12. * 60.+tminMANUAL[1]/60., 15. * 60.+tminMANUAL[1]/60., 18. * 60.+tminMANUAL[1]/60., 21 * 60.+tminMANUAL[1]/60., 24 * 60.+tminMANUAL[1]/60.])
            ax.set_xticklabels(["3:00", "6:00", "9:00", "12:00", "15:00", "18:00", "21:00", "24:00", "3:00"],
                               fontsize="x-small")
            ax.set_xlabel('Time (h)', fontsize="small")
            ax.set_xlim(tmin / 60., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)

    if ModifyKset == "No":
        PlotAndSave(fig1, "HeatShockResponse_" + self.name, "PS", 1, 0)
    elif ModifyKset == "Yes":
        MyString = ""
        for key in ParamsToBeModif:
            MyString = MyString + key + str(ParamsToBeModif[key])
        for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
            ax.axvline(x=((TimeOfMod - vorl) / 60.), linestyle='--', color='black')
        PlotAndSave(fig1, "HeatShockResponse_"  + MyString + "_" + self.name, "PS", 1, 0)


def FuncTimeRun(self, DirectControlnuPp, AvoidPlots, ModifyKset="No", ParamsToBeModif={}, TimeOfMod=vorl, Add3rdmRNA_ARS="No", ZoomInPanelA="No", tminMANUAL=("No",123456789)):
    """Solve ODEs system for given parameter values, plot time behaviour of T and concentrations"""
    #print("!!! In FuncTimeRun   " + str(DirectControlnuPp))
    #print("Using 12 equations!")
    SolveODEsystem(self, DirectControlnuPp, ModifyKset, ParamsToBeModif, TimeOfMod, Add3rdmRNA_ARS)
    if AvoidPlots == "No":
        self.PlotTemperature(tminMANUAL)
        PlotTrajectories(self, ModifyKset, ParamsToBeModif, TimeOfMod, ZoomInPanelA, tminMANUAL)
    elif AvoidPlots == "Yes":
        pass
    else:
        print("Error in FuncTimeRun")


def FuncTimeRun10eqs(self, DirectControlnuPp, AvoidPlots):
    """Solve ODEs system for given parameter values, plot time behaviour of T and concentrations"""
    #print("!!! In FuncTimeRun   " + str(DirectControlnuPp))
    #print("Using 10 equations!")
    SolveODEsystem10eqs(self, DirectControlnuPp)
    if AvoidPlots == "No":
        self.PlotTemperature()
        PlotTrajectories(self)
    elif AvoidPlots == "Yes":
        pass
    else:
        print("Error in FuncTimeRun10eqs")


def FuncTimeRun9eqs(self, DirectControlnuPp, AvoidPlots):
    """Solve ODEs system for given parameter values, plot time behaviour of T and concentrations"""
    #print("!!! In FuncTimeRun   " + str(DirectControlnuPp))
    #print("Using 9 equations!")
    SolveODEsystem9eqs(self, DirectControlnuPp)
    if AvoidPlots == "No":
        self.PlotTemperature()
        PlotTrajectories(self)
    elif AvoidPlots == "Yes":
        pass
    else:
        print("Error in FuncTimeRun10eqs")


def FuncFeedingExperimentStaur2Dplots(self, ParamName, ListOf3Values):
    """ Plots mRNAf and Fs for different values of ParamName, reproduces staurosporine plots """

    ListRF, ListFs = [], []
    i = 0
    for value in ListOf3Values:
        self.model.ParamsSetRATES.UpdatePar({ParamName: ListOf3Values[i]})
        SolveODEsystem(self, ("No",0.))
        ListRF.append(np.copy(self.RF))
        ListFs.append(np.copy(self.Fs))
        i = i + 1

    self.PlotTemperature()

    fig2 = figure()
    plt.plot(self.t, ListRF[0], 'r', linewidth=1., label=r"k$_S$' = %r s$^{-1}$" % ListOf3Values[0])
    plt.plot(self.t, ListRF[1], 'b', linewidth=1., label=r"k$_S$' = %r s$^{-1}$" % ListOf3Values[1])
    plt.plot(self.t, ListRF[2], 'black', linewidth=1., label=r"k$_S$' = %r s$^{-1}$" % ListOf3Values[2])
    plt.xlim(0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    plt.xlabel('Time (min)', fontsize=18)
    plt.ylabel(r'mRNA$_{F}$ ($\mu$M)', fontsize=18)
    plt.legend(loc="upper right", fontsize="medium", fancybox=True)

    MakeGrayBackgroundTemperature(plt, self.timeset, self.TparamSet)
    PlotAndSave(fig2, "FeedingExperimentRF_" + self.name, "PS", 1, 0)

    fig3 = figure()
    plt.plot(self.t, ListFs[0], 'r', linewidth=1., label=r"k$_S$' = %r s$^{-1}$" % ListOf3Values[0])
    plt.plot(self.t, ListFs[1], 'b', linewidth=1., label=r"k$_S$' = %r s$^{-1}$" % ListOf3Values[1])
    plt.plot(self.t, ListFs[2], 'black', linewidth=1., label=r"k$_S$' = %r s$^{-1}$" % ListOf3Values[2])
    plt.xlim(0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    plt.xlabel('Time (min)', fontsize=18)
    plt.ylabel(r'phosphorylated heatshock factor ($\mu$M)', fontsize=18)
    plt.legend(loc="upper right", fontsize="medium", fancybox=True)

    MakeGrayBackgroundTemperature(plt, self.timeset, self.TparamSet)
    PlotAndSave(fig3, "FeedingExperimentF_" + self.name, "PS", 1, 0)

    self.model.ParamsSetRATES.RestoreDefaultParams()


def FuncFeedingExperiment3Dplots(self, switchHPmRNAF, ParamName, ParamLabel, ParMin, ParMAX, ParStepsNumber, Camera):
    """ Make a 3D plot of some concentrations as functions of time while variating ParamName between ParMin and ParMAX """
    # switchHPmRNAF = 0 --> HP, switchHPmRNAF = 1 --> RF

    ParStep = (ParMAX - ParMin) / ParStepsNumber
    ListOfNValues = []
    for i in range(0, ParStepsNumber + 1):
        ListOfNValues.append(ParMin + i * ParStep)

    ZvalList = []
    tis0index = 123456789
    for val in ListOfNValues:
        self.model.ParamsSetRATES.UpdatePar({ParamName: val})
        SolveODEsystem(self, ("No",0.))
        if switchHPmRNAF == 0:  # Selects which variable to plot on Z
            Zval = np.copy(self.HP)
            ConcentrationLabel = "HP (mM)"
        elif switchHPmRNAF == 1:
            Zval = np.copy(self.RF)
            ConcentrationLabel = r"mR$_{F}$ ($\mu$M)"
        else:
            print("Index error in FeedingExperiment3Dplots!")
        if tis0index == 123456789:  # Cut Z vector to have only times >= 0 for plot
            for i in range(1, len(self.t)):
                if self.t[i] >= 0:
                    tis0index = i
                    break
        ZvalCUT = Zval[tis0index:].copy()
        ZvalList.append(ZvalCUT)

    # Cut the time vector to have only times >= 0 for plot
    tCUT = self.t[tis0index:].copy()

    fig4 = plt.figure()
    axf4 = fig4.add_subplot(111, projection='3d')

    # ---> x (obtain oriantation of x vector required)
    a_list_of_lists = tCUT.tolist()
    a_flattened = [val for sublist in a_list_of_lists for val in sublist]
    X_timeList = []
    for val in ListOfNValues:
        X_timeList.append(a_flattened)
    X_time = np.array(X_timeList)

    # |
    # V  y (obtain oriantation of y vector required)
    ListOfNLists = []
    for j in range(0, len(ListOfNValues)):
        ListOfNLists.append([])

    i = 0
    while i < len(a_flattened):
        for j in range(0, len(ListOfNValues)):
            ListOfNLists[j].append(ListOfNValues[j])
        i = i + 1

    Y_k = np.array(ListOfNLists)

    # ---> x
    # |
    # V  y (obtain oriantation of x and y in a matrix to finally add Z)
    a_list_of_lists = []
    a_flattened = []
    for k in range(0, len(ListOfNValues)):
        a_list_of_lists.append(ZvalList[k].tolist())
        a_flattened.append([val for sublist in a_list_of_lists[k] for val in sublist])
    Z_Concentration = np.array(a_flattened)

    # axf4.plot_wireframe(X_time, Y_k, Z_Concentration, rstride=1, cstride=25)
    axf4.plot_surface(X_time, Y_k, Z_Concentration, cmap=cm.jet, rstride=1, cstride=5, linewidth=0.2)

    axf4.set_xlabel("time (min)")
    axf4.set_ylabel(ParamLabel)
    axf4.set_zlabel(ConcentrationLabel)
    axf4.view_init(elev=Camera[0], azim=Camera[1])  # Set camera angle for optimal view

    PlotAndSave(fig4, "FeedingExperiment3D_F_" + self.name, "PS", 0, 0)

    self.model.ParamsSetRATES.RestoreDefaultParams()


def FuncFeedingExperimentPlotsVsData(self, ParamName, ListOfValues,
                                     ModelLegend, YVariabToPlot, DataFileName, DataLegend, Ylabel, NameOfFigure,
                                     LegendPosition, ColumnNumber, MultipleT="No", TParamsToUpdate={}):
    """ Simulate feeding experiment by variating a k between two values, and compare with literature data """

    Yval, YvalNORM, ListForPlotting = [], [], []
    i = 0
    for val in ListOfValues:
        if MultipleT == "No":
            self.model.ParamsSetRATES.UpdatePar({ParamName: ListOfValues[i]})
            SolveODEsystem(self, ("No",0.))
            Yval.append(np.copy(YVariabToPlot))
        if MultipleT == "Yes":
            SolveODEsystem(self, ("No",0.), ModifyKset="Yes", ParamsToBeModif={ParamName: ListOfValues[i]}, TimeOfMod=vorl)
            Yval.append(np.copy(YVariabToPlot))
            self.TparamSet.UpdatePar(TParamsToUpdate)  # Update T settings for HS
            SolveODEsystem(self, ("No",0.), ModifyKset="Yes", ParamsToBeModif={ParamName: ListOfValues[i]}, TimeOfMod=vorl)
            Yval.append(np.copy(YVariabToPlot))
            self.TparamSet.RestoreDefaultParams()  # T settings for HS back to default
        i = i + 1

    self.PlotTemperature()

    Max = max(np.max(Yval[i]) for i in range(len(Yval)))

    for j in range(len(Yval)):
        YvalNORM.append(np.asarray(Yval[j]) * 100. / Max)
        ListForPlotting.append([ModelLegend[j], YvalNORM[j]])  ### ListForPlotting.append([ModelLegend[j],YvalNORM[j]])

    ListOfOutputArrays = []
    FromDataFileToArrays(DataFileName, ColumnNumber, ListOfOutputArrays)  # Read data file, put in list of arrays
    DataTimeArMOD = ListOfOutputArrays[0] + ((self.TparamSet.CurrentParams["ta"] - vorl) / 60.)  # Translate time values

    fig6 = figure()

    print()
    print(DataTimeArMOD)
    print()
    print(ListOfOutputArrays)
    print()

    ax1 = plt.subplot(121)
    SubPlot(ax1, self.t, ListForPlotting, 'Time (min)', 0.,
            (self.timeset.CurrentParams["t_stop"] - vorl) / 60., Ylabel, 0, 0, LegendPosition, "A",
            Legendfontsize="small", Legendfancybox=True)

    ax2 = plt.subplot(122)
    DataSubPlot(ax2, DataTimeArMOD, ListOfOutputArrays, 'Time (min)', 0.,
                (self.timeset.CurrentParams["t_stop"] - vorl) / 60., Ylabel, 0., 100., LegendPosition, DataLegend, "B",
                Legendfontsize="small", Legendfancybox=True)

    for ax in [ax1, ax2]:
        MakeGrayBackgroundTemperature(ax, self.timeset, self.TparamSet)
    PlotAndSave(fig6, "FeedingExperimentVsData_" + NameOfFigure + '_' + self.name, "PS", 1, 0)

    self.model.ParamsSetRATES.RestoreDefaultParams()


def FuncTimeCourseVsDataPlot(self):
    """ Compare time course with literature data """

    SolveODEsystem(self, ("No",0.))
    #Y = np.copy(self.HP)
    # RESCALING FOR PLOTTING
    Y = np.copy(self.HP)
    Yrescaled = []
    i = 0
    for element in self.t:
        Yrescaled.append(self.HP[i]/SetOfReferenceValuesForPlotting["HSP"])
        i=i+1

    self.PlotTemperature()

    ListOfDataArrays = []
    FromDataFileToArrays("DataFiles/DataMuehlhaus2011Fig3ABC.dat", 19, ListOfDataArrays)  # Charge data file
    DataTimeMOD = ListOfDataArrays[0] + ((self.TparamSet.CurrentParams["ta"] - vorl) / 60.)  # Translate time values

    fig7 = figure()

    ax1 = plt.subplot(221)
    #ax1.plot(self.t, Y, 'blue', linewidth=1., label="HP")
    ax1.plot(self.t, Yrescaled, 'blue', linewidth=1., label="HP")
    ax1.set_xlim(0. - 5., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    ax1.set_xlabel('Time (min)', fontsize="small")
    ax1.set_ylim(0., 6.5/SetOfReferenceValuesForPlotting["HSP"]) # 5.5
    #ax1.set_ylabel(r'Heatshock protein (mM)', fontsize="small")
    ax1.set_ylabel(r'Heatshock protein (a.u.)', fontsize="small")
    ax1.legend(loc='lower right', fontsize="small", fancybox=True)
    AddLetterToSubplot(ax1, "A", -0.25, 1.055)

    ax2 = plt.subplot(222)
    ax2.errorbar(DataTimeMOD, ListOfDataArrays[1], marker="s", yerr=[ListOfDataArrays[2], ListOfDataArrays[3]],
                 color='green', label="HSP70A, western blot")
    ax2.errorbar(DataTimeMOD, ListOfDataArrays[4], marker="s", yerr=[ListOfDataArrays[5], ListOfDataArrays[6]],
                 color='red', label="HSP70A, mass spectrometry")
    ax2.set_xlim(0. - 5., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    ax2.set_ylim(-1.6, 1.6)
    ax2.set_xlabel('Time (min)', fontsize="small")
    ax2.set_ylabel('z-score', fontsize="small")
    ax2.legend(loc='lower right', numpoints=1, fontsize="small", fancybox=True)
    AddLetterToSubplot(ax2, "B", -0.25, 1.055)

    ax3 = plt.subplot(223)
    ax3.errorbar(DataTimeMOD, ListOfDataArrays[7], marker="s", yerr=[ListOfDataArrays[8], ListOfDataArrays[9]],
                 color='green', label="HSP70B, western blot")
    ax3.errorbar(DataTimeMOD, ListOfDataArrays[10], marker="s", yerr=[ListOfDataArrays[11], ListOfDataArrays[12]],
                 color='red', label="HSP70B, mass spectrometry")
    ax3.set_xlim(0. - 5., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    ax3.set_ylim(-1.6, 1.6)
    ax3.set_xlabel('Time (min)', fontsize="small")
    ax3.set_ylabel('z-score', fontsize="small")
    ax3.legend(loc='lower right', numpoints=1, fontsize="small", fancybox=True)
    AddLetterToSubplot(ax3, "C", -0.25, 1.055)

    ax4 = plt.subplot(224)
    ax4.errorbar(DataTimeMOD, ListOfDataArrays[13], marker="s", yerr=[ListOfDataArrays[14], ListOfDataArrays[15]],
                 color='green', label="HSP90A, western blot")
    ax4.errorbar(DataTimeMOD, ListOfDataArrays[16], marker="s", yerr=[ListOfDataArrays[17], ListOfDataArrays[18]],
                 color='red', label="HSP90A, mass spectrometry")
    ax4.set_xlim(0. - 5., (self.timeset.CurrentParams["t_stop"] - vorl) / 60.)
    ax4.set_ylim(-1.6, 1.6)
    ax4.set_xlabel('Time (min)', fontsize="small")
    ax4.set_ylabel('z-score', fontsize="small")
    ax4.legend(loc='lower right', numpoints=1, fontsize="small", fancybox=True)
    AddLetterToSubplot(ax4, "D", -0.25, 1.055)

    for ax in [ax1, ax2, ax3, ax4]:
        MakeGrayBackgroundTemperature(ax, self.timeset, self.TparamSet)
    PlotAndSave(fig7, "TimeCoruseData_" + self.name, "PS", 1, 0)

    self.model.ParamsSetRATES.RestoreDefaultParams()


def FuncTimeRunPlusARS(self):
    """ Compare with data for single HS experiment (using ARS enzyme) """

    SolveODEsystem(self, ("No",0.), Add3rdmRNA_ARS="Yes")
    self.PlotTemperature()

    ListOfOutputArrays = []
    FromDataFileToArrays("DataFiles/DataShroda2000ARSFig6b.dat", 3,
                         ListOfOutputArrays)  # Read data file, put in list of arrays
    DataTimeArMOD = ListOfOutputArrays[0] + ((self.TparamSet.CurrentParams["ta"] - vorl) / 60.)  # Translate time values
    Ylabel = ["ARS mRNA levels (arbitrary)", r"ARS enzyme activity"]
    #DataLegend = ["ARS mRNA", "ARS activity"]
    DataLegend = ["", ""]

    # RESCALING FOR PLOTTING
    print("started ARS...")
    RHP_ARS_rescaled, HP_ARS_rescaled = [], []
    i = 0
    Max_RHP_ARS = max(self.RHP_ARS)
    Max_HP_ARS = max(self.HP_ARS)
    print("started for loop...")
    for element in self.t:
        RHP_ARS_rescaled.append(self.RHP_ARS[i]/Max_RHP_ARS)
        HP_ARS_rescaled.append(self.HP_ARS[i]/Max_HP_ARS)
        i=i+1
    print("ended for loop...")
    fig8 = figure()

    ax1 = plt.subplot(221)
    SubPlot(ax1, self.t, [["", RHP_ARS_rescaled]],
            'Time (min)', 0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60., 'mRNA for ARS enzyme ($\mu$M)', 0.,
            0., "center right", "A")
    #SubPlot(ax1, self.t, [["", self.RHP_ARS]],
    #        'Time (min)', 0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60., 'mRNA for ARS enzyme ($\mu$M)', 0.,
    #        0., "center right", "A")

    ax2 = plt.subplot(222)
    #SubPlot(ax2, self.t, [[r"", self.HP_ARS]],
    #        'Time (min)', 0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60., 'ARS enzyme ($\mu$M)', 0., 0.,
    #        "center right", "B")
    SubPlot(ax2, self.t, [[r"", HP_ARS_rescaled]],
            'Time (min)', 0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60., 'ARS enzyme ($\mu$M)', 0., 0.,
            "center right", "B")

    ax3 = plt.subplot(223)
    DataSubPlot(ax3, DataTimeArMOD, [DataTimeArMOD, ListOfOutputArrays[1]], 'Time (min)', 0.,
                (self.timeset.CurrentParams["t_stop"] - vorl) / 60., Ylabel[0], 0., 1200., "center right",
                [DataLegend[0]], "C", Legendfontsize="x-small")

    ax4 = plt.subplot(224)
    DataSubPlot(ax4, DataTimeArMOD, [DataTimeArMOD, ListOfOutputArrays[2]], 'Time (min)', 0.,
                (self.timeset.CurrentParams["t_stop"] - vorl) / 60., Ylabel[1], 0., 2000., "center right",
                [DataLegend[1]], "D", Legendfontsize="x-small")

    for ax in [ax1, ax2, ax3, ax4]:
        MakeGrayBackgroundTemperature(ax, self.timeset, self.TparamSet)

    PlotAndSave(fig8, "HeatShockARSExp_" + self.name, "PS", 1, 1)

    PlotTrajectories(self)


def FuncTimeRunPlusARSdoubleHS(self, EmptyListToBeFilled, AvoidPlots):
    """ Compare with data for double HS experiments (using ARS enzyme) """

    Yval, ListForPlotting = [], []

    Legend = ["single HS", "2nd HS 2 h after 1st", "2nd HS 3 h after 1st", "2nd HS 4 h after 1st",
              "2nd HS 5 h after 1st"]

    # Solve ODEs system with no 2nd HS
    SolveODEsystem(self, ("No",0.), Add3rdmRNA_ARS="Yes")
    Yval.append(np.copy(self.HP_ARS))

    # Solve ODEs system multiple times with 2nd HS after, 2h, 3h, 4h, 5h.
    HSduration = self.TparamSet.CurrentParams["tb"] - self.TparamSet.CurrentParams["ta"]
    for val in [2., 3., 4., 5.]:
        self.TparamSet.UpdatePar({"Ttype": 4, "tc": val * 3600. + HSduration + vorl,
                                  "td": (val * 3600. + HSduration) + HSduration + vorl})  # Update T
        SolveODEsystem(self, ("No",0.), Add3rdmRNA_ARS="Yes")
        Yval.append(deepcopy(self.HP_ARS))

    for j in range(len(Yval)):
        ListForPlotting.append(deepcopy([Legend[j], Yval[j]]))

    if AvoidPlots == "Yes":
        pass
    elif AvoidPlots == "No":
        # Commands to plot from data file
        ListOfOutputArrays = []
        FromDataFileToArrays("DataFiles/DataShroda2000ARSFig7b.dat", 6,
                             ListOfOutputArrays)  # Read data file, put in list of arrays
        DataTimeArMOD = ListOfOutputArrays[0] + ((self.TparamSet.CurrentParams["ta"] - vorl) / 60.)  # Translate time values

        fig9 = figure()

        ax1 = plt.subplot(121)
        SubPlot(ax1, self.t, ListForPlotting,
                'Time (min)', 0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60., 'ARS enzyme concentration (a.u.)', 0.,
                0., "lower right", "A")
                #'Time (min)', 0., (self.timeset.CurrentParams["t_stop"] - vorl) / 60., 'ARS enzyme concentration ($\mu$M)', 0.,
                #0., "lower right", "A")

        ax2 = plt.subplot(122)
        DataSubPlot(ax2, DataTimeArMOD, ListOfOutputArrays, 'Time (min)', 0.,
                    (self.timeset.CurrentParams["t_stop"] - vorl)/60.,r"ARS enzyme activity (nmol $\alpha$-naphtol/mg protein $\cdot$ h)", 0., 3500., "lower right",
                    Legend, "B", Legendfontsize="x-small")

        PlotAndSave(fig9, "HeatShockARSExpDoubleHS_" + self.name, "PS", 1, 1)

        self.PlotTemperature()
        PlotTrajectories(self)
    else:
        print("Error in AvoidPlots ARS time run!")

    # In case in addition to the plots you also would like the numbers, take them out as imput of the function
    OutputList = []
    OutputList.append(deepcopy(self.t))
    OutputList.append(deepcopy(ListForPlotting))

    #print("HERE\n" + str(OutputList) + "\n")
    EmptyListToBeFilled.append(deepcopy(OutputList))




