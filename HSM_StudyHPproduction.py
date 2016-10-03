
import matplotlib.pyplot as plt

from HSM_SimulateClass import *

def StudyHPproductionForDifferentTemeperaturesDurationsHS(HSModel,TEMPERATUREstart,TEMPERATUREstop,NstepsTEMPERATURE,
                                                          DURATIONstart, DURATIONstop, NstepsDURATION, FigureName):
    """Study the production of HS protein for different HS temeperatures and HS durations, and plot in a 2D color/comtours map"""

    ############### PART I: COMPUTE HP BY RUNNING A SIMULATION ###############

    DTEMPERATURE = (TEMPERATUREstop-TEMPERATUREstart) / (NstepsTEMPERATURE - 1)
    DDURATION = (math.log10(DURATIONstop)-math.log10(DURATIONstart)) / (NstepsDURATION - 1)

    # Prepare list and arrays to host values of HS duration (x), HS temperature (y) and HP concentration at HS's end (z).
    x = [None] * (NstepsDURATION)
    y = []
    zpart = [None] * (NstepsDURATION)
    z_flattened = []
    z_flattenedNotArray = []
    z_list = []

    # compute z on a rectangular grid in x (loop over j) and y (loop over i). 
    for i in range(NstepsTEMPERATURE):
        TEMPERATURE = TEMPERATUREstart + i*DTEMPERATURE

        for j in range(NstepsDURATION):
            DURATION = DURATIONstart * math.pow(10,j*DDURATION)

            # for each couple x,y run a simulation...
            Tset = ParametersSet({"Ttype": 1, "Tin": 20., "Tup": TEMPERATURE, "tau": 5., "ta": vorl})
            timeset = ParametersSet({"t_start": 0., "t_stop": DURATION + vorl, "delta_t": 5.0})
            Simulation = Simulate(HSModel, timeset, Tset, "Name_Not_Used")
            Simulation.TimeRun(AvoidPlots="Yes")

            # ...and store the values of time (duration of HS) and HP concentration, and...
            x[j]=Simulation.t[-1] 
            zpart[j]=Simulation.HP[-1]

            # print(str(i) + "  " + str(j) + "  " + str(TEMPERATURE) + "  " + str(DURATION) + "  " 
            #       + str(Simulation.t[-1]) + "  " + str(Simulation.Tplot[-1]) + "  " + str(Simulation.HP[-1]) )

        # ...store the values of Temperature as well.
        y.append(Simulation.Tplot[-1])
        # Flatten z and combine the parts into a single list
        zpart_flattened = [val for sublist in zpart for val in sublist]
        z_list.append(zpart_flattened)

    # Provide to x, y, z the vector/matricial form appropriate to use them with plt.contour
    z_flattened = np.vstack(z_list)
    x_flattened = np.array([math.log10(val) for sublist in x for val in sublist])
    y_flattened = np.array([val for sublist in y for val in sublist])
    
    # print("\n"+str(x_flattened)+"\n"+str(y_flattened)+"\n"+str(z_flattened)
    #       +"\n"+"\n"+str(x_flattened.shape)+"\n"+str(y_flattened.shape)+"\n"+str(z_flattened.shape)+"\n")

    ############### PART II: PLOT THE RESULTS ###############

    # RESCALING FOR PLOTTING
    HP_rescaled = []
    i = 0
    for element in z_flattened:
        HP_rescaled.append(z_flattened[i]/SetOfReferenceValuesForPlotting["HSP"])
        i=i+1


    fig = figure()
    
    # Set axis limits and labels
    plt.xlim(math.log10(DURATIONstart/60.),math.log10(DURATIONstop/60.))
    plt.ylim(TEMPERATUREstart,TEMPERATUREstop)
    plt.xlabel('HS Duration (minutes)', fontsize = 18)
    plt.ylabel('HS Temeperature ($^\circ$C)', fontsize = 18)

    # draw contour lines
    #CS = plt.contour(x_flattened, y_flattened, z_flattened, 20, linewidths=1, colors='Black')
    CS = plt.contour(x_flattened, y_flattened, HP_rescaled, 20, linewidths=1, colors='Black')
    plt.clabel(CS, inline=1, fontsize=10 ,fmt = '%0.2f')

    # draw color map on the background
    #CS = plt.contourf(x_flattened, y_flattened, z_flattened, 100, cmap=plt.cm.rainbow, )
    CS = plt.contourf(x_flattened, y_flattened, HP_rescaled, 100, cmap=plt.cm.rainbow, )

    # draw colorbar
    cbar = plt.colorbar()
    #cbar.set_label('HP concentration at the end of HS (mM)', fontsize="medium")
    cbar.set_label('HP concentration at the end of HS (a.u.)', fontsize="medium")

    # Set axis' ticks at proper places (convert the axis from linear in the exponents to log)
    plt.tick_params(axis='x', labelsize=14)
    ax=plt.gca()
    #ax.set_xticks([0,1,2,3])
    #ax.set_xticklabels(["1","10","100","1000"])
    #ax.xaxis.set_ticks([math.log10(2),math.log10(3),math.log10(4),math.log10(5),math.log10(6),math.log10(7),math.log10(8),math.log10(9),
    #math.log10(20),math.log10(30),math.log10(40),math.log10(50),math.log10(60),math.log10(70),math.log10(80),math.log10(90),
    #math.log10(200),math.log10(300),math.log10(400),math.log10(500),math.log10(600),math.log10(700),math.log10(800),math.log10(900)], minor = True)

    ax.set_xticks([1,2,3])
    ax.set_xticklabels(["10","100","1000"])
    ax.xaxis.set_ticks([math.log10(20),math.log10(30),math.log10(40),math.log10(50),math.log10(60),math.log10(70),math.log10(80),math.log10(90),
    math.log10(200),math.log10(300),math.log10(400),math.log10(500),math.log10(600),math.log10(700),math.log10(800),math.log10(900)], minor = True)

    plt.show()

    PlotAndSave(fig, FigureName, "PS", 1, 0)
