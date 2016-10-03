from matplotlib.pylab import *


############ TEMPERATURE ############

# IMPORTANT: note difference with the original T(t) using steps, due to .ode crashing for t>220s.
# Once the new shape for the temperature is used .ode works fine, and this is more physical.
# But be aware of the difference!!!


def T(t, TparamsSet):
    Ttype = TparamsSet.CurrentParams["Ttype"]  # 1 --> Tin, Tup   2 --> Tin, Tup, Tin   3 --> Tin, Tup, Tin, Tup
    Tin = TparamsSet.CurrentParams["Tin"]

    # REMARK!!! tau should be ways smaller than (ti-tj) approx few 1000s.
    # Around 1s to 10s. Smaller -> .ode crashes, bigger -> to big difference with step T(t).
    # Best value would be 7s for all the code (but then canavnine gives problems)



    if Ttype == 0:
        TT = Tin

    elif Ttype == 1:
        Tup = TparamsSet.CurrentParams["Tup"]
        DT = Tup - Tin
        ta = TparamsSet.CurrentParams["ta"]
        tau = TparamsSet.CurrentParams["tau"]
        if t <= ta:
            TT = Tin
        elif t > ta:
            TT = Tin + DT * (1. - exp(-(t - ta) / tau))
        else:
            print("Error in Temperature function!")

    elif Ttype == 2:
        Tup = TparamsSet.CurrentParams["Tup"]
        DT = Tup - Tin
        ta = TparamsSet.CurrentParams["ta"]
        tb = TparamsSet.CurrentParams["tb"]
        tau = TparamsSet.CurrentParams["tau"]
        if t <= ta:
            TT = Tin
        elif t > ta and t <= tb:
            TT = Tin + DT * (1. - exp(-(t - ta) / tau))
        elif t > tb:
            TT = Tin + DT * (exp(-(t - tb) / tau))
        else:
            print("Error in Temperature function!")

    elif Ttype == 3:
        Tup = TparamsSet.CurrentParams["Tup"]
        DT = Tup - Tin
        ta = TparamsSet.CurrentParams["ta"]
        tb = TparamsSet.CurrentParams["tb"]
        tc = TparamsSet.CurrentParams["tc"]
        tau = TparamsSet.CurrentParams["tau"]
        if t <= ta:
            TT = Tin
        elif t > ta and t <= tb:
            TT = Tin + DT * (1. - exp(-(t - ta) / tau))
        elif t > tb and t <= tc:
            TT = Tin + DT * (exp(-(t - tb) / tau))
        elif t > tc:
            TT = Tin + DT * (1. - exp(-(t - tc) / tau))
        else:
            print("Error in Temperature function!")

    elif Ttype == 4:
        Tup = TparamsSet.CurrentParams["Tup"]
        DT = Tup - Tin
        ta = TparamsSet.CurrentParams["ta"]
        tb = TparamsSet.CurrentParams["tb"]
        tc = TparamsSet.CurrentParams["tc"]
        td = TparamsSet.CurrentParams["td"]
        tau = TparamsSet.CurrentParams["tau"]
        if t <= ta:
            TT = Tin
        elif t > ta and t <= tb:
            TT = Tin + DT * (1. - exp(-(t - ta) / tau))
        elif t > tb and t <= tc:
            TT = Tin + DT * (exp(-(t - tb) / tau))
        elif t > tc and t <= td:
            TT = Tin + DT * (1. - exp(-(t - tc) / tau))
        elif t > td:
            TT = Tin + DT * (exp(-(t - td) / tau))
        else:
            print("Error in Temperature function!")

    elif Ttype == 5:
        Tup1 = TparamsSet.CurrentParams["Tup1"]
        Tup2 = TparamsSet.CurrentParams["Tup2"]
        DT1 = Tup1 - Tin
        DT2 = Tup2 - Tup1
        t1 = TparamsSet.CurrentParams["t1"]
        t2 = TparamsSet.CurrentParams["t2"]
        tau = TparamsSet.CurrentParams["tau"]
        if t <= t1:
            TT = Tin
        elif t > t1 and t <= t2:
            TT = Tin + DT1 * (1. - exp(-(t - t1) / tau))
        elif t > t2:
            TT = Tup1 + DT2 * (1. - exp(-(t - t2) / tau))
        else:
            print("Error in Temperature function!")

    elif Ttype == 6:
        Tup = TparamsSet.CurrentParams["Tup"]
        DT = Tup - Tin
        tUP = TparamsSet.CurrentParams["Dt"]
        NumHS = TparamsSet.CurrentParams["NumHS"]
        tau = TparamsSet.CurrentParams["tau"]
        TT = Tin

        number = int((t-vorl)/tUP)
        if t >= vorl:
            if number == 0:
              pass
            elif number % 2 == 0:
                TT = Tin + DT * (exp(-(t - (vorl + (number) * tUP)) / tau))
            else:
                TT = Tin + DT * (1. - exp(-(t - (vorl + (number) * tUP)) / tau))
        #for i in range(0, NumHS):
        #    if t > vorl + (1 + 2 * i) * tUP and t <= vorl + (2 + 2 * i) * tUP:
        #        TT = Tin + DT * (1. - exp(-(t - (vorl + (1 + 2 * i) * tUP)) / tau))
        #    elif t > vorl + (2 + 2 * i) * tUP and t <= vorl + (3 + 2 * i) * tUP:
        #        TT = Tin + DT * (exp(-(t - (vorl + (2 + 2 * i) * tUP)) / tau))

    elif Ttype == 7:
        Tup = TparamsSet.CurrentParams["Tup"]
        Period = TparamsSet.CurrentParams["Period"]
        # Phase = TparamsSet.CurrentParams["Phase"]
        DT = Tup - Tin
        if t < vorl:
            TT = Tin  # Tin+DT/2.+DT/2.*sin( (-Phase)/Period*2*3.14152 - 3.14152/2.)
        elif t >= vorl:
            TT = Tin + DT / 2. + DT / 2. * sin((t - vorl) / Period * 2 * 3.14152 - 3.14152 / 2.)
        else:
            print("Error in Temperature function!")

    else:
        print("Error, wrong Temperature index!")

    #print( str(Ttype) + "   " + str((t-vorl)/60.) + "   " + str(TT) ) 

    return TT


vorl = (2500.) * 60.  # (changed from 1 hour=60*60 in Alex's code!!!)
# time after which we start to perturb/look at the system (i.e. we assume steady state has been reached)


def MakeGrayBackgroundTemperature(ax, timeset, TparamsSet):
    """ Put gray background to the plots when the temperature is high (Heat Shock)  """

    tend = (timeset.CurrentParams["t_stop"] - vorl) / 60.
    Ttype = TparamsSet.CurrentParams["Ttype"]

    ColorOfFace = (255. / 255, 220. / 255., 203. / 255.)  # (0.8, 0.8, 0.8)
    ColorOfEdge = (255. / 255, 210. / 255., 200. / 255.)  # (0.8, 0.8, 0.8)

    if Ttype == 3:
        ta = (TparamsSet.CurrentParams["ta"] - vorl) / 60.
        tb = (TparamsSet.CurrentParams["tb"] - vorl) / 60.
        tc = (TparamsSet.CurrentParams["tc"] - vorl) / 60.
        ax.axvspan(ta, tb, facecolor=ColorOfFace, edgecolor=ColorOfEdge, linewidth=0.8, alpha=1)
        ax.axvspan(tc, tend, facecolor=ColorOfFace, edgecolor=ColorOfEdge, linewidth=0.8, alpha=1)

    elif Ttype == 2:
        ta = (TparamsSet.CurrentParams["ta"] - vorl) / 60.
        tb = (TparamsSet.CurrentParams["tb"] - vorl) / 60.
        ax.axvspan(ta, tb, facecolor=ColorOfFace, edgecolor=ColorOfEdge, linewidth=0.8, alpha=1)

    elif Ttype == 1:
        ta = (TparamsSet.CurrentParams["ta"] - vorl) / 60.
        ax.axvspan(ta, tend, facecolor=ColorOfFace, edgecolor=ColorOfEdge, linewidth=0.8, alpha=1)

    elif Ttype == 4:
        ta = (TparamsSet.CurrentParams["ta"] - vorl) / 60.
        tb = (TparamsSet.CurrentParams["tb"] - vorl) / 60.
        tc = (TparamsSet.CurrentParams["tc"] - vorl) / 60.
        td = (TparamsSet.CurrentParams["td"] - vorl) / 60.
        ax.axvspan(ta, tb, facecolor=ColorOfFace, edgecolor=ColorOfEdge, linewidth=0.8, alpha=1)
        ax.axvspan(tc, td, facecolor=ColorOfFace, edgecolor=ColorOfEdge, linewidth=0.8, alpha=1)

    elif Ttype == 0 or Ttype == 5 or Ttype == 6 or Ttype == 7:
        pass

    else:
        print("Error, wrong Temperature index in function PlotTrajectories!")
