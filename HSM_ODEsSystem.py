from HSM_Reactions import *


########## RIGHT MEMBERS OF ODEs ##############

def f(t, y, ksetDict, TparamSet, REACparamSet, DirectControlnuPp, Add3rdmRNA_ARS="No"): # , RemoveTwoRedundantEquations="No"
    P = y[0]
    Ph = y[1]
    S = y[2]
    Ss = y[3]
    F = y[4]
    Fs = y[5]
    G = y[6]
    FsG = y[7]
    FG = y[8]
    RF = y[9]
    RHP = y[10]
    HP = y[11]
    if Add3rdmRNA_ARS == "Yes":
        RHP_ARS = y[12]
        HP_ARS = y[13]

    kP0 = ksetDict["kP0"]
    kP0p = ksetDict["kP0p"]
    kS = ksetDict["kS"]
    kSp0 = ksetDict["kSp0"]
    kFp0 = ksetDict["kFp0"]
    kF0 = ksetDict["kF0"]
    kFpi0 = ksetDict["kFpi0"]
    kFGp = ksetDict["kFGp"]
    kFG = ksetDict["kFG"]
    ketaF = ksetDict["ketaF"]
    kFsG = ksetDict["kFsG"]
    kFsGp = ksetDict["kFsGp"]
    kFsp = ksetDict["kFsp"]
    kFs = ksetDict["kFs"]
    kpiRF = ksetDict["kpiRF"]
    kpiRH = ksetDict["kpiRH"]
    kpiHP = ksetDict["kpiHP"]
    ketaHP = ksetDict["ketaHP"]
    ketaRF = ksetDict["ketaRF"]
    ketaRHP = ksetDict["ketaRHP"]
    if Add3rdmRNA_ARS == "Yes":
        kpiRHP_ARS = 9. * 4.            # F*G: mRHP    (s^-1)
        kpiHP_ARS = 0.5 * 5. / 20 / 1000.      # mRHP: HP     (s^-1)
        ketaHP_ARS = 0.00005 / 8. * 2   # HP-->        (s^-1)
        ketaRHP_ARS = 0.006 / 4. * 0.8  # mRHP-->      (s^-1)

    n1 = REACparamSet["n1"]
    n2 = REACparamSet["n2"]
    P0const = REACparamSet["P0const"]
    I = REACparamSet["I"]
    T0const = REACparamSet["T0const"]
    piRFconst = REACparamSet["piRFconst"]
    piRHPconst = REACparamSet["piRHPconst"]

    system = [
        nuP(Ph, HP, kP0) - nuPp(P, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp),                      # P
        - nuP(Ph, HP, kP0) + nuPp(P, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp),                    # Ph
        nuS(Ss, kS) - nuSp(S, Ph, kSp0, n2, P0const),                                     # S
        - nuS(Ss, kS) + nuSp(S, Ph, kSp0, n2, P0const),                                   # Ss
        nuF(I, Fs, kF0) + piF(RF, kFpi0) + nuFGp(FG, kFGp) - nuFG(G, F, kFG) - nuFp(F, Ss, kFp0) - etaF(F, ketaF),  # F
        - nuF(I, Fs, kF0) + nuFp(F, Ss, kFp0) + nuFsGp(FsG, kFsGp) - nuFsG(G, Fs, kFsG),  # Fs
        nuFsGp(FsG, kFsGp) + nuFGp(FG, kFGp) - nuFG(G, F, kFG) - nuFsG(G, Fs, kFsG),      # G
        nuFsG(G, Fs, kFsG) + nuFs(FG, kFs) - nuFsp(FsG, I, kFsp) - nuFsGp(FsG, kFsGp),    # FsG
        nuFsp(FsG, I, kFsp) + nuFG(G, F, kFG) - nuFGp(FG, kFGp) - nuFs(FG, kFs),          # FG
        piRF(FsG, kpiRF) + piRFAddConst(piRFconst) - etaRF(RF, ketaRF),                   # RF Added const to Alex model
        piRHP(FsG, kpiRH) + piRHPAddConst(piRHPconst) - etaRHP(RHP, ketaRHP),             # RHP Aded const to Alex model
        piHP(RHP, kpiHP) - etaHP(HP, ketaHP)]                                             # HP
    # Notice presence of nuFG() in line of F, presence of nuFsG() in that of Fs, absence of pi in that of FsG.

    if Add3rdmRNA_ARS == "Yes":
        system.append(piRHP_ARS(FsG, kpiRHP_ARS) - etaRHP_ARS(RHP_ARS, ketaRHP_ARS))     # mRNA_HP_ARS (RHPARS)
        system.append(piHP_ARS(RHP_ARS, kpiHP_ARS) - etaHP_ARS(HP_ARS, ketaHP_ARS))      # HP_ARS (HPARS)

    #print("The ODEs system has been called!!!   " + str(system) )

    return system

