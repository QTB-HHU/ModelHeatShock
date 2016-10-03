from HSM_Reactions import *


########## RIGHT MEMBERS OF ODEs, rewritten with only 10 equations to isolate those that are independent ##############

def f10eqs(t, y, ksetDict, TparamSet, REACparamSet, DirectControlnuPp, IC_PplusPp, IC_SplusSs):
    #P = y[0]
    Ph = y[0]
    #S = y[2]
    Ss = y[1]
    F = y[2]
    Fs = y[3]
    G = y[4]
    FsG = y[5]
    FG = y[6]
    RF = y[7]
    RHP = y[8]
    HP = y[9]

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

    n1 = REACparamSet["n1"]
    n2 = REACparamSet["n2"]
    P0const = REACparamSet["P0const"]
    I = REACparamSet["I"]
    T0const = REACparamSet["T0const"]
    piRFconst = REACparamSet["piRFconst"]
    piRHPconst = REACparamSet["piRHPconst"]

    PplusPpCONST = IC_PplusPp           # (microM) Initial Condition protein P
    SplusSsCONST = IC_SplusSs           # (microM) Initial Condition stresskinease S

    system = [
        #nuP(Ph, HP, kP0) - nuPp(P, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp),                  # P
        - nuP(Ph, HP, kP0) + nuPp(PplusPpCONST - Ph, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp), # Ph
        #nuS(Ss, kS) - nuSp(S, Ph, kSp0, n2, P0const),                                                    # S
        - nuS(Ss, kS) + nuSp(SplusSsCONST - Ss, Ph, kSp0, n2, P0const),                                   # Ss
        nuF(I, Fs, kF0) + piF(RF, kFpi0) + nuFGp(FG, kFGp) - nuFG(G, F, kFG) - nuFp(F, Ss, kFp0) - etaF(F, ketaF),  # F
        - nuF(I, Fs, kF0) + nuFp(F, Ss, kFp0) + nuFsGp(FsG, kFsGp) - nuFsG(G, Fs, kFsG),  # Fs
        nuFsGp(FsG, kFsGp) + nuFGp(FG, kFGp) - nuFG(G, F, kFG) - nuFsG(G, Fs, kFsG),      # G
        nuFsG(G, Fs, kFsG) + nuFs(FG, kFs) - nuFsp(FsG, I, kFsp) - nuFsGp(FsG, kFsGp),    # FsG
        nuFsp(FsG, I, kFsp) + nuFG(G, F, kFG) - nuFGp(FG, kFGp) - nuFs(FG, kFs),          # FG
        piRF(FsG, kpiRF) + piRFAddConst(piRFconst) - etaRF(RF, ketaRF),                   # RF Added const to Alex model
        piRHP(FsG, kpiRH) + piRHPAddConst(piRHPconst) - etaRHP(RHP, ketaRHP),             # RHP Aded const to Alex model
        piHP(RHP, kpiHP) - etaHP(HP, ketaHP)]                                             # HP
    # Notice presence of nuFG() in line of F, presence of nuFsG() in that of Fs, absence of pi in that of FsG.

    return system



########## RIGHT MEMBERS OF ODEs, rewritten with only 9 equations to isolate those that are independent ##############

def f9eqs(t, y, ksetDict, TparamSet, REACparamSet, DirectControlnuPp, IC_PplusPp, IC_SplusSs, IC_GplusFsGplusFG):
    #P = y[0]
    Ph = y[0]
    #S = y[2]
    Ss = y[1]
    F = y[2]
    Fs = y[3]
    #G = y[4]
    FsG = y[4]
    FG = y[5]
    RF = y[6]
    RHP = y[7]
    HP = y[8]

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

    n1 = REACparamSet["n1"]
    n2 = REACparamSet["n2"]
    P0const = REACparamSet["P0const"]
    I = REACparamSet["I"]
    T0const = REACparamSet["T0const"]
    piRFconst = REACparamSet["piRFconst"]
    piRHPconst = REACparamSet["piRHPconst"]

    PplusPpCONST = IC_PplusPp           # (microM) Initial Condition protein P
    SplusSsCONST = IC_SplusSs           # (microM) Initial Condition stresskinease S
    GplusFsGplusFG = IC_GplusFsGplusFG # (microM) Initial Condition gene G

    G = GplusFsGplusFG - FsG - FG

    system = [
        #nuP(Ph, HP, kP0) - nuPp(P, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp),                  # P
        - nuP(Ph, HP, kP0) + nuPp(PplusPpCONST - Ph, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp), # Ph
        #nuS(Ss, kS) - nuSp(S, Ph, kSp0, n2, P0const),                                                    # S
        - nuS(Ss, kS) + nuSp(SplusSsCONST - Ss, Ph, kSp0, n2, P0const),                                   # Ss
        nuF(I, Fs, kF0) + piF(RF, kFpi0) + nuFGp(FG, kFGp) - nuFG(G, F, kFG) - nuFp(F, Ss, kFp0) - etaF(F, ketaF),  # F
        - nuF(I, Fs, kF0) + nuFp(F, Ss, kFp0) + nuFsGp(FsG, kFsGp) - nuFsG(G, Fs, kFsG),  # Fs
        #nuFsGp(FsG, kFsGp) + nuFGp(FG, kFGp) - nuFG(G, F, kFG) - nuFsG(G, Fs, kFsG),      # G
        nuFsG(G, Fs, kFsG) + nuFs(FG, kFs) - nuFsp(FsG, I, kFsp) - nuFsGp(FsG, kFsGp),    # FsG
        nuFsp(FsG, I, kFsp) + nuFG(G, F, kFG) - nuFGp(FG, kFGp) - nuFs(FG, kFs),          # FG
        piRF(FsG, kpiRF) + piRFAddConst(piRFconst) - etaRF(RF, ketaRF),                   # RF Added const to Alex model
        piRHP(FsG, kpiRH) + piRHPAddConst(piRHPconst) - etaRHP(RHP, ketaRHP),             # RHP Aded const to Alex model
        piHP(RHP, kpiHP) - etaHP(HP, ketaHP)]                                             # HP
    # Notice presence of nuFG() in line of F, presence of nuFsG() in that of Fs, absence of pi in that of FsG.

    return system







