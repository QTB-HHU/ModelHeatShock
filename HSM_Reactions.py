from HSM_Temperature import *


############ REACTIONS ############

def nuP(Ph, HP, kP0):  # P#-->P
    nuP = kP0 * Ph * HP  # //(2000+HP)
    return nuP


def nuPp(P, t, kP0p, n1, T0const, TparamSet, DirectControlnuPp):  # P-->P#
    # the DirectControl optional argument serves to switch between the normal nuPp and the nuPp that we change directly.
    if DirectControlnuPp[0] == "No":
        #n1 and T0const are not used anymore

        R = 8.3144598 # Ideal Gas Constant (J mol^-1 K^-1)
        Ea = 174440. # (J mol^-1) Activation energy (J mol^-1)

        B = Ea/R # B = 20980.330
        A = kP0p * 9.431831774375398 * pow(10.,28) # 9.27*pow(10,30)  (kP0p was 98.28419570824974)

        KelvinToCelsius = 273.15 
        TinKelvin = T(t, TparamSet) + KelvinToCelsius

        nuPp = P  *  A*math.exp(-B/TinKelvin) # P * (Arrenius Equation for the Temp dependent k)

        #nuPp = kP0p * P * T(t, TparamSet) ** n1 / (T0const ** n1 + T(t, TparamSet) ** n1)

        #print( "!!!   " + str((t-vorl)/60.) + "   " + str(T(t, TparamSet)) + "   " + str(nuPp)) 
        #print()
        #print(DirectControlnuPp[0]+"   "+str(nuPp))
    elif DirectControlnuPp[0] == "Yes":
        nuPp = DirectControlnuPp[1] * P
        #print( "XXXXXXXXXXXXXXXXXXXXXXX   " + str((t-vorl)/60.) + "   " + str(T(t, TparamSet)) + "   " + str(nuPp)) 
        #print()
        #print(DirectControlnuPp[0]+"   "+str(nuPp))
    else:
        print("Error in nuPp in HSM_Reactions.py")
    return nuPp



def nuS(Ss, kS):  # S*-->S
    nuS = kS * Ss
    return nuS


def nuSp(S, Ph, kSp0, n2, P0const):  # S-->S*
    nuSp = kSp0 * S * pow(Ph, n2) / (pow(P0const, n2) + pow(Ph, n2))
    return nuSp


def nuFp(F, Ss, kFp0):  # F-->F*
    nuFp = kFp0 * F * Ss / (1. + Ss)
    return nuFp


def nuF(I, Fs, kF0):  # F*-->F //wird wohl michaelis Menten sein ?
    nuF = kF0 * I * Fs
    return nuF


def piF(RF, kFpi0):  # mRF: F + mRF
    piF = kFpi0 * RF
    return piF


def nuFGp(FG, kFGp):  # FG -> F + G
    nuFGp = kFGp * FG
    return nuFGp


def nuFG(G, F, kFG):  # F +G -> FG  //PROBLEM!!! WIE REAKTION (= how reactive) MIT 2 REAKTANTEN?
    nuFG = kFG * G * F
    return nuFG


def etaF(F, ketaF):  # F-->
    etaF = ketaF * F
    return etaF


def nuFsG(G, Fs, kFsG):  # F* + G -> F*G //PROBLEM!!! sehe oben (=see above)!?
    nuFsG = kFsG * G * Fs
    return nuFsG


def nuFsGp(FsG, kFsGp):  # F*G->F* + G
    nuFsGp = kFsGp * FsG
    return nuFsGp


def nuFsp(FsG, I, kFsp):  # F*G->FG
    nuFsp = kFsp * FsG * I
    return nuFsp


def nuFs(FG, kFs):  # FG->F*G
    nuFs = kFs * FG
    return nuFs


def piRF(FsG, kpiRF):  # F*G: mRF
    piRF = kpiRF * FsG
    return piRF


def piRHP(FsG, kpiRH):  # F*G: mRHP
    piRHP = kpiRH * FsG
    return piRHP


def piRFAddConst(piRFconst):
    return piRFconst


def piRHPAddConst(piRHPconst):
    return piRHPconst


def piHP(RHP, kpiHP):  # mRHP: HP + mRHP
    piHP = kpiHP * RHP
    return piHP


def etaHP(HP, ketaHP):  # HP-->
    etaHP = ketaHP * HP
    return etaHP


def etaRF(RF, ketaRF):  # mRF-->
    etaRF = ketaRF * RF
    return etaRF


def etaRHP(RHP, ketaRHP):  # mRHP-->
    etaRHP = ketaRHP * RHP
    return etaRHP


# The following 4 reactions, not present in the original model by Alexander, can be added to the system
# to reproduce the results of the experiments performed in Schroda et al. 2000
def piRHP_ARS(FsG, kpiRH_ARS):  # F*G: mRHP_ARS
    piRHP_ARS = kpiRH_ARS * FsG
    return piRHP_ARS


def piHP_ARS(RHP_ARS, kpiHP_ARS):  # mRHP_ARS: HP_ARS + mRHP_ARS
    piHP_ARS = kpiHP_ARS * RHP_ARS
    return piHP_ARS


def etaHP_ARS(HP_ARS, ketaHP_ARS):  # HP_ARS-->
    etaHP_ARS = ketaHP_ARS * HP_ARS
    return etaHP_ARS


def etaRHP_ARS(RHP_ARS, ketaRHP_ARS):  # mRHP_ARS-->
    etaRHP_ARS = ketaRHP_ARS * RHP_ARS
    return etaRHP_ARS
