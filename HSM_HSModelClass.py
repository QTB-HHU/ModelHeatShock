from scipy.integrate import ode
import matplotlib.pyplot as plt
import numpy as np

from HSM_ODEsSystem import *
from HSM_ParametersClass import *


class HeatShockModel:
    """ Class to model the Heat shock model """

    def __init__(self, ParametersSetIC, ParametersSetRATES, ParamSetForREACTIONS):
        """ Set the values of the parameters, i.e. the rates """

        self.ParamsSetIC = ParametersSetIC
        self.ParamsSetRATES = ParametersSetRATES
        self.ParamsSetForREAC = ParamSetForREACTIONS
