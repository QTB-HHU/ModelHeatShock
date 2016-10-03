class ParametersSet:
    """ A class to create objects containing a parameter set for the model """

    def __init__(self, defaultparams={}):

        """ Set the values of the parameters"""

        self.DefaultParams = {}
        self.CurrentParams = {}

        for key in defaultparams.keys():
            self.DefaultParams.setdefault(key, defaultparams[key])
            self.CurrentParams.setdefault(key, defaultparams[key])

    def UpdatePar(self, ParamsToUpdate={}):

        """ Update paramter set by changing one or more parameters """

        for key, val in ParamsToUpdate.items():
            self.CurrentParams[key] = val

    def RestoreDefaultParams(self):

        """ Restore paramter set to default """

        self.CurrentParams = self.DefaultParams.copy()
