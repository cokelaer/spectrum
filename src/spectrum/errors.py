



class SpectrumError(Exception):
    """A Base class error
    
    """
    #def __init__(self, value=None):
    #    self.value = value
    def __str__(self):
        return "Unknown Spectrum Error:"
    
    
class SpectrumModifiedError(SpectrumError):
    def __str__(self):
        msg = """Data has been modified (or NFFT, or sampling, or N) but 
            the PSD has not be (re)-computed. Call it using name.()"""
        return msg
    
class SpectrumChoiceError(SpectrumError):
    def __init__(self, value, valid):
        #super(SpectrumChoiceError, self).__init__()
        self.value = value
        self.valid = valid
    def __str__(self):
        msg = """Spectrum error: invalid choice (%s). Possible values are %s.""" %\
            (self.value, self.valid)
        return msg
    
class SpectrumPSDError(SpectrumError):
    def __str__(self):
        msg = """PSD not computed. Populate the attribute psd."""
        return msg
    
class SpectrumARMAError(SpectrumError):
    def __str__(self):
        msg = """Either AR parameter and/or MA must be provided"""
        return msg

class SpectrumARError(SpectrumError):
    def __str__(self):
        msg = """AR order must be stricly positive"""
        return msg    
    
class SpectrumMAError(SpectrumError):
    def __str__(self):
        msg = """AR order must be stricly positive"""
        return msg
    
class SpectrumOrder(SpectrumError):
    def __str__(self):
        msg = """AR or MA order must be stricly positive"""
        return msg

class SpectrumNPSD(SpectrumError):
    def __str__(self):
        msg = """NPSD must be stricly positive"""
        return msg
        
def is_positive_integer(order, name="wrong argument. "):
    if order < 0:
        raise SpectrumOrder
    if type(order) != int:
        raise TypeError("%s must be an integer" % name)
    return True