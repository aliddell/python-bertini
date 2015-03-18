class NAGlibBaseException(Exception):
    """
    Base exception type
    """
    def __init__(self, message):
        self.message = message
    
    def __str__(self):
        return self.message
    
class AffineInfinityException(NAGlibBaseException):
    """
    """
    def __init__(self, message):
        super(AffineInfinityException, self).__init__(message)
    
class BertiniError(Exception):
    """
    BertiniError
    
    Raise BertiniError when Bertini returns nonzero exit status
    """
    def __init__(self, message):
        super(BertiniError, self).__init__(message)
    
class ExitSpaceError(NAGlibBaseException):
    """
    """
    def __init__(self, message):
        super(ExitSpaceError, self).__init__(message)

class NoBertiniException(NAGlibBaseException):
    """
    NoBertiniException
    
    Raise NoBertiniException when Bertini can't be located on the system
    """
    def __init__(self):
        message = "you don't seem to have Bertini installed anywhere I can find it"
        super(NoBertiniException, self).__init__(message)
    
class NonPolynomialException(NAGlibBaseException):
    """
    NonPolynomialException
    
    Raise NonPolynomialException when user attempts to instantiate a
    PolynomialSystem or LinearSystem object with non-polynomials
    """
    def __init__(self, message):
        super(NonPolynomialException, self).__init__(message)
    
class NonHomogeneousException(NAGlibBaseException):
    """
    NonHomogeneousException
    
    Raise NonHomogeneousException when user attempts to instantiate
    a projective system with nonhomogeneous polynomials
    """
    def __init__(self, message):
        super(NonHomogeneousException, self).__init__(message)
        
class UnclassifiedException(NAGlibBaseException):
    """
    """
    def __init__(self, message):
        super(UnclassifiedException, self).__init__(message)
        
class WitnessDataException(NAGlibBaseException):
    """
    """
    def __init__(self, message):
        super(WitnessDataException, self).__init__(message)