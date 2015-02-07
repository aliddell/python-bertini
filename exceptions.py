class BertiniError(Exception):
    """BertiniError
    
    Raise BertiniError when Bertini returns nonzero exit status
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

class NoBertiniException(Exception):
    """
    NoBertiniException
    
    Raise NoBertiniException when Bertini can't be located on the system
    """
    def __init__(self):
        self.message = "You don't seem to have Bertini installed anywhere " \
                       "I can find it."
    def __str__(self):
        return self.message
    
class NonPolynomialException(Exception):
    """
    NonPolynomialException
    
    Raise NonPolynomialException when user attempts to instantiate a
    PolynomialSystem or LinearSystem object with non-polynomials
    """
    def __init__(self, function):
        self.message = 'function {0} is not a polynomial'.format(function)
    
    def __str__(self):
        return self.message
    
class NonLinearException(Exception):
    """
    NonLinearException
    
    Raise NonLinearException when user attempts to instantiate a
    LinearSystem object with nonlinear polynomials
    """
    def __init__(self, polynomial):
        self.message = 'polynomial {0} is nonlinear'.format(polynomial)
        
    def __str__(self):
        return self.message