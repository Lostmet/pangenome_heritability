class PangenomeHeritabilityError(Exception):
    """Base exception for the package"""
    pass

class InputError(PangenomeHeritabilityError):
    """Raised when input files are invalid"""
    pass

class AlignmentError(PangenomeHeritabilityError):
    """Raised when MUSCLE alignment fails"""
    pass

class WindowError(PangenomeHeritabilityError):
    """Raised when window generation fails"""
    pass

class ConversionError(PangenomeHeritabilityError):
    """Raised when PLINK conversion fails"""
    pass

class MemoryError(PangenomeHeritabilityError):
    """Raised when memory usage exceeds limits"""
    pass