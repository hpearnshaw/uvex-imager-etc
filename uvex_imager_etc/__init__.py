from importlib.metadata import version as _version, PackageNotFoundError
try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass

# Suppress certain warnings
import warnings
from erfa import ErfaWarning
warnings.simplefilter('ignore', ErfaWarning)
