########################################################################################################################

__doc__ = \
    """
    See GitHub documentation
    """
__author__ = "Kate Fieseler. [Github](https://github.com/kate-fie)"
__email__ = "kate.k.fieseler@gmail.com"
__date__ = "2024"
__license__ = "MIT"
__citation__ = ""
from .version import __version__

##################################################################################################################

if __name__ == "__main__":
    from .cli import main
    main()