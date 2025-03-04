from .ticlust import *
try:
    from ._version import version as __version__
except ImportError:
    from setuptools_scm import get_version
    __version__ = get_version(root="..", relative_to=__file__)


__all__ = ['__version__', 'ticlust']
