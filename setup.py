import setuptools
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='motifsearch',
    version='1.0',
    packages=setuptools.find_packages(),
)