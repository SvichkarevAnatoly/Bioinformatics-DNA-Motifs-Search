try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='motifsearch',
    version='1.0',
    py_modules=[
        'src/lib',
        'src/utils/bed_center_extender',
        'src/utils/logo_generator',
        'src/utils/pattern_matching',
        'src/utils/pwm_generator',
    ],
)