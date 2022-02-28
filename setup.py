from setuptools import setup

VERSION = '0.1.0'
DESCRIPTION = 'Python package for symbolic circuit analysis'
LONG_DESCRIPTION = 'Python package for symbolic circuit analysis'

setup(
        name="symcirc", 
        version=VERSION,
        author="Matyas Vasek",
        author_email="<matyas.vasek@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=['symcirc'],
        install_requires=['sympy'],
        license="MIT license",
)
