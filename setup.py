from setuptools import setup

VERSION = '1.0.0' 
DESCRIPTION = 'My first Python package'
LONG_DESCRIPTION = 'My first Python package with a slightly longer description'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="symcirc", 
        version=VERSION,
        author="name",
        author_email="<youremail@email.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=['symcirc'],
        install_requires=['sympy'],
)
