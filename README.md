# SymCirc
[![pypi version](https://img.shields.io/pypi/v/symcirc.svg)](https://pypi.python.org/pypi/symcirc)
[![Downloads](https://static.pepy.tech/badge/symcirc)](https://pepy.tech/project/symcirc)
[![Downloads](https://static.pepy.tech/badge/symcirc/month)](https://pepy.tech/project/symcirc) 
[![Powered by SymPy](https://img.shields.io/badge/powered%20by-SymPy-Green.svg?style=flat&colorA=gray&colorB=green)](https://www.sympy.org/)

A lightweight python library for symbolic circuit analysis.

## Documentation

Documentation is available on [GitHub Pages](https://matyasvasek.github.io/SymCirc/).

## What can it do?

**SymCirc** currently offers symbolic and semisymbolic **DC**, **AC** and **transient** small frequency circuit analysis.
It supports the following ideal circuit elements: **resistors, inductors, capacitors, independent sources, controlled sources, ideal operational amplifiers and coupled inductors**.
Transient simulation allows for initial conditions of capacitors and standard/coupled inductors.

## Install

**SymCirc** runs on [SymPy](https://www.sympy.org/), which should get automatically installed with SymCirc. Use this command to install via pip
```
pip install symcirc
```
## Build

Or you can build with this command:
```
python setup.py sdist bdist_wheel
```

## Simulation example

``` python
from symcirc import *

# Insert your netlist
netlist = """
CIRCUIT NAME - First line is always the circuit name
# This is a comment
R1 1 0 2k
R2 3 0 (1/G)
V1 2 1 dc 1 ac 1
R3 3 4 6k
C1 3 4 1n
R4 4 0 10k
V2 4 0 dc 5
I1 3 2 dc 1m ac 0
"""

# Execute netlist simulation
analysis_type = "DC"  # or "AC", "TF", "tran"
method = "tableau"  # Default is "two_graph_node", which is faster, but currently lacks coupled inductors.
symbolic = True  # If set to False, only elements which have no numeric value are left as symbolic. In this case only R2 stays symbolic.

circuit = AnalyseCircuit(netlist, analysis_type, symbolic=True, method=method)

all_values = circuit.component_values("all")
latex_formatted_values = to_latex(all_values)

print(all_values)
```

## Netlist syntax
There are many example netlists available in the [netlists](tests/netlists) folder. This folder contains a database of netlists which are used to automatically test every new version of SymCirc.

### Elements
* Resistor:          **RXXX N+ N-** _RESISTANCE_
* Capacitor:         **CXXX N+ N-** _CAPACITY IC=<INIT_VOLTAGE>_
* Inductor:          **LXXX N+ N-** _INDUCTANCE IC=<INIT_CURRENT>_
* Coupling:          **KXXX LYYY LZZZ K**
* Voltage source:    **VXXX N+ N-** _dc VOLTAGE ac AMPLITUDE <PHASE>_
* Current source:    **IXXX N+ N-** _dc CURRENT ac AMPLITUDE <PHASE>_
* IOAMP:             **AXXX Nout1 Nout2 Nin+ Nin-**
* VCVS:              **EXXX N1 N2 Ncontrol+ Ncontrol-** _GAIN_
* CCCS:              **FXXX N1 N2 VSENSE** _GAIN_
* VCCS:              **GXXX N1 N2 Ncontrol+ Ncontrol-** _GAIN_
* CCVS:              **HXXX N1 N2 VSENSE** _GAIN_

### Units
* terra: T
* giga: G
* mega: meg
* kilo: k
* mili: m
* micro: u
* nano: n
* pico: p
* femto: f


