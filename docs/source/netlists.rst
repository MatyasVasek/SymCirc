Netlists
^^^^^^^^

There are many example netlists available in the
`netlists <https://github.com/MatyasVasek/SymCirc/tree/main/tests/netlists>`_ folder. This folder contains a database of
netlists which are used to automatically test every new version of
**SymCirc**.

Elements
--------

* **Resistor**:
  ``RXXX N+ N-`` *RESISTANCE*

* **Capacitor**:
  ``CXXX N+ N-`` *CAPACITY* ``IC=<INIT_VOLTAGE>``

* **Inductor**:
  ``LXXX N+ N-`` *INDUCTANCE* ``IC=<INIT_CURRENT>``

* **Coupling**:
  ``KXXX LYYY LZZZ K``

* **Voltage source**:
  ``VXXX N+ N-`` *dc VOLTAGE ac AMPLITUDE <PHASE>*

* **Current source**:
  ``IXXX N+ N-`` *dc CURRENT ac AMPLITUDE <PHASE>*

* **IOAMP**:
  ``AXXX Nout1 Nout2 Nin+ Nin-``

* **VCVS**:
  ``EXXX N1 N2 Ncontrol+ Ncontrol-`` *GAIN*

* **CCCS**:
  ``FXXX N1 N2 VSENSE`` *GAIN*

* **VCCS**:
  ``GXXX N1 N2 Ncontrol+ Ncontrol-`` *GAIN*

* **CCVS**:
  ``HXXX N1 N2 VSENSE`` *GAIN*

Units
-----

* terra: ``T``
* giga: ``G``
* mega: ``meg``
* kilo: ``k``
* mili: ``m``
* micro: ``u``
* nano: ``n``
* pico: ``p``
* femto: ``f``
