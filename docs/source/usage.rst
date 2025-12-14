Usage
=====

.. _installation:

Install
^^^^^^^

Use the following command to install **SymCirc** via ``pip``:

.. code-block:: console

   pip install symcirc


Hard dependencies
^^^^^^^^^^^^^^^^^

SymPy
-----

**SymCirc** is a lightweight package. It requires only the
`SymPy <https://www.sympy.org/>`_ library, which is used for symbolic
computations.

When installing **SymCirc** via ``pip``, **SymPy** is installed
automatically as a dependency.


Optional dependencies
^^^^^^^^^^^^^^^^^^^^^

gmpy2
-----

To achieve better performance, install the
`gmpy2 <https://gmpy2.readthedocs.io/en/latest/>`_ package.

If **gmpy2** is available, **SymPy** automatically uses it for integer
operations. This significantly improves the performance of
semi-symbolic analysis.

symengine
---------

`SymEngine <https://symengine.org>`_ is an optional symbolic backend for
**SymPy**.

To enable **SymEngine**, either:

* Run your script with the environment variable:

  .. code-block:: console

     USE_SYMENGINE=1

* Or pass the optional argument ``use_symengine=True`` to
  ``AnalyseCircuit``:

  .. code-block:: python

     symcirc.AnalyseCircuit(netlist, use_symengine=True)


numpy and matplotlib
--------------------

These packages are required for Bode plots and other graphing utilities.

