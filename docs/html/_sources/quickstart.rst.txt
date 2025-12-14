Quickstart
==========

Use the following command to install **SymCirc** via ``pip``:

.. code-block:: console

   pip install symcirc


.. code-block:: python

    from symcirc import *

    # Insert your netlist
    netlist = """CIRCUIT NAME - First line is always the circuit name
    * This is a comment
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

    dc_analysis = AnalyseCircuit(netlist, analysis_type, symbolic=True, method=method)

    all_values = dc_analysis.all_component_values()

    for key in all_values:
        print(f"{key}: {all_values[key]}")

    latex_formatted_values = to_latex(all_values)   # Format results into latex


