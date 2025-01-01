import sympy


class Component:
    """
    Parent component class

    :param str name: component id
    :param str type: component type id
    :param str node1: first node id
    :param str node2: second node id
    :param sympy_expression sym_value: first node id
    :param sympy_expression value: first node id
    """
    def __init__(self, name=None, type=None, node1=None, node2=None, sym_value=None, value=None):
        self.name = name
        self.type = type
        self.node1 = node1
        self.node2 = node2
        self.sym_value = sym_value
        self.value = value


class Resistor(Component):
    """
        Resistor component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: first node id
        :param sympy_expression value: first node id
    """
    def __init__(self, name, type, node1, node2, sym_value, value):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.netlist_keywords = ["R", "r"]


class Capacitor(Component):
    """
        Capacitor component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: first node id
        :param sympy_expression value: first node id
        :param int init_cond: initial condition
    """
    def __init__(self, name, type, node1, node2, sym_value, value, init_cond=0):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.init_cond = init_cond
        self.netlist_keywords = ["C", "c"]


class Inductor(Component):
    """
        Inductor component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: symbolic value of the component
        :param sympy_expression value: numeric value of the component
        :param int init_cond: initial condition
    """
    def __init__(self, name, type, node1, node2, sym_value, value, init_cond=0, coupling=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.init_cond = init_cond
        self.coupling = coupling
        self.netlist_keywords = ["L", "l"]


class VoltageSource(Component):
    """
        Voltage source component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: symbolic value of the component
        :param sympy_expression dc_value: dc value of the component
        :param sympy_expression ac_value: ac value of the component
        :param sympy_expression tran_value: tran value of the component
        :param int position: this element causes equation matrix expansion and needs the row/col index saved
    """
    def __init__(self, name, type, node1, node2, sym_value, dc_value=None, ac_value=None, tran_value=None,
                 position=None, shorted_nodes=None):
        super().__init__(name, type, node1, node2, sym_value, value=dc_value)
        if shorted_nodes is None:
            self.shorted_nodes = []
        else:
            self.shorted_nodes = shorted_nodes
        self.dc_value = dc_value
        self.ac_value = ac_value
        self.tran_value = tran_value
        self.position = position
        self.netlist_keywords = ["V", "v", "U", "u"]


class CurrentSource(Component):
    """
        Current source component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: symbolic value of the component
        :param sympy_expression dc_value: dc value of the component
        :param sympy_expression ac_value: ac value of the component
        :param sympy_expression tran_value: tran value of the component

    """
    def __init__(self, name, type, node1, node2, sym_value, dc_value=None, ac_value=None, tran_value=None):
        super().__init__(name, type, node1, node2, sym_value, value=dc_value)
        self.dc_value = dc_value
        self.ac_value = ac_value
        self.tran_value = tran_value
        self.netlist_keywords = ["I", "i"]


class OperationalAmplifier(Component):
    """
        Operational amplifier component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param str node2: third node id
        :param str node2: fourth node id
        :param sympy_expression sym_value: symbolic value of the component
        :param int position: this element causes equation matrix expansion and needs the row/col index saved
    """
    def __init__(self, name, type, node1, node2, node3, node4, sym_value, position = None):
        super().__init__(name, type, node1, node2, sym_value)
        self.node3 = node3
        self.node4 = node4
        self.position = position
        self.netlist_keywords = ["A", "a"]


class CurrentControlledSource(Component):
    """
        Current controlled source component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: symbolic value of the component
        :param sympy_expression value: numeric value of the component
        :param str current_sensor: id of the element across which is the controlling current
        :param int position: this element causes equation matrix expansion and needs the row/col index saved

    """
    def __init__(self, name, type, node1, node2, sym_value, current_sensor, position, value=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.current_sensor = current_sensor
        self.position = position
        self.netlist_keywords = ["F", "f", "H", "h"]


class VoltageControlledSource(Component):
    """
        Voltage controlled source component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param str node3: third node id
        :param str node4: fourth node id
        :param sympy_expression sym_value: symbolic value of the component
        :param sympy_expression value: numeric value of the component
        :param int position: this element causes equation matrix expansion and needs the row/col index saved

    """
    def __init__(self, name, type, node1, node2, node3, node4, sym_value, position=None, value=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.node3 = node3
        self.node4 = node4
        self.position = position
        self.netlist_keywords = ["G", "g", "E", "e"]


class Coupling(Component):
    def __init__(self, name, type, L1, L2, sym_value, value):
        super().__init__(name, type, sym_value, value)
        self.L1 = L1
        self.L2 = L2
        self.sym_value = sym_value
        self.value = value
        self.netlist_keywords = ["K", "k"]

class Subcircuit():
    def __init__(self, name, model_id, node_list, param_dict):
        self.name = name
        self.model_id = model_id
        self.node_list = node_list
        self.param_dict = param_dict
        self.netlist_keywords = ["X", "x"]

class SubcktModel():
    def __init__(self, model_id, node_list, param_dict):
        self.model_id = model_id
        self.node_list = node_list
        self.param_dict = param_dict
        self.elements = []

class PeriodicSwitch(Component):
    # periodic switch used for SC/SI analysis
    def __init__(self, name, type, node1, node2, phase):
        super().__init__(name, type, node1, node2)
        self.phase = phase
        self.netlist_keywords = ["S", "s"]

class Short(Component):
    pass


