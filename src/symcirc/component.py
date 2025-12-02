import sympy
from symcirc.utils import s

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
    def __init__(self, name=None, node1=None, node2=None, sym_value=None, value=None):
        self.name = name
        self.node1 = node1
        self.node2 = node2

        if sym_value is None:
            self.sym_value = sympy.Symbol(name)
        else:
            self.sym_value = sym_value

        if value is None:
            self.value = 0
        else:
            self.value = value

    def nodes(self):
        ret = set()
        if self.node1 is not None:
            ret.add(self.node1)
        if self.node2 is not None:
            ret.add(self.node2)
        return ret


class Resistor(Component):
    """
        Resistor component class

        :param str name: component id
        :param str type: component type id
        :param str node1: first node id
        :param str node2: second node id
        :param sympy_expression sym_value: a sympy expression used in symbolic
        :param sympy_expression value: a numeric value used in semisymbolic
    """
    def __init__(self, name, node1, node2, sym_value=None, value=None):
        super().__init__(name, node1, node2, sym_value, value)
        self.netlist_keywords = ["R", "r"]
        self.type = "r"


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
    def __init__(self, name, node1, node2, sym_value=None, value=None, init_cond=0):
        super().__init__(name, node1, node2, sym_value, value)
        self.init_cond = init_cond
        self.netlist_keywords = ["C", "c"]
        self.type = "c"


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
    def __init__(self, name, node1, node2, sym_value=None, value=None, init_cond=0, coupling=None):
        super().__init__(name, node1, node2, sym_value, value)
        self.init_cond = init_cond
        self.coupling = coupling
        self.netlist_keywords = ["L", "l"]
        self.type = "l"


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
    def __init__(self, name, node1, node2,
                 dc_num=0, dc_sym=None,
                 ac_num=0, ac_sym=None, ac_phase=0,
                 tran_num=0, tran_sym=None,
                 position=None, shorted_nodes=None):
        super().__init__(name, node1, node2, sym_value=dc_sym, value=dc_num)
        if shorted_nodes is None:
            self.shorted_nodes = []
        else:
            self.shorted_nodes = shorted_nodes

        symbol = sympy.Symbol(name)

        self.dc_num = dc_num
        if dc_sym is None:
            self.dc_sym = symbol
        else:
            self.dc_sym = dc_sym

        self.ac_num = ac_num
        if ac_sym is None:
            self.ac_sym = symbol
        else:
            self.ac_sym = ac_sym
        self.ac_phase = ac_phase

        self.tran_num = tran_num

        if tran_sym is None:
            self.tran_sym = symbol/s
        else:
            self.tran_sym = tran_sym

        self.position = position
        self.netlist_keywords = ["V", "v", "U", "u"]
        self.type = "v"


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
    def __init__(self, name, node1, node2,
                 dc_num=0, dc_sym=None,
                 ac_num=0, ac_sym=None, ac_phase=0,
                 tran_num=0, tran_sym=None,
                 position=None, shorted_nodes=None):
        super().__init__(name, node1, node2, sym_value = dc_sym, value=dc_num)
        if shorted_nodes is None:
            self.shorted_nodes = []
        else:
            self.shorted_nodes = shorted_nodes
        self.dc_num = dc_num
        self.dc_sym = dc_sym
        self.ac_num = ac_num
        self.ac_sym = ac_sym
        self.ac_phase = ac_phase
        self.tran_num = tran_num
        self.tran_sym = tran_sym
        self.position = position

        self.netlist_keywords = ["I", "i"]
        self.type = "i"


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
    def __init__(self, name, node1, node2, node3, node4, sym_value, position = None):
        super().__init__(name, node1, node2, sym_value)
        self.node3 = node3
        self.node4 = node4
        self.position = position
        self.netlist_keywords = ["A", "a"]
        self.type = "a"

    def nodes(self):
        return {self.node1, self.node2, self.node3, self.node4}


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
        super().__init__(name, node1, node2, sym_value, value)
        self.current_sensor = current_sensor
        self.position = position
        self.type = type
        self.netlist_keywords = ["F", "f", "H", "h"]

    def nodes(self):
        return {self.node1, self.node2, self.node3, self.node4}


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
        super().__init__(name, node1, node2, sym_value, value)
        self.node3 = node3
        self.node4 = node4
        self.position = position
        self.type = type
        self.netlist_keywords = ["G", "g", "E", "e"]

    def nodes(self):
        return {self.node1, self.node2, self.node3, self.node4}


class Coupling(Component):
    def __init__(self, name, L1, L2, sym_value, value):
        super().__init__(name, None, None, sym_value, value)
        self.L1 = L1
        self.L2 = L2
        self.sym_value = sym_value
        self.value = value
        self.netlist_keywords = ["K", "k"]
        self.type = "k"


class Subcircuit():
    def __init__(self, name, model_id, node_list, param_dict):
        self.name = name
        self.model_id = model_id
        self.node_list = node_list
        self.param_dict = param_dict
        self.netlist_keywords = ["X", "x"]
        self.type = "x"


class SubcktModel():
    def __init__(self, model_id, node_list, param_dict):
        self.model_id = model_id
        self.node_list = node_list
        self.param_dict = param_dict
        self.elements = []


class DiodeModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["a", "k"]
        super().__init__(model_id, node_list, param_dict)
        self.elements = []
        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        params = param_dict.keys()
        if "gd" in params:
            self.elements.append(f"rd a k 1/{param_dict['gd']}")
        else:
            self.elements.append("rd a k")


class NPNModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["c", "b", "e"]
        super().__init__(model_id, node_list, param_dict)
        self.elements = []

        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        params = param_dict.keys()

        # Resistors
        if "gpi" in params:
            self.elements.append(f"rpi b e 1/{param_dict['gpi']}")
        else:
            self.elements.append(f"rpi b e 1/gpi")
        if "go" in params:
            self.elements.append(f"ro c e 1/{param_dict['go']}")
        if "gx" in params:
            self.elements.append(f"ro b e 1/{param_dict['gx']}")

        # VCCS
        if "gm" in params:
            self.elements.append(f"gm c e b e {param_dict['gm']}")
        else:
            self.elements.append(f"gm c e b e")

        # Parasitic caps
        if "cjc" in params:
            self.elements.append(f"cmu b c {param_dict['cjc']}")
        if "cje" in params:
            self.elements.append(f"cpi b e {param_dict['cje']}")

        print(self.elements)


class PNPModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["c", "b", "e"]
        super().__init__(model_id, node_list, param_dict)
        self.elements = []
        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        params = param_dict.keys()

        # Resistors
        if "gpi" in params:
            self.elements.append(f"rpi b e 1/{param_dict['gpi']}")
        else:
            self.elements.append(f"rpi b e 1/gpi")
        if "go" in params:
            self.elements.append(f"ro c e 1/{param_dict['go']}")
        if "gx" in params:
            self.elements.append(f"ro b e 1/{param_dict['gx']}")

        # VCCS
        if "gm" in params:
            self.elements.append(f"gm c e e b {param_dict['gm']}")
        else:
            self.elements.append(f"gm c e b e")

        # Parasitic caps
        if "cjc" in params:
            self.elements.append(f"cmu b c {param_dict['cjc']}")
        if "cje" in params:
            self.elements.append(f"cpi b e {param_dict['cje']}")


class NFETModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["d", "g", "s", "b"]
        super().__init__(model_id, node_list, param_dict)
        self.elements = []
        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        params = param_dict.keys()

        # VCCS
        if "gm" in params:
            self.elements.append(f"gm d s g s {param_dict['gm']}")
        else:
            self.elements.append("gm d s g s")

        if "gmb" in params:  # body-effect VCCS
            self.elements.append(f"gmb d s b s {param_dict['gmb']}")

        # Parasitic caps
        if "cbd" in params:  # parasitic cap bulk-drain
            self.elements.append(f"cbd b d {param_dict['cbd']}")
        if "cbs" in params:  # parasitic cap bulk-source
            self.elements.append(f"cbs b s {param_dict['cbs']}")
        if "cgb" in params:  # parasitic cap gate-bulk
            self.elements.append(f"cgb g b {param_dict['cgb']}")
        if "cgs" in params:  # parasitic cap gate-source
            self.elements.append(f"cgs g s {param_dict['cgs']}")
        if "cgd" in params:  # parasitic cap gate-drain
            self.elements.append(f"cgd g d {param_dict['cgd']}")
        if "cds" in params:  # parasitic cap drain-source
            self.elements.append(f"cds d s {param_dict['cds']}")

        # Parasitic res
        if "gds" in params:  # parasitic res drain-source
            self.elements.append(f"rds d s 1/{param_dict['gds']}")
        if "gbs" in params:  # parasitic res bulk-source
            self.elements.append(f"rbs b s 1/{param_dict['gbs']}")
        if "gbd" in params:  # parasitic res bulk-drain
            self.elements.append(f"rbd b d 1/{param_dict['gbd']}")


class PFETModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["d", "g", "s", "b"]
        super().__init__(model_id, node_list, param_dict)
        self.elements = []
        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        params = param_dict.keys()

        # VCCS
        if "gm" in params:  # main VCCS
            self.elements.append(f"gm s d g s {param_dict['gm']}")
        else:
            self.elements.append("gm s d g s")

        if "gmb" in params: # body-effect VCCS
            self.elements.append(f"gmb s d b s {param_dict['gmb']}")

        # Parasitic caps
        if "cbd" in params:  # parasitic cap bulk-drain
            self.elements.append(f"cbd b d {param_dict['cbd']}")
        if "cbs" in params:  # parasitic cap bulk-source
            self.elements.append(f"cbs b s {param_dict['cbs']}")
        if "cgb" in params:  # parasitic cap gate-bulk
            self.elements.append(f"cgb g b {param_dict['cgb']}")
        if "cgs" in params:  # parasitic cap gate-source
            self.elements.append(f"cgs g s {param_dict['cgs']}")
        if "cgd" in params:  # parasitic cap gate-drain
            self.elements.append(f"cgd g d {param_dict['cgd']}")
        if "cds" in params:  # parasitic cap drain-source
            self.elements.append(f"cds d s {param_dict['cds']}")

        # Parasitic res
        if "gds" in params:  # parasitic res drain-source
            self.elements.append(f"rds d s 1/{param_dict['gds']}")
        if "gbs" in params:  # parasitic res bulk-source
            self.elements.append(f"rbs b s 1/{param_dict['gbs']}")
        if "gbd" in params:  # parasitic res bulk-drain
            self.elements.append(f"rbd b d 1/{param_dict['gbd']}")


class PeriodicSwitch(Component):
    # periodic switch used for SC/SI analysis
    def __init__(self, name, type, node1, node2, phase):
        super().__init__(name, type, node1, node2)
        self.phase = phase
        self.netlist_keywords = ["S", "s"]


class Short(Component):
    pass


