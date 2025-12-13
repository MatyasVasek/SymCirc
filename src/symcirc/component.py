from typing import List, Set, Dict, Union, Any
import sympy
from symcirc.utils import s


class Component:
    """
    Parent component class

    :param name: component id
    :param node1: first node id
    :param node2: second node id
    :param sym_value: symbolic value
    :param value: numeric value
    """
    def __init__(self, name:str=None, node1:str=None, node2:str=None, sym_value:sympy.Symbol|sympy.Expr=None, value:float=None):
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

    def nodes(self) -> Set[str]:
        """
        :return: set of component nodes
        """
        ret = set()
        if self.node1 is not None:
            ret.add(self.node1)
        if self.node2 is not None:
            ret.add(self.node2)
        return ret


class Coupling(Component):
    def __init__(self, name, L1, L2, sym_value, value):
        super().__init__(name, None, None, sym_value, value)
        self.L1 = L1
        self.L2 = L2
        self.sym_value = sym_value
        self.value = value
        self.netlist_keywords = ["K", "k"]
        self.type = "k"

class Resistor(Component):
    """
        Resistor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
    """
    def __init__(self, name:str, node1:str, node2:str, sym_value:sympy.Symbol|sympy.Expr=None, value:float=None):
        super().__init__(name, node1, node2, sym_value, value)
        self.netlist_keywords = ["R", "r"]
        self.type = "r"


class Capacitor(Component):
    """
        Capacitor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
        :param init_cond: initial condition in volts
    """
    def __init__(self, name:str, node1:str, node2:str, sym_value:sympy.Symbol|sympy.Expr=None, value:float=None,
                 init_cond:sympy.Symbol|sympy.Expr=0):
        super().__init__(name, node1, node2, sym_value, value)
        self.init_cond = init_cond
        self.netlist_keywords = ["C", "c"]
        self.type = "c"


class Inductor(Component):
    """
        Inductor component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param sym_value: symbolic value
        :param value: numeric value used in semisymbolic analyes
        :param init_cond: initial condition in amps
        :param coupling: a coupling object defining the couple
    """
    def __init__(self, name:str, node1:str, node2:str, sym_value:sympy.Symbol|sympy.Expr=None, value:float=None,
                 init_cond:sympy.Symbol|sympy.Expr=0, coupling:Coupling|None=None):
        super().__init__(name, node1, node2, sym_value, value)
        self.init_cond = init_cond
        self.coupling = coupling
        self.netlist_keywords = ["L", "l"]
        self.type = "l"


class VoltageSource(Component):
    """
        Independent voltage source component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param dc_num: numeric dc value used in semisymbolic dc analyes
        :param dc_sym: symbolic dc value
        :param ac_num: numeric ac value used in semisymbolic ac analysis
        :param ac_sym: symbolic ac value
        :param ac_phase: ac phase value
        :param tran_num: numeric transient value used in semisymbolic transient analysis
        :param tran_sym: symbolic transient value
    """
    def __init__(self, name:str, node1:str, node2:str,
                 dc_num: float=0, dc_sym: sympy.Symbol|sympy.Expr=None,
                 ac_num: float=0, ac_sym: sympy.Symbol|sympy.Expr=None, ac_phase: sympy.Symbol|sympy.Expr=0,
                 tran_num: float=0, tran_sym: sympy.Symbol|sympy.Expr=None):
        super().__init__(name, node1, node2, sym_value=dc_sym, value=dc_num)

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

        self.netlist_keywords = ["V", "v", "U", "u"]
        self.type = "v"


class CurrentSource(Component):
    """
        Independent current source component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param dc_num: numeric dc value used in semisymbolic dc analyes
        :param dc_sym: symbolic dc value
        :param ac_num: numeric ac value used in semisymbolic ac analysis
        :param ac_sym: symbolic ac value
        :param ac_phase: ac phase value
        :param tran_num: numeric transient value used in semisymbolic transient analysis
        :param tran_sym: symbolic transient value
    """
    def __init__(self, name:str, node1:str, node2:str,
                 dc_num: float=0, dc_sym: sympy.Symbol|sympy.Expr=None,
                 ac_num: float=0, ac_sym: sympy.Symbol|sympy.Expr=None, ac_phase: float=0,
                 tran_num: float=0, tran_sym: sympy.Symbol|sympy.Expr=None):

        super().__init__(name, node1, node2, sym_value = dc_sym, value=dc_num)
        self.dc_num = dc_num
        self.dc_sym = dc_sym
        self.ac_num = ac_num
        self.ac_sym = ac_sym
        self.ac_phase = ac_phase
        self.tran_num = tran_num
        self.tran_sym = tran_sym

        self.netlist_keywords = ["I", "i"]
        self.type = "i"


class IdealOperationalAmplifier(Component):
    """
        Operational amplifier component class

        :param name: component id
        :param node1: first node id
        :param node2: second node id
        :param node3: third node id
        :param node4: fourth node id
    """
    def __init__(self, name:str, node1:str, node2:str, node3:str, node4:str):
        super().__init__(name, node1, node2)
        self.node3 = node3
        self.node4 = node4
        self.netlist_keywords = ["A", "a"]
        self.type = "a"

    def nodes(self) -> Set[str]:
        """
        :return: set of component nodes
        """
        return {self.node1, self.node2, self.node3, self.node4}


class CurrentControlledSource(Component):
    """
        Current controlled source component class

        :param name: component id
        :param type: component type id
        :param node1: first node id
        :param node2: second node id
        :param current_sensor: id of the element across which is the controlling current
        :param sym_value: symbolic value of the component
        :param value: numeric value of the component
    """
    def __init__(self, name:str, type:str, node1:str, node2:str, current_sensor:str,
                 sym_value:sympy.Symbol|sympy.Expr, value: float=None):
        super().__init__(name, node1, node2, sym_value, value)
        self.current_sensor = current_sensor
        self.type = type
        self.netlist_keywords = ["F", "f", "H", "h"]

    def nodes(self):
        return {self.node1, self.node2, self.node3, self.node4}


class VoltageControlledSource(Component):
    """
        Voltage controlled source component class

        :param name: component id
        :param type: component type id
        :param node1: first node id
        :param node2: second node id
        :param node3: third node id
        :param node4: fourth node id
        :param sym_value: symbolic value of the component
        :param value: numeric value of the component
    """
    def __init__(self, name:str, type:str, node1:str, node2:str, node3:str, node4:str,
                 sym_value:sympy.Symbol|sympy.Expr, value:float=None):
        super().__init__(name, node1, node2, sym_value, value)
        self.node3 = node3
        self.node4 = node4
        self.type = type
        self.netlist_keywords = ["G", "g", "E", "e"]

    def nodes(self):
        return {self.node1, self.node2, self.node3, self.node4}


class SubcktModel:
    """
        Voltage controlled source component class

        :param model_id: component id
        :param node_list: subcircuit interface node list
        :param param_dict: subcircuit parameters dictionary
    """
    def __init__(self, model_id:str, node_list:List[str], param_dict:Dict[str, sympy.Symbol|sympy.Expr]):
        self.model_id = model_id
        self.node_list = node_list
        self.param_dict = self._sanitize_param_dict(param_dict)
        self.elements = []

    def build_instance(self, operational_point: Any) -> List[str]:
        return self.elements

    @staticmethod
    def _sanitize_param_dict(param_dict):
        sanitized_param_dict = {}
        for key in param_dict:
            sanitized_param_dict[key.lower()] = param_dict[key]
        param_dict = sanitized_param_dict
        return param_dict


class DiodeModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["a", "k"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: dict[str: Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        if "gd" in params:
            elements.append(f"rd a k 1/{param_dict['gd']}")
        else:
            elements.append("rd a k")

        return elements


class NPNModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["c", "b", "e"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op:dict[str: Any]=None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # Resistors
        if "gpi" in params:
            elements.append(f"rpi b e 1/{param_dict['gpi']}")
        else:
            elements.append(f"rpi b e")
        if "go" in params:
            elements.append(f"ro c e 1/{param_dict['go']}")
        elif "vaf" in params:
            elements.append(f"ro c e")

        if "gx" in params:
            elements.append(f"rx b e 1/{param_dict['gx']}")

        # VCCS
        if "gm" in params:
            elements.append(f"gm c e b e {param_dict['gm']}")
        else:
            elements.append(f"gm c e b e")

        # Parasitic caps
        if "cjc" in params:
            elements.append(f"cmu b c {param_dict['cjc']}")
        if "cje" in params:
            elements.append(f"cpi b e {param_dict['cje']}")

        return elements


class PNPModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["c", "b", "e"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: dict[str: Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # Resistors
        if "gpi" in params:
            elements.append(f"rpi b e 1/{param_dict['gpi']}")
        else:
            elements.append(f"rpi b e")
        if "go" in params:
            elements.append(f"ro c e 1/{param_dict['go']}")
        elif "vaf" in params:
            elements.append(f"ro c e")
        if "gx" in params:
            elements.append(f"rx b e 1/{param_dict['gx']}")

        # VCCS
        if "gm" in params:
            elements.append(f"gm c e e b {param_dict['gm']}")
        else:
            elements.append(f"gm c e e b")

        # Parasitic caps
        if "cjc" in params:
            elements.append(f"cmu b c {param_dict['cjc']}")
        if "cje" in params:
            elements.append(f"cpi b e {param_dict['cje']}")

        return elements


class NFETModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["d", "g", "s", "b"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: dict[str: Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # VCCS
        if "gm" in params:
            elements.append(f"gm d s g s {param_dict['gm']}")
        else:
            elements.append("gm d s g s")

        if "gmb" in params:  # body-effect VCCS
            elements.append(f"gmb d s b s {param_dict['gmb']}")

        # Parasitic caps
        if "cbd" in params:  # parasitic cap bulk-drain
            elements.append(f"cbd b d {param_dict['cbd']}")
        if "cbs" in params:  # parasitic cap bulk-source
            elements.append(f"cbs b s {param_dict['cbs']}")
        if "cgb" in params:  # parasitic cap gate-bulk
            elements.append(f"cgb g b {param_dict['cgb']}")
        if "cgs" in params:  # parasitic cap gate-source
            elements.append(f"cgs g s {param_dict['cgs']}")
        if "cgd" in params:  # parasitic cap gate-drain
            elements.append(f"cgd g d {param_dict['cgd']}")
        if "cds" in params:  # parasitic cap drain-source
            elements.append(f"cds d s {param_dict['cds']}")

        # Parasitic res
        if "gds" in params:  # parasitic res drain-source
            elements.append(f"rds d s 1/{param_dict['gds']}")
        if "gbs" in params:  # parasitic res bulk-source
            elements.append(f"rbs b s 1/{param_dict['gbs']}")
        if "gbd" in params:  # parasitic res bulk-drain
            elements.append(f"rbd b d 1/{param_dict['gbd']}")

        return elements



class PFETModelAC(SubcktModel):
    def __init__(self, model_id, param_dict):
        node_list = ["d", "g", "s", "b"]
        super().__init__(model_id, node_list, param_dict)

    def build_instance(self, op: dict[str: Any] = None):
        elements = []
        if op is not None:
            param_dict = {**self.param_dict, **self._sanitize_param_dict(op)}
        else:
            param_dict = self.param_dict
        params = param_dict.keys()

        # VCCS
        if "gm" in params:  # main VCCS
            elements.append(f"gm d s g s {param_dict['gm']}")
        else:
            elements.append("gm d s g s")

        if "gmb" in params: # body-effect VCCS
            elements.append(f"gmb d s b s {param_dict['gmb']}")

        # Parasitic caps
        if "cbd" in params:  # parasitic cap bulk-drain
            elements.append(f"cbd b d {param_dict['cbd']}")
        if "cbs" in params:  # parasitic cap bulk-source
            elements.append(f"cbs b s {param_dict['cbs']}")
        if "cgb" in params:  # parasitic cap gate-bulk
            elements.append(f"cgb g b {param_dict['cgb']}")
        if "cgs" in params:  # parasitic cap gate-source
            elements.append(f"cgs g s {param_dict['cgs']}")
        if "cgd" in params:  # parasitic cap gate-drain
            elements.append(f"cgd g d {param_dict['cgd']}")
        if "cds" in params:  # parasitic cap drain-source
            elements.append(f"cds d s {param_dict['cds']}")

        # Parasitic res
        if "gds" in params:  # parasitic res drain-source
            elements.append(f"rds d s 1/{param_dict['gds']}")
        if "gbs" in params:  # parasitic res bulk-source
            elements.append(f"rbs b s 1/{param_dict['gbs']}")
        if "gbd" in params:  # parasitic res bulk-drain
            elements.append(f"rbd b d 1/{param_dict['gbd']}")

        return elements

class PeriodicSwitch(Component):
    # periodic switch used for SC/SI analysis
    def __init__(self, name, type, node1, node2, phase):
        super().__init__(name, type, node1, node2)
        self.phase = phase
        self.netlist_keywords = ["S", "s"]


class Short(Component):
    pass


