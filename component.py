import sympy


class Component:
    def __init__(self, name=None, type=None, node1=None, node2=None, node3=None, node4=None, sym_value=None,
                 init_cond=None, position=None, value=None, control_voltage=None):
        self.name = name
        self.type = type
        self.value = value
        self.sym_value = sym_value
        self.init_cond = init_cond
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.position = position
        self.control_voltage = control_voltage


