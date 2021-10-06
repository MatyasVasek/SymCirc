import sympy


class Component:
    def __init__(self, name=None, type=None, node1=None, node2=None, sym_value=None, value=None):
        self.name = name
        self.type = type
        self.node1 = node1
        self.node2 = node2
        self.sym_value = sym_value
        self.value = value


class Resistor(Component):
    def __init__(self, name, type, node1, node2, sym_value, value):
        super().__init__(name, type, node1, node2, sym_value, value)


class Capacitor(Component):
    def __init__(self, name, type, node1, node2, sym_value, value, init_cond=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.init_cond = init_cond


class Inductor(Component):
    def __init__(self, name, type, node1, node2, sym_value, value, init_cond=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.init_cond = init_cond


class VoltageSource(Component):
    def __init__(self, name, type, node1, node2, sym_value, dc_value=None, ac_value=None, position=None):
        super().__init__(name, type, node1, node2, sym_value, value=dc_value)
        self.dc_value = dc_value
        self.ac_value = ac_value
        self.position = position


class CurrentSource(Component):
    def __init__(self, name, type, node1, node2, sym_value, dc_value=None, ac_value=None):
        super().__init__(name, type, node1, node2, sym_value, value=dc_value)
        self.dc_value = dc_value
        self.ac_value = ac_value


class OperationalAmplifier(Component):
    def __init__(self, name, type, node1, node2, node3, node4, sym_value, position):
        super().__init__(name, type, node1, node2, sym_value)
        self.node3 = node3
        self.node4 = node4
        self.position = position


class CurrentControlledSource(Component):
    def __init__(self, name, type, node1, node2, sym_value, control_voltage, position, value=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.control_voltage = control_voltage
        self.position = position


class VoltageControlledSource(Component):
    def __init__(self, name, type, node1, node2, node3, node4, sym_value, position=None, value=None):
        super().__init__(name, type, node1, node2, sym_value, value)
        self.node3 = node3
        self.node4 = node4
        self.position = position
