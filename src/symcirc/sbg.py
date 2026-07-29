from symcirc.component import ParallelAdmittance, SerialAdmittance
from symcirc.analysis import ACNumeric, Circuit
from symcirc.utils import mag, s, j, pi, evalf, xpoints
from symcirc.simplification import build_signalflow_graph
import time


def generate_reference(circuit, node, control_freqs, err_metric):
    # Semisymbolic
    '''
    full_ac_analysis = AC(circuit, method="two_graph_node", symbolic=False)
    vout = evalf(full_ac_analysis.get_node_voltage(node))
    global_ctrl_values = ypoints(vout, control_freqs, f)
    '''
    # Numeric
    ac_numeric = ACNumeric(circuit, method="two_graph_node", symbolic=False)
    results = ac_numeric.run(control_freqs)
    global_ctrl_values = results[f"v({node})"]

    if err_metric == "mag":
        global_ctrl_values = [mag(x, do_evalf=False) for x in global_ctrl_values]

    return global_ctrl_values

def greedy_grow(circuit: Circuit, control_freqs, node, fixed_comps=None, max_err=0.1, err_metric="s_plane_dist", relative_err=False):
    t1 = time.time()
    if fixed_comps is None:
        fixed_comps = []
    global_ctrl_values = generate_reference(circuit, node, control_freqs, err_metric)
    global_effect = sum(global_ctrl_values) / len(global_ctrl_values)
    print(global_effect)

    t_ref_gen = time.time()
    print(f"Reference gen time: {t_ref_gen-t1}")
    sf_graph = build_signalflow_graph(circuit, "Vi")
    t_graph_built = time.time()
    print(f"Graph build time: {t_graph_built - t_ref_gen}")
    #sf_graph.visualize()
    t_graph_visu = time.time()
    print(f"Graph visualisation time: {t_graph_visu - t_graph_built}")

    # TODO: make this general
    t_search_start = time.time()
    path = sf_graph.find_most_impactful_path("nodein_negv", f"{node}v", control_freqs)
    print(path)
    print(f"Time: {time.time()-t_search_start}; path: {path}")
    '''paths = sf_graph.find_all_paths("nodeiv", f"{node}v")
    total_gain = 0
    t_graph_paths = time.time()
    print(f"Path enumeration time: {t_graph_paths - t_graph_visu}")'''
    paths = [path]
    # Remove path components
    path_comps = set()
    for path in paths:
        path_comps.update(path["elements"])
    comps_to_be_removed = path_comps
    for comp in comps_to_be_removed:
        circuit.pop(comp.name)

    # Grow paths
    for path in paths:
        num_gain_list = []
        nodes = path["path"]
        elements = path["elements"]
        num_gain = path["num_gain"]

        # Approximate mean effect
        for freq in control_freqs:
            eval_num_gain = evalf(num_gain, subs={s: 2*pi*freq*j})
            num_gain_list.append(mag(eval_num_gain))
        effect = sum(num_gain_list) / len(num_gain_list)
        path["effect"] = effect

    # Sort by effect strength (highest first)
    paths.sort(key=lambda p: p["effect"], reverse=True)

    for path in paths:
        comps = path["elements"]
        for comp in comps:
            if comp.name not in circuit.components:
                circuit.add(comp)

        err = greedy_prune(circuit, control_freqs, node, fixed_comps, 0.01, err_metric=err_metric,
                           relative_err=True)

        ac_numeric = ACNumeric(circuit, method="two_graph_node", symbolic=False)

        results = ac_numeric.run(control_freqs)
        volt = results[f"v({node})"]

        if err_metric == "mag":
            volt = [mag(x, do_evalf=False) for x in volt]

        if relative_err:
            err = rel_err(global_ctrl_values, volt)
        else:
            err = abs_err(global_ctrl_values, volt)

        print(f"ERR: {err}")

        if err < max_err:
            break

    print(f"Greedy growth time: {time.time()-t1}")
    err = rel_err(global_ctrl_values, volt)
    return err


def greedy_prune(circuit, control_freqs, node, fixed_comps=None, max_err=0.1, err_metric="s_plane_dist", relative_err=False):
    if fixed_comps is None:
        fixed_comps = []
    global_ctrl_values = generate_reference(circuit, node, control_freqs, err_metric)

    comp_name_list = list(circuit.components.keys())

    ''' SORT COMPONENTS '''
    '''order = ["c", "r", "l", "C", "L", "R"]
    priority = {letter: i for i, letter in enumerate(order)}
    comp_name_list = sorted(comp_name_list, key=lambda x: priority.get(x[0], float('inf')))'''

    err = 0
    volt = None
    last_volt = global_ctrl_values
    for comp_name in comp_name_list:
        short = None
        try:
            if comp_name in fixed_comps:
                continue
            comp = circuit.pop(comp_name)

            if comp.type == "v":
                short = comp.short()
                circuit.add(short)

            # Analyze in control freqs
            '''
            partial_ac_analysis = AC(circuit, method="two_graph_node", symbolic=False)
            vout = evalf(partial_ac_analysis.get_node_voltage(node))

            volt = ypoints(vout, control_freqs, f)
            '''

            ''' ------RUN NUMERIC------ '''
            ac_numeric = ACNumeric(circuit, method="two_graph_node", symbolic=False)
            results = ac_numeric.run(control_freqs)
            volt = results[f"v({node})"]

            if err_metric == "mag":
                volt = [mag(x, do_evalf=False) for x in volt]

            if relative_err:
                err = rel_err(global_ctrl_values, volt)
            else:
                err = abs_err(global_ctrl_values, volt)

            if err > max_err:
                circuit.add(comp)
                if short is not None:
                    circuit.pop(short.name)
            else:
                last_volt = volt
        except:
            circuit.add(comp)
            continue

    err = rel_err(global_ctrl_values, last_volt)
    return err

def prune_components(circuit, control_freqs, threshold=0.1):
    prunable_comp_list = []
    for component in circuit.components.values():
        if isinstance(component, (ParallelAdmittance, SerialAdmittance)):
            prunable_comp_list.append(component.name)

    for comp_name in prunable_comp_list:
        c = circuit.pop(comp_name)
        prune_comp(c, control_freqs, threshold)
        circuit.add(c)

def prune_comp(comp, control_freqs, threshold=0.1):
    if isinstance(comp, ParallelAdmittance):
        val = comp.y()
        parallel = True
    elif isinstance(comp, SerialAdmittance):
        val = comp.z()
        parallel = False
    else:
        return
    val_points = []
    for freq in control_freqs:
        val_points.append(evalf(val, {s: j * 2 * pi * freq}))
    comp_name = comp.name
    component_names = list(comp.subcomponents.keys()).copy()
    for subcomp_name in component_names:
        subcomp = comp.subcomponents[subcomp_name]
        if parallel:
            subcomp_val = subcomp.y()
        else:
            subcomp_val = subcomp.z()
        i = 0
        ratio_max = 0
        for freq in control_freqs:
            val_sub_point = evalf(subcomp_val, {s: j * 2 * pi * freq})
            ratio = abs(val_sub_point / val_points[i])
            ratio_max = max(ratio, ratio_max) # TODO: update metric
            i += 1
        if ratio_max < threshold:
            comp.remove(subcomp)

def abs_err(global_ctrl_values, volt):
    res = []
    largest_dist = 0
    for x, y in zip(global_ctrl_values, volt):
        distance = abs(x - y)
        res.append(distance)
        if distance > largest_dist:
            largest_dist = distance
    return largest_dist

def rel_err(global_ctrl_values, volt):
    res = []
    largest_dist = 0
    for x, y in zip(global_ctrl_values, volt):
        if x == 0:
            if y == 0:
                continue
            else:
                distance = 1000000000
        else:
            distance = abs(x - y)/abs(x)
        res.append(distance)
        if distance > largest_dist:
            largest_dist = distance
    return largest_dist


def run_optimization(circuit, node_of_interest, max_err, bw_start, bw_end, ctrl_points, is_err_relative=True, metric="mag", fixed_comps=None, method="greedy_prune", solver="gauss"):
    log = {"title": f"SBG optimization log", "circuit": circuit,
           "node": node_of_interest, "BW": (bw_start, bw_end),
           "points": ctrl_points, "metric": metric,
           "relative_err": is_err_relative, "fixed_comps": fixed_comps,
           "method": method}

    control_freqs = xpoints(bw_start, bw_end, ctrl_points, log=True)

    start_comp_no = 0
    for comp in circuit.components.values():
        if comp.type not in ["w"]:
            start_comp_no += 1

    log["comp_no_start"] = len(circuit.components)

    t_tmp = time.time()
    if method == "greedy_grow":
        err = greedy_grow(circuit, control_freqs, node_of_interest, fixed_comps, max_err, err_metric=metric,
                           relative_err=is_err_relative)
    elif method == "greedy_prune":
        err = greedy_prune(circuit, control_freqs, node_of_interest, fixed_comps, max_err, err_metric=metric,
                           relative_err=is_err_relative)
    else:
        raise(NotImplementedError(f"Unknown circuit model optimization method: {method}; Allowed are: ['greedy_prune', 'greedy_grow']"))
    t_sbg = time.time() - t_tmp

    sbg_comp_no = 0
    for comp in circuit.components.values():
        if comp.type not in ["w"]:
            sbg_comp_no += 1

    t = time.time()
    log["comp_no_after"] = len(circuit.components)
    log["max_err"] = err

    return log

