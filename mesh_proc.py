import numpy as np

nodal_coords = np.loadtxt("input_node_2.inp", delimiter=',')
nodal_coords = nodal_coords[:, 1:]
num_nodes = nodal_coords.shape[0]

elem_conn = np.loadtxt("input_elmt_2.inp", delimiter=',', dtype=int)
elem_conn = elem_conn[:, 1:]
elem_conn = elem_conn - 1
num_elems = elem_conn.shape[0]

top_nodes = []
bottom_nodes = []
top_centre_node = -1
bottom_centre_node = -1
bottom_rc_node = -1
for ind in np.arange(num_nodes):
    nodal_coord = nodal_coords[ind, :]
    if abs(nodal_coord[1] - 1.5) < 1e-4:    # top nodes
        if abs(nodal_coord[0] - 0) < 1e-4 and abs(nodal_coord[2] - 1.5) < 1e-4:   # top centre node
            top_centre_node = ind + 1
        else:
            top_nodes.append(ind + 1)
    elif abs(nodal_coord[1] + 1.5) < 1e-4:  # bottom nodes
        bottom_nodes.append(ind + 1)
        if (bottom_centre_node != -1) and (bottom_rc_node == -1):   # bottom right centre node
            bottom_rc_node = ind + 1
        if abs(nodal_coord[0] - 0) < 1e-4 and abs(nodal_coord[2] - 1.5) < 1e-4:   # bottom centre node
            bottom_centre_node = ind + 1

print(top_nodes)
print(top_centre_node)
print(bottom_nodes)
print(bottom_centre_node)
print(bottom_rc_node)
print(nodal_coords[bottom_rc_node - 1])
print(nodal_coords[bottom_centre_node - 1])
print(nodal_coords[top_centre_node - 1])