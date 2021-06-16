node_ref = []
with open('Input_node.inp','r') as f:
    lines = f.readlines() # lines is a list
    for line in lines:
        line = line.strip('\n')
        line = line.split(',')
        node_ref.append(( int(line[0]), line[1], line[2], line[3] ))

node_BC = []
with open('INIT_BC') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip('\n')
        line = line.split()
        node_BC.append(( int(line[0]), line[1] ))

with open('frame_levelset_BC.vtk', 'w') as f:
    print>>f, '# vtk DataFile Version 3.0'
    print>>f, '3D shear band'
    print>>f, 'ASCII'
    print>>f, 'DATASET UNSTRUCTURED_GRID' 
    print>>f, 'POINTS %d float' % len(node_BC)
    for node in node_BC:
        node_num = node[0]
        print>>f, '%8s %8s %8s'  %node_ref[node_num-1][1:]
    print>>f, 'POINT_DATA %d' % len(node_BC)
    print>>f, 'SCALARS levelset float'
    print>>f, 'LOOKUP_TABLE default'
    for node in node_BC:
        print>>f, '%8s' %node[1]
