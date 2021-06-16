import sys, glob
## Read nodal coordinates in the reference configuration
target_elmt = int(sys.argv[1])

node_ref = []
with open('Input_node.inp','r') as f:
    lines = f.readlines() # lines is a list
    for line in lines:
        line = line.strip('\n')
        line = line.split(',')
        node_ref.append(( int(line[0]), line[1], line[2], line[3] ))

node_num = len(node_ref)

## Read elmt connectivity from the input file
elmts_VU1 = []
with open('Input_elmt.inp','r') as f:
    lines = f.readlines() # lines is a list
    for line in lines:
        line = line.strip('\n')
        line = line.split(',')
        elmts_VU1.append(( int(line[0]), int(line[1]), int(line[2]), 
                           int(line[3]), int(line[4]), int(line[5]),
                           int(line[6]), int(line[7]), int(line[8]) ))

elmts_num = len(elmts_VU1)

## Organize the elmt number into ascending order
elmts_VU1_ascend = range(elmts_num)
for elmt in elmts_VU1:
    label = elmt[0]
    elmts_VU1_ascend[label-1] = elmt[1:]

elmts_VU1 = elmts_VU1_ascend


elmt_sb = range(elmts_num) # elmt shear band status
elmt_SDV = range(elmts_num) # elmt shear band status
node_coords = range(node_num)
node_theta = range(node_num)
target_elmt_SDV = []
for file in glob.glob("kinc_*_elmts.dat"):
    kinc_num = file.strip('kinc_').strip('_elmts.dat') 

    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            line = line.strip()
            line = line.split()
            label = int(line[0])
            elmt_sb[label-1] = line

    coords_file = 'kinc_' + kinc_num + '_coords.dat'
    with open(coords_file,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            line = line.strip()
            line = line.split()
            label = int(line[0])
            node_coords[label-1] = line

    theta_file = 'kinc_' + kinc_num + '_theta.dat'
    with open(theta_file,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            line = line.strip()
            line = line.split()
            label = int(line[0])
            node_theta[label-1] = line

    SDV_file = 'kinc_' + kinc_num + '_SDV.dat'
    with open(SDV_file,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            line = line.strip()
            line = line.split()
            label = int(line[0])
            elmt_SDV[label-1] = line

    if int(kinc_num) < 10:
        kinc_num = '00000'+ kinc_num
    elif int(kinc_num) < 100:
        kinc_num = '0000'+ kinc_num
    elif int(kinc_num) < 1000:
        kinc_num = '000'+ kinc_num
    elif int(kinc_num) < 10000:
        kinc_num = '00'+ kinc_num
    else:
        kinc_num = '0'+ kinc_num

    file = 'frame' + kinc_num + '.vtk'
    with open(file, 'w') as f:
        print>>f, '# vtk DataFile Version 3.0'
        print>>f, '3D shear band'
        print>>f, 'ASCII'
        print>>f, 'DATASET UNSTRUCTURED_GRID' 
        print>>f, 'POINTS %d float' % node_num
        for i in range(node_num):
            print>>f, '%8s %8s %8s' \
                   %(node_coords[i][1],
                     node_coords[i][2],
                     node_coords[i][3])

        print>>f, 'CELLS %d %d' \
            %(elmts_num, 9*elmts_num)
        for i in range(elmts_num):
            print>>f, '%6d %6d %6d %6d %6d %6d %6d %6d %6d' \
                %(8, elmts_VU1[i][0]-1, elmts_VU1[i][1]-1, elmts_VU1[i][2]-1,
                     elmts_VU1[i][3]-1, elmts_VU1[i][4]-1, elmts_VU1[i][5]-1,
                     elmts_VU1[i][6]-1, elmts_VU1[i][7]-1)

        print>>f, 'CELL_TYPES %d' %elmts_num
        for i in range(elmts_num):
            print>>f, '12'

        print>>f, 'CELL_DATA %d' %elmts_num
        print>>f, 'FIELD SDVn 11'

        print>>f, 'Shearband 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%1s' %elmt_sb[i][1]

        print>>f, 'SDV1 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][1]

        print>>f, 'SDV2 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][2]

        print>>f, 'SDV3 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][3]

        print>>f, 'SDV4 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][4]

        print>>f, 'SDV5 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][5]

        print>>f, 'SDV6 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][6]

        print>>f, 'SDV7 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][7]

        print>>f, 'SDV8 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][8]

        print>>f, 'SDV9 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][9]

        print>>f, 'SDV10 1 %d float' %elmts_num
        for i in range(elmts_num):
            print>>f, '%8s' %elmt_SDV[i][10]

        print>>f, 'POINT_DATA %d' %node_num
        print>>f, 'SCALARS levelset float'
        print>>f, 'LOOKUP_TABLE default'
        for i in range(node_num):
            print>>f, '%8s' %node_theta[i][1]

    target_elmt_SDV.append(( int(kinc_num), elmt_SDV[(target_elmt)-1][1],
                                    elmt_SDV[(target_elmt)-1][2],
                                    elmt_SDV[(target_elmt)-1][3],
                                    elmt_SDV[(target_elmt)-1][4],
                                    elmt_SDV[(target_elmt)-1][5],
                                    elmt_SDV[(target_elmt)-1][6],
                                    elmt_SDV[(target_elmt)-1][7],
                                    elmt_SDV[(target_elmt)-1][8],
                                    elmt_SDV[(target_elmt)-1][9],
                                    elmt_SDV[(target_elmt)-1][10]   ))


## sort elmt_xx.hist according to increments
sorted_elmt_SDV = sorted(target_elmt_SDV, key=lambda a: a[0])
target_elmt_file = 'elmt_' + str(target_elmt) + '.hist'
with open(target_elmt_file, 'w') as f:
    for line in sorted_elmt_SDV:
        print>>f, '%8d %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s' \
                   %line


## write vtk file for mesh
with open('frame_mesh'+'.vtk', 'w') as f:
    print>>f, '# vtk DataFile Version 3.0'
    print>>f, '3D mesh'
    print>>f, 'ASCII'
    print>>f, 'DATASET UNSTRUCTURED_GRID'
    print>>f, 'POINTS %d float' % node_num
    for i in range(node_num):
        print>>f, '%8s %8s %8s' \
               %(node_ref[i][1],
                 node_ref[i][2],
                 node_ref[i][3])

    print>>f, 'CELLS %d %d' \
            %(elmts_num, 9*elmts_num)
    for i in range(elmts_num):
        print>>f, '%6d %6d %6d %6d %6d %6d %6d %6d %6d' \
                %(8, elmts_VU1[i][0]-1, elmts_VU1[i][1]-1, elmts_VU1[i][2]-1,
                     elmts_VU1[i][3]-1, elmts_VU1[i][4]-1, elmts_VU1[i][5]-1,
                     elmts_VU1[i][6]-1, elmts_VU1[i][7]-1)

    print>>f, 'CELL_TYPES %d' %elmts_num
    for i in range(elmts_num):
        print>>f, '12'
