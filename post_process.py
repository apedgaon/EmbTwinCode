import os
import glob
import numpy as np

# Load mesh
elem_conn2 = np.loadtxt('input_elmt_1000.inp', delimiter=',', dtype=int)
elem_conn = elem_conn2[:, 1:]
elem_conn = elem_conn - 1

nodal_coords = np.loadtxt('input_node_1000.inp', delimiter=',')
nodal_coords = nodal_coords[:, 1:]

# Load time and create its dictionary
time = np.loadtxt('kinc_time.dat')
num_frames = len(time)
time_dict = dict()
for i in np.arange(num_frames):
    time_dict[int(time[i, 0])] = time[i, 1]

# Loop over vuel element status files and generate vtk for all results
directory = os.getcwd()
for file in glob.glob(directory + os.sep + "kinc_elmt_stat_*.dat"):
    file_name_components = file.split('_')
    sval = file_name_components[-1].split('.')[0]

    # Open vtk file with increment number and start writing
    output_file = os.path.join(directory, 'twin_' + sval + '.vtu')
    if os.path.exists(output_file):
        continue

    print(sval)
    vtk_stream = open(output_file, "w", newline="\n")

    # Header
    vtk_stream.write("<?xml version=\"1.0\"?>\n")
    vtk_stream.write("<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n")
    vtk_stream.write("<UnstructuredGrid>\n")

    # Write time for this frame
    vtk_stream.write("<FieldData>\n")
    vtk_stream.write("<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">\n")
    vtk_stream.write(str(time_dict[int(sval)]) + "\n")
    vtk_stream.write("</DataArray>\n")
    vtk_stream.write("</FieldData>\n")

    # Write mesh info
    num_elems = elem_conn.shape[0]
    num_nodes = nodal_coords.shape[0]
    dim_num = nodal_coords.shape[1]
    vtk_stream.write("<Piece  NumberOfPoints=\"" + str(num_nodes) + "\"  NumberOfCells=\"" + str(num_elems) + "\">\n")

    # Write nodal coords
    vtk_stream.write("<Points>\n")
    vtk_stream.write("<DataArray  type=\"Float32\"  NumberOfComponents=\"" + str(dim_num) + "\"  format=\"ascii\">\n")
    for coords in nodal_coords:
        for coord in coords:
            vtk_stream.write(str(coord) + "\t")

        vtk_stream.write("\n")

    vtk_stream.write("</DataArray>\n")
    vtk_stream.write("</Points>\n")     # End of mesh nodal data

    # Write element connectivity, offsets and type
    vtk_stream.write("<Cells>\n")
    vtk_stream.write("<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">\n")
    for conn in elem_conn:
        for node_id in conn:
            vtk_stream.write(str(node_id) + "\t")
        
        vtk_stream.write("\n")
    vtk_stream.write("</DataArray>\n")

    # Offsets
    vtk_stream.write("<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">\n")
    offset = 0
    for conn in elem_conn:
        offset = offset + conn.shape[0]
        vtk_stream.write(str(offset) + "\n")

    vtk_stream.write("</DataArray>\n")

    # Element Types
    vtk_stream.write("<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n")
    for conn in elem_conn:
        elem_order = conn.shape[0]
        if ((dim_num == 2) and (elem_order == 3)):
            elem_type = 5
        elif ((dim_num == 2) and (elem_order == 4)):
            elem_type = 9
        elif ((dim_num == 2) and (elem_order == 9)):
            elem_type = 28
        elif ((dim_num == 3) and (elem_order == 4)):
            elem_type = 10
        elif ((dim_num == 3) and (elem_order == 8)):
            elem_type = 12
        else:
            raise Exception("Unknown Element Type!")
        
        vtk_stream.write(str(elem_type) + "\n")
    
    vtk_stream.write("</DataArray>\n")
    vtk_stream.write("</Cells>\n")  # End of mesh element data
    
    # Start writing results

    # Start writing elemental results
    vtk_stream.write("<CellData  Scalars=\"u\">\n")

    # Load Element status
    elem_stat = np.loadtxt(file, delimiter=',')
    scalar_name = 'Element Status'

    # Write element status
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(elem_stat[i]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Load stress data
    input_file = os.path.join(directory, 'kinc_elmt_stress_' + sval + '.dat')
    sdata  = np.loadtxt(input_file)

    # Write Von Mises Stress in the untwinned region
    stress = sdata.reshape(num_elems, 12)
    stress_vm = np.sqrt(1.0/2 * ((stress[:,0] - stress[:,1])**2 + (stress[:,1] - stress[:,2])**2 + (stress[:,2] - stress[:,0])**2)
                    + 3 * (stress[:,3]**2 + stress[:,4]**2 + stress[:,5]**2))
    scalar_name = 'Von Mises Stress U'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(stress_vm[i]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write Von Mises Stress in the twinned region
    stress_vm = np.sqrt(1.0/2 * ((stress[:,6] - stress[:,7])**2 + (stress[:,7] - stress[:,8])**2 + (stress[:,8] - stress[:,6])**2)
                    + 3 * (stress[:,9]**2 + stress[:,10]**2 + stress[:,11]**2))
    scalar_name = 'Von Mises Stress Z'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(stress_vm[i]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Load SDV
    input_file = os.path.join(directory, 'kinc_elmt_SDV_' + sval + '.dat')
    SDV = np.loadtxt(input_file)

    # Write temperature in the untwinned region
    scalar_name = 'Temperature U'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 0]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write cumulative dislocation density in the untwinned region
    scalar_name = 'Dislocation Density U'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 4]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write effective temperature in the untwinned region
    scalar_name = 'Effective Temperature U'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 5]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Wrtie misorientation angle in the untwinned region
    scalar_name = 'Misorientation Angle U'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 6]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write temperature in the twinned region
    scalar_name = 'Temperature Z'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 7]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write cumulative dislocation densitiy in the twinned region
    scalar_name = 'Dislocation Density Z'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 11]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write effective temperature in the twinned region
    scalar_name = 'Effective Temperature Z'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 12]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Write misorientation angle in the twinned region
    scalar_name = 'Misorientation Angle Z'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(SDV[i, 13]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    # Load twin nucleation probability
    input_file = os.path.join(directory, 'kinc_nuc_prob_' + sval + '.dat')
    nuc_prob = np.loadtxt(input_file)

    # Write twin nucleation probability
    scalar_name = 'Nucleation Probability'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + scalar_name + 
                     "\"  NumberOfComponents=\"" + str(1) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_elems):
        vtk_stream.write(str(nuc_prob[i]) + "\n")
    
    vtk_stream.write("</DataArray>\n")

    vtk_stream.write("</CellData> \n")  # End of elemental results

    # Start writing the nodal results
    vtk_stream.write("<PointData  Scalars=\"u\">\n")

    # Load displacements
    input_file = os.path.join(directory, 'kinc_disp_' + sval + '.dat')
    disp = np.loadtxt(input_file)

    # Write displacements
    vector_name = 'Displacement'
    vtk_stream.write("<DataArray  type=\"Float32\"  Name=\"" + vector_name +  "\"  NumberOfComponents=\"" + str(dim_num) + "\"  format=\"ascii\">\n")
    for i in np.arange(num_nodes):
        for j in np.arange(dim_num):
            vtk_stream.write(str(disp[i, j]) + "\t")
        
        vtk_stream.write("\n")

    vtk_stream.write("</DataArray>\n")

    vtk_stream.write("</PointData> \n")     # End of nodal results

    # VTK file footer
    vtk_stream.write("</Piece> \n")
    vtk_stream.write("</UnstructuredGrid> \n")
    vtk_stream.write("</VTKFile>")
    vtk_stream.close()