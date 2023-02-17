# Call library odbAccess which gives 
# tools to open odb files (must be opened from Abaqus)
# Call numpy (an old numpy version)
from odbAccess import*
import numpy as np

def post_process(name_file,name_output,study_case_number):
    """
    The following functions aims to open the .odb file 
    from Abaqus to obtain:
    name_output+'Int': Electric flux on each node
                        of "nodeset" in the format
                        [label,X_coord,Y_coord,Value]
                        for each "time_step" 
    name_output+'For': Force at the closest node
                        from coord_scn for each "time_step"
    name_output+'Dis': Displacement at the closest node
                    from coord_scn for each "time_step"
    """
    # Name of the ODB file
    OdbName=name_file   
    # Open the Obd file by using an ABAQUS library
    odb=openOdb(OdbName+'.odb',readOnly=False) 
    # Name of the step 
    steptoread=odb.steps['Step-1'] 
    # Number of the time steps
    Nts=len(odb.steps['Step-1'].frames)-1 
    # Name of the nodeset to integrate the current density 
    # It can be see in abaqus 
    # Output Databases -> name of the odb
    # ->Node Sets-> Select the nodeset
    nodeset='ELECTRODE2' 
    # Node set information (coordinate,instance,label)
    nodal_set=odb.rootAssembly.nodeSets[nodeset].nodes[0]
    # Obtain the nnumber of nodes in nodeset 
    num_elec_node=len(nodal_set)       
    # Select the study case 
    if study_case_number == 1:
        # Coordinate where the load will be measured
        coord_scn1 = np.array([0.0,0.2,0.005])
        # Element information
        geometry = odb.rootAssembly.instances['PART-1-1'].elementSets['VISUALIZATION']
        # Number of elements
        num_ele = len(geometry.elements)
        # Number of nodes
        num_nodes = len(odb.rootAssembly.instances['PART-1-1'].nodeSets['SET-1'].nodes) 
        # Nodes information
        nodes_info = odb.rootAssembly.instances['PART-1-1'].nodeSets['SET-1'].nodes 
        # Idx where coord_scn1 is located (closest)
        idx_nodes = find_node_by_coord(nodes_info,coord_scn1)
    elif study_case_number == 2:
        # Coordinate where the load will be measured
        coord_scn2 = np.array([0.0,0.2,0.005])
        # Element information
        geometry = odb.rootAssembly.instances['PART-1-1'].elementSets['VISUALIZATION']
        # Number of elements
        num_ele = len(geometry.elements)
        # Number of nodes
        num_nodes = len(odb.rootAssembly.instances['PART-1-1'].nodeSets['SET-1'].nodes)  
        # Nodes information
        nodes_info = odb.rootAssembly.instances['PART-1-1'].nodeSets['SET-1'].nodes 
        # Idx where coord_scn1 is located (closest)
        idx_nodes = find_node_by_coord(nodes_info,coord_scn2)
    elif study_case_number == 3:
        # Coordinate where the load will be measured
        coord_scn3 = np.array([0.0,0.2,0.005])
        # Number of the instances where the visualization is made
        dict_key = odb.rootAssembly.instances.keys()[2]
        # Element information
        geometry = odb.rootAssembly.instances[dict_key].elementSets['VISUALIZATION']
        # Number of elements
        num_ele = len(geometry.elements)
        # Nodes information
        nodes_info = odb.rootAssembly.instances[dict_key].nodeSets['SET-1'].nodes
        # Number of nodes
        num_nodes = len(odb.rootAssembly.instances[dict_key].nodeSets['SET-1'].nodes)  
        # Idx where coord_scn1 is located (closest)
        idx_nodes = find_node_by_coord(nodes_info,coord_scn3)     
    else:
        # Print an error
        print("Error study case number not found")
    # Current data to integrate "nodeset"
    Int_info=[[] for i in range(Nts)]
    # Force data at the point idx_nodes
    For_info=[[] for i in range(Nts)]
    # Displacement data at the point idx
    Dis_info=[[] for i in range(Nts)]
    # Extraction of the connectivity  matrix
    # Not necesary we have to do this, we may just impose the equallity    
    connectivity=[[] for i in range(num_ele)]

    # Node load displacement
    node_load_dis=idx_nodes

    # Connectivity matrix extraction
    for con in range(num_ele):
        connectivity[con]=list(geometry.elements[con].connectivity)
    # Let us extract the current, force and displacement information  
    # For each time step "time_step"
    # Select from the frame 1 to the end
    for time_step in range(Nts): 
        # Frame information (description, time, etc...)
        frametoread=steptoread.frames[time_step+1] 
        # Field varible Electrical current for each element
        odbSelectResults=frametoread.fieldOutputs['SDV_JZ'] 
        # Extraction of the load information
        odbResult_force= steptoread.frames[time_step+1].fieldOutputs['RF'].values[node_load_dis].data[1]
        # Extraction of the displacement information
        odbResult_dis= steptoread.frames[time_step+1].fieldOutputs['U'].values[node_load_dis].data[1]
        # This is the element data:
        # Element label (which element correspond that node) SDV_JZ
        element=odbSelectResults.bulkDataBlocks[0].elementLabels
        # Nodal value of each node SDV_JZ
        Field_el=odbSelectResults.bulkDataBlocks[0].data 
        # Extract the values from the integration points VEP
        count=0
        VEP=np.zeros((num_ele,8))
        nodal_val_el= [[] for i in range(num_nodes)]
        nodal_val= [[] for i in range(num_nodes)]
        # For each element we collect all the values 
        # In a form VEP[element, node] 
        for ele in range(num_ele):  
            # For each node in that element       
            for i in range(8):
                # Extract the values of the integration points
                # from the Field_el and save it in a matrix with
                # dimenssion [number of element, nodes per element]
                VEP[ele,i]=Field_el[count][0]
                # add 1 the count to visits all nodes in that element
                count+=1
                # Now we save those values by introducing it, in its
                # corresponding nodal_val_el 
                # It may be introduced more than one value for each node
                # (One node may be in more than one element)
                nodal_val_el[(connectivity[ele][i])-1].append(VEP[ele,i])  
        # If there is more than one nodal value in "nodal_val_el"
        # the mean of those value are calculated
        for val in range(num_nodes):
            nodal_val[val]=sum(nodal_val_el[val])/len(nodal_val_el[val])    
        # Set the count to 0    
        count=0
        # Create a list with the len of the nodeset
        Int_info_step=[[] for i in range(num_elec_node)]
        #for each node in the nodeset selected
        for node in nodal_set:
            # We extract the [node label,X Coord, Y coord, Value] 
            Int_info_step[count]=[node.label,node.coordinates[0],node.coordinates[1],nodal_val[node.label-1]]
            count+=1
        # Save the data for the time step "time_step" 
        # in the list:
        # Int_info: that contains the electrical current info in "nodeset" section
        # For_info: contains the force information in the idx_nodes
        # Dis_info: contains the displacement information in the idx_nodes 
        Int_info[time_step]=  Int_info_step
        For_info[time_step]=  odbResult_force
        Dis_info[time_step]=  odbResult_dis
    
    data  = np.array(Int_info)
    data1 = np.array(For_info)
    data2 = np.array(Dis_info)
    np.savez(name_output+'Int', data)
    np.savez(name_output+'For', data1)
    np.savez(name_output+'Dis', data2)



def find_node_by_coord(odb_nodes_set,coord):
    """Enter the coord in a [x,y,z] format
    and you revieve the closet node of the 
    odb_nodes_set"""
    number_of_nodes=len(odb_nodes_set)
    #Iterates throughout all the nodes
    for i in range(number_of_nodes):
        #Extract the coord_i
        coord_i=np.array(odb_nodes_set[i].coordinates)
        if i==0:
            old_min = np.linalg.norm(coord_i-coord)
            old_idx = i
        else:
            new_min = np.linalg.norm(coord_i-coord)
            new_idx = i
            if new_min <= old_min:
                old_min = new_min
                old_idx = new_idx
    return old_idx

def trap_np(Y,X):
    integral=0
    for i in range(len(X)-1):
        integral = integral+Y[i]*(X[i+1]-X[i])+(Y[i+1]-Y[i])*(X[i+1]-X[i])*0.5
    return integral     


