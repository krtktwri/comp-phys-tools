# angles subtended on a particle by a node
function angle_subtended(particle_location, node)
    
    # unpacking particle location
    x, y, z = particle_location
    
    # unpacking node boundaries
    xMin, xMax, yMin, yMax, zMin, zMax = node.bounds[1][1], node.bounds[8][2], node.bounds[1][3], node.bounds[8][4], node.bounds[1][5], node.bounds[8][6]

    d = sqrt((xMax - xMin)^2 + (yMax - yMin)^2 + (zMax - zMin)^2)        # distance between the two corners of the node
    xMid, yMid, zMid = (xMin + xMax)/2, (yMin + yMax)/2, (zMin + zMax)/2 #  center of the node

    r = sqrt((x - xMid)^2 + (y - yMid)^2 + (z - zMid)^2) # distance between the center of the node and the particle

    # calculating the angle subtended by the node on the particle
    angle = 2*atan(d/r)

    return angle*180/pi # output is in degrees
end;   


##############################
######## COM LIST ############
##############################

# Center of Mass computing function for equal mass particles
function compute_COM(particle_list)

    N = length(particle_list)      # total number of particles 
    vecAdd = [0,0,0]               # initiating array to compute center of mass
    i = 1                          # initiating iterating index

    while i <= N
        vecAdd += particle_list[i] # adding each particle's contribution to the total
        i += 1 
    end
 
    return vecAdd/N                 # dividing my total number of particles (note it only works with equal mass particles)
end;

# recursive function to get the centers of mass of all nodes
function get_COMs(tree, COM_list)
    # if node has only one particle storing particle position as COM of the node
    if length(tree.bodies) == 1 
        push!(COM_list, tree.bodies[1])
    end

    push!(COM_list, tree.COM)                 # appending the center of mass position to the list
    for i in 1:8                              # iterating over the 8 octants
        if typeof(tree.nodes[i]) == TreeNode  # convoluted way of checking if the node is empty
            get_COMs(tree.nodes[i], COM_list) # recursive step
        end
    end
end;

# Recursively tracing the tree to find the centers of mass of the non-empty nodes
# Input: tree structure
# Output: list of center of mass locations for all "non-empty" nodes

function list_COMs(tree)

    COM_list = []   # intiating an empty list to store the centers of mass

    get_COMs(tree, COM_list)                          # updating list
    COM_list = [x for x in COM_list if !isnan(x[1])]  # removing the NaNs from the list

    return COM_list
end;

################################
######### MASS LIST ############
################################

# recursive function to get the masses of all nodes
function get_mass(tree, mass_list)

    # if node has only one particle, storing particle mass as mass of the node
    if length(tree.bodies) == 1 
        push!(mass_list, tree.mass)
    end

    push!(mass_list, tree.mass)                # appending the center of mass position to the list
    for i in 1:8                               # iterating over the 8 octants
        if typeof(tree.nodes[i]) == TreeNode   # convoluted way of checking if the node is empty
            get_mass(tree.nodes[i], mass_list) # recursive step
        end
    end
end;

# Recursively tracing the tree to find the centers of mass of the non-empty nodes
# Input: tree structure
# Output: list of masses of corresponding to the COMs generated by list_COM() function

function list_mass(tree)

    mass_list = []   # intiating an empty list to store the centers of mass

    get_mass(tree, mass_list)                          # updating list
    mass_list = [x for x in mass_list if x > 0]        # removing the NaNs from the list

    return mass_list
end;


################################
######### NODE LIST ############
################################

# recusrively finding angle subtended by each non-empty node on a particle
function get_nodes(particle, tree, theta_cutoff, outLists)

    # base case: angle_subtended is less than theta_cutoff
    if angle_subtended(particle, tree) <= theta_cutoff || length(tree.bodies) == 1
        return push!(outLists, tree.COM)
    end

    # recursively finding the centers of mass of the nodes to be considered
    for i in 1:8     # iterating through octants
        if typeof(tree.nodes[i]) == TreeNode
            push!(outLists, get_nodes(particle, tree.nodes[i], theta_cutoff, outLists))
        end
    end

    return outLists
end;

function list_nodes(p, cutoff, tree) 

    consider_COM = []    # initiating a list of COM positions to consider for the given particle 
    consider_COM = get_nodes(particles[p], tree, cutoff, consider_COM) # populating the list

    # removing non-vector elements from the list
    consider_COM = [x for x in consider_COM if typeof(x) == Vector{Float64} && !isnan(x[1])];

    return consider_COM
end;