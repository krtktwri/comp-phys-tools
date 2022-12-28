###############################
### treeGenerator.jl Module ###
###############################
"For turning n-body position data into octrees"

# creating custom julia structure to store tree type data 
struct TreeNode
    bounds # stores the spatial bounds of the node
    bodies # stores the particles in the node
    COM    # stores center of mass of the node
    mass   # stores the mass of the node (in solar masses)
    nodes  # stores the child nodes
end;

# function to divide the node into 8 equal octants
function canvas_split(bound)

    # unpacking bounds of the full region
    xMin, xMax, yMin, yMax, zMin, zMax = bound

    # Calculating midpoints of the total region
    xMid = (xMax + xMin)/2.0
    yMid = (yMax + yMin)/2.0
    zMid = (zMax + zMin)/2.0

    # Defining the new bounds of the 8 octants
    new_bounds = [(xMin, xMid, yMin, yMid, zMin, zMid), # label (1, 1, 1)
                  (xMid, xMax, yMin, yMid, zMin, zMid), # label (2, 1, 1)
                  (xMin, xMid, yMid, yMax, zMin, zMid), # label (1, 2, 1)
                  (xMin, xMid, yMin, yMid, zMid, zMax), # label (1, 1, 2)
                  (xMid, xMax, yMid, yMax, zMin, zMid), # label (2, 2, 1)
                  (xMid, xMax, yMin, yMid, zMid, zMax), # label (2, 1, 2)
                  (xMin, xMid, yMid, yMax, zMid, zMax), # label (1, 2, 2)
                  (xMid, xMax, yMid, yMax, zMid, zMax)  # label (2, 2, 2)
    ]   

    return new_bounds # output is a list of 8 tuples
end;


# function to locate which octant a particle belongs to
function which_node(particle_location, nodes_boundaries)

    xMid = nodes_boundaries[1][2]   # extracting node boundary information
    yMid = nodes_boundaries[1][4]   # extracting node boundary information
    zMid = nodes_boundaries[1][6]   # extracting node boundary information

    x, y, z = particle_location  # unpacking particle locations
 
    xPass = x < xMid    # checking if x-coordinate passes the y-z plane
    yPass = y < yMid    # checking if y-coordinate passes the x-z plane
    zPass = z < zMid    # checking if z-coordinate passes the x-y plane

    identifier = (xPass, yPass, zPass) # creating an identifier to check which octant the particle is in

    if identifier == (true, true, true)
        n = 1
    elseif identifier == (false, true, true)
        n = 2
    elseif identifier == (true, false, true)
        n = 3
    elseif identifier == (true, true, false)
        n = 4
    elseif identifier == (false, false, true)
        n = 5
    elseif identifier == (false, true, false)
        n = 6
    elseif identifier == (true, false, false)
        n = 7
    else
        n = 8
    end

    return n
end;

# The following function builds the tree recursively 
# Input parameters: parent node, particles in the parent node, boundaries of the parent node, max allowed recursion depth
# Output: a tree structure with the parent node as the root

function build_tree(particles, nodes_boundaries, depth)

    # Base case: if there is only one particle in the node or the depth of the tree is 0
    if length(particles) <= 1 || depth == 0
        return TreeNode(nodes_boundaries, particles, compute_COM(particles), length(particles), nodes_boundaries)
    end

    # keeping track of which particle belongs to which octant
    child = [which_node(particles[i], nodes_boundaries) for i in eachindex(particles)]

    # Recursively building the tree
    new_nodes = [build_tree(particles[child .== i], canvas_split(nodes_boundaries[i]), depth-1) for i in 1:8]

    return TreeNode(nodes_boundaries, particles, compute_COM(particles), length(particles), new_nodes)
end;