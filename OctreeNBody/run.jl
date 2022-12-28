include("treeGenerator.jl");
include("treeAttributes.jl");

using CSV # importing libraries to read the data

# reading correlated data from a csv file
pos = CSV.read("nbody.csv", NamedTuple, header=true);    
particles = [[pos[1][i], pos[2][i], pos[3][i]] for i in eachindex(pos[1])]; # creating a list of particles positions

# size of the parent node
init_bound = (minimum(pos.X), maximum(pos.X), minimum(pos.Y), maximum(pos.Y), minimum(pos.Z), maximum(pos.Z));

first_div = canvas_split(init_bound)                                                               # splitting the parent node into 8 equal octants
top_leaf = TreeNode(init_bound, particles, compute_COM(particles), length(particles), first_div);  # initiating the top_leaf (parent node)
tree = build_tree(particles, top_leaf.nodes, 200);                                                 # building the tree starting with recursive depth of 200


COM_list = list_COMs(tree);  # list of centers of mass of all non-empty nodes
mass_list = list_mass(tree); # list of masses of all non-empty nodes
out = [(COM_list[i], mass_list[i]) for i in eachindex(mass_list)] # list of centers of mass and masses of all nodes

using DelimitedFiles
writedlm("out.txt", out)     # writes the output to a text file ([COM location, mass] for each node)
length(mass_list)

# uses Plotly for interactive plots in 3D
# In case any of the dependencies are missing
# uncomment following line to install PlotlyJS, DataFrames or CSV 
# import Pkg; Pkg.add(["PlotlyJS", "DataFrames", "CSV"])

using Plots            # importing library to plot the data
plotlyjs();            # for interactive plots

# function to plot nodes_boundaries
function plot_nodes(nodes_boundaries, color, alpha)
    for i in 1:8
        xMin, xMax, yMin, yMax, zMin, zMax = nodes_boundaries[i]
        plot!([xMin, xMax], [yMin, yMin], [zMin, zMin], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMax], [yMax, yMax], [zMin, zMin], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMax], [yMin, yMin], [zMax, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMax], [yMax, yMax], [zMax, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMin], [yMin, yMax], [zMin, zMin], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMin], [yMin, yMax], [zMax, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMax, xMax], [yMin, yMax], [zMin, zMin], color=color, alpha=alpha, linewidth=2)
        plot!([xMax, xMax], [yMin, yMax], [zMax, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMin], [yMin, yMin], [zMin, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMin, xMin], [yMax, yMax], [zMin, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMax, xMax], [yMin, yMin], [zMin, zMax], color=color, alpha=alpha, linewidth=2)
        plot!([xMax, xMax], [yMax, yMax], [zMin, zMax], color=color, alpha=alpha, linewidth=2)
    end
end;

# # illustrative visualization of subtended angle calculation for a single particle and single node
# p, n = 1, 1   #particle index, #node index

# node = tree.nodes[5].nodes[3].nodes[n]
# xMin, xMax, yMin, yMax, zMin, zMax = node.bounds[1][1], node.bounds[8][2], node.bounds[1][3], node.bounds[8][4], node.bounds[1][5], node.bounds[8][6];

# angle = round(angle_subtended(particles[p], node))

# scatter(pos.X, pos.Y, pos.Z, legend = false, markersize=0.6, alpha=0.8, label="All Particles");

# # plotting subtended angle
# scatter!((particles[p][1], particles[p][2], particles[p][3]), color="red", legend = false, markersize=2) # plotting the particle
# plot!([xMin, xMax], [yMin, yMax], [zMin, zMax], color="blue", linewidth=3) # plotting the node diagonal 

# # plotting lines connecting the particle to edges of the node diagonal 
# plot!([particles[p][1], xMax], [particles[p][2], yMax], [particles[p][3], zMax], color="blue", linewidth=3) # plotting the subtended angle
# plot!([particles[p][1], xMin], [particles[p][2], yMin], [particles[p][3], zMin], color="blue", linewidth=3) # plotting the subtended angle


# # plotting the node boundaries
# plot_nodes(tree.bounds, "red", 0.4);
# plot_nodes(tree.nodes[5].bounds, "red", 0.2);
# plot_nodes(tree.nodes[5].nodes[3].bounds, "black", 0.3);
# display(plot!(dpi=200, size=(900, 800), camera=(15, 25, 250), background_color="white", title="Canvas Spliting (node depth = 3) and Subtended Angle ~ $angle",
#                 xlabel="x", ylabel="y", zlabel="z"))
# png("canvas_split")


# # plotting dependance on tree parameter θ
# treeParam_range = 1:5:65
# scatter()
# for p in 1:24:length(particles)
#     no_of_nodes = [length(list_nodes(p, theta, tree)) for theta in treeParam_range]
#     scatter!(treeParam_range, no_of_nodes, mode="markers+lines", label="Particle $p", alpha=0.5)
# end
# display(scatter!(xlabel="Tree Parameter θ", ylabel="Number of interactions to consider"))
# png("no_of_nodes")

# scatter(pos.X, pos.Y, pos.Z, legend = false, markersize=0.6, alpha=0.7, label="All Particles");