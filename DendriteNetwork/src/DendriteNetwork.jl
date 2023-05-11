module DendriteNetwork

using OrdinaryDiffEq, Graphs, Flux, NetworkDynamics, Plots 

include("utils.jl")
include("dynamics.jl")
include("graphs.jl")
include("plotting.jl")
include("auxiliary_networks.jl")

# Can be used to manipulate dynamics:
export smoothstep
export hat
export pulse
# Useful for dealing with the functions of NetworkDynamics.jl:
export sorted_symbols
export find_neuron_number
export grouped_symbols
export filter_variables
export layer_symbols
export idx_of_symbol
export number_from_symbol
# Finding spikes:
export spike_finder

# Dynamics:
export soma!
export input_soma!
export soma_output!
export synapse!
export electric_coupling!
# Callbacks:
export spiking_cb
export plateau_cb
export test_spiking_cb, test_plateau_cb

# Creating auxiliary network:
export Θ
export build_binary_network
export build_binary_network2

# Building graphs:
export create_network_graph

# Help with plotting:
export plotsym

end
