# To create the graphs:

"""
    create_network_graph(weights::Array)

Create a weighted, directed graph from an array of weights.
The weights must be matrices specifying the connection strenghts and size of the
layers.

We shall use the convention here that ``(W_{i})^{m \times n}``, specifying the weights
from layer ``i`` to layer ``i+1``, means layer ``i`` has ``n`` neurons and layer ``i+1``
has ``m`` neurons, and ``(W_i)_{kn}`` is the weight element from ``n`` to ``k``. 
(In a linear layer this implies you left multiply the weights
with the input, i.e. ``y = Wx``.)

Returns `[g, edge_weights, vertex_numbers]`,
where `g` is the graph object, `edge_weights` the list of weights created
from the `weights` argument, vertex_numbers the integers of the nodes for each layer.
"""
function create_network_graph(weights::Array)
    sizes = vcat([reverse(size(weights[1]))...], [size(W, 1) for W in weights[2:end]])

    startidxs = [sum(sizes[1:(i-1)]) + 1 for i in 1:length(sizes)]
    endidxs = [sum(sizes[1:i]) for i in 1:length(sizes)]
    vertexnums = [b:e for (b, e) in zip(startidxs, endidxs)]
    
    edge_weights = []
    g = SimpleGraph(sum(sizes))

    # Iterates over all (source, destination) vertex pairs,
    # as well as the corresponding indices of these pairs for the weights matrix.
    # `begin_vertices` are the vertex numbers in the source layers
    # `end_vertices` are the vertex numbers in the destination layer
    # `begin_offset` is used to turn the vertex numbers into indices (the first vertex becomes 1)
    # I.e. `i = v - begin_offset`, with `i` the index and `v` the vertex number.
    # Similarly, `end_offset` is used to turn the destination vertices to indices.
    for (i, (begin_vertices, end_vertices, begin_offset, end_offset)) in enumerate(zip(vertexnums[1:(end-1)], vertexnums[2:end], (startidxs .- 1)[1:(end-1)], endidxs[1:(end-1)]))
        for (source, source_index) in zip(begin_vertices, begin_vertices .- begin_offset)
            for (destination, destination_index) in zip(end_vertices, end_vertices .- end_offset)
                weight = weights[i][destination_index, source_index]
                if weight != 0.
                    add_edge!(g, source, destination) 
                    push!(edge_weights, weight)
                else
                    nothing
                end
            end
        end
    end

    return g, edge_weights, vertexnums
end
