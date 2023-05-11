# Some useful functions
smoothstep(x, x̃) = 0.5*(tanh(1e1 * (x - x̃)) + 1)
hat(x, x1, x2) = smoothstep(x, x1) - smoothstep(x, x2)
pulse(t, t1) = 2. * hat(t, t1-0.5, t1+0.5)

## Get the right symbols / indices of system variables:

"""
    sorted_symbols(symbol::Symbol, network)
    
Finds all the variable symbols containing `symbol` in `network` and sorts them based on the vertex number.
Useful if you want to apply a callback to a certain type of vertex.
"""
function sorted_symbols(symbol::Symbol, network)
    # Get the variables in the network containing `symbol`:
    symbols = syms_containing(network, symbol)
    # Get to which vertices these variables belong:
    symbol_numbers = parse.(Int, (last.(split.(String.(symbols), "_"))))
    sorted_sym_numbers = sort(collect(enumerate(symbol_numbers)); by=x->x[2])
    [symbols[el[1]] for el in sorted_sym_numbers]
end

"""
    find_neuron_number(symbol::Symbol, network, i::Integer)

Finds all the variable symbols in `network`, of vertex `i`, that contain `symbol`
"""
function find_neuron_number(symbol::Symbol, network, i::Integer)
    # Get the variables in the network containing `symbol`:
    symbols = syms_containing(network, symbol)
    # Get to which vertices these variables belong:
    symbol_numbers = parse.(Int, (last.(split.(String.(symbols), "_"))))
    @assert i ∈ symbol_numbers "The combination of variable and vertex number was not found."
    sorted_sym_numbers = sort(collect(enumerate(symbol_numbers)); by=x->x[2])
    idxs = [x[1] for x in sorted_sym_numbers[findall(x->x[2]==i, sorted_sym_numbers)]]
    symbols[idxs]
end

"""
    grouped_symbols(symbol::Symbol, network)

Returns a vector of all the symbols in the network containing 'symbol'.
If two symbols belong to the same vertex, they are returned in a vector together.
"""
function grouped_symbols(symbol::Symbol, network)
    # Get the variables in the network containing `symbol`:
    symbols = syms_containing(network, symbol)
    # Get to which vertices these variables belong:
    symbol_numbers = parse.(Int, (last.(split.(String.(symbols), "_")))) 
    indices = [symbol_numbers[1]]
    for number in symbol_numbers
        number ∉ indices ? push!(indices, number) : nothing
    end
    [(x) for x in find_neuron_number.(symbol, network, indices)]
end

"""
    filter_variables(vars, sym) 

Filters out all variables in `vars` that contain `sym`,
and verifies that this is only one variable.
"""
function filter_variables(vars, sym) 
    symbol = [s for s in vars if occursin("$sym", string(s))]
    @assert length(symbol) == 1 "Did not find one variable."
    return symbol[1]
end

"""
    function layer_symbols(vertexnumbers, network)

Given an iterable range of `vertexnumbers` and a `network`,
return all symbols pertaining to those vertices.
"""
function layer_symbols(vertexnumbers, network)
    symbollist = Vector{Symbol}[]
    for n in vertexnumbers
        vertexlist = Symbol[]
        symbols = syms_containing(network, "_$n")
        for sym in symbols
            # Check that the symbol was detected because it belong to vertex n
            # not because the number accidentally appeared elsewhere in the symbol
            if parse(Int, last(split(String(sym), "_"))) == n
                push!(vertexlist, sym)
            else
                nothing
            end
        end
        push!(symbollist, vertexlist)
    end
    symbollist
end

"""
    idx_of_symbol(nd, expr)

Find all indices of variables with symbols equal to the string, regex or symbol `expr`
"""
function idx_of_symbol(nd, expr)
    [i for (i, s) in enumerate(nd.syms) if isequal(expr, string(s))]
end
idx_of_symbol(nd, expr::Symbol) = idx_of_symbol(nd, string(expr))

number_from_symbol(sym) = parse(Int, last(split(string(sym), "_")))

## Finding spikes:

"""
    spike_finder(network, sol, V_symbols; threshold=1.)

Computes where in time the variables denoted by `V_symbols` cross `threshold` upwards.
Returns a dictionary with `variable_symbol => spike_time` of all the neurons that spike and their corresponding spike times.
"""
function spike_finder(network, sol, V_symbols; threshold=0.9)
	spikedict = Dict()
    V_indices = vcat(idx_of_symbol.(network, V_symbols)...)
	voltages = sol[V_indices, :]
	
	V_storage = zeros(size(voltages))
	V_storage[voltages .> threshold] .= 1
    V_storage = diff(V_storage, dims=2)
	
	for coord in findall(x -> x == 1, V_storage)
		if haskey(spikedict, V_symbols[coord[1]])
			push!(spikedict[V_symbols[coord[1]]], sol.t[coord[2]])
		else
			spikedict[V_symbols[coord[1]]] = [sol.t[coord[2]]]
		end
	end
	sort(spikedict; by = x -> number_from_symbol(x))
end
