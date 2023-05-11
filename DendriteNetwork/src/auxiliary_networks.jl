layernum_to_param_idxs(i::Int) = (2i-1, 2i)

## For use with Flux.jl:

"""
    Θ(x)

Heaviside stepfunction: 0 for x < 0, 0.5 for x == 0, 1 for x > 0.
"""
Θ(x) = 0.5 * (sign(x) + 1) # Heaviside step function

"""
    build_binary_network(params)

Construct a binary ANN from `params`.
Requires that params is a list of `[W, b, W, b, ...]`.
The activations of all layers are a heaviside function,
    except for the last layer which has a sigmoid.
"""
function build_binary_network(params)
	k = length(params)

	activ_f1 = Θ
	activ_f2 = sigmoid
	
	layers = []
	for i in 1:2:k
		dims = reverse(size(params[i]))
		i < (k-1) ? push!(layers, Dense(dims..., activ_f1)) :
		push!(layers, Dense(dims..., activ_f2))
	end

	m = Chain(layers...)
	
	ps = let
		ps = []
		# Ensure that weights are clamped between [0, ∞]
		for (i, W, b) in zip(1:2:k, params[1:2:k], params[2:2:k])
			if i < (k-1)
				push!(ps, clamp.(W, 0f0, Inf32))
				push!(ps, b)
			else
				push!(ps, W)
				push!(ps, b)
			end
		end
		Flux.params(ps)
	end

	Flux.loadparams!(m, ps)
	m
end


function build_binary_network2(params)
	k = length(params)

	layers = []
	for i in 1:2:k
		dims = reverse(size(params[i]))
		push!(layers, Dense(dims..., Θ))
	end

    ps = let
        ps = Array{Float32}[]
        for (i, W, b) in zip(1:(k÷2), params[1:2:k], params[2:2:k])
            if iseven(i)
                push!(ps, clamp.(W, 0f0, Inf32))
                push!(ps, b)
            elseif isodd(i)
                push!(ps, W)
                push!(ps, b)
            end
        end
        Flux.params(ps)
    end

	m = Chain(layers...)
	Flux.loadparams!(m, ps)

	m
end
