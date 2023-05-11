
function soma_function!(dv, v, edges, P, t)
    # dv[1] = dV, dv[2] = ds, dv[3] = dg_r, dv[4] = dT
    V, s, g_r, T = v
    τ_m, τ_s = P

    I_in = 0.
    for e in edges
        I_in += (e[1] / τ_m)
    end


    dv[1] = if T < 0.
        (-V + I_in - g_r * V) / τ_m # dV
    elseif T ≥ 0.
        0. # No integration during refractory period
    end
    dv[2] = -s / τ_s # ds
    dv[3] = -g_r / (10*τ_m) # dg_r
    dv[4] = -smoothstep(T, 0.) # dT

end
soma! = ODEVertex(; f=soma_function!, dim=4, sym=[:V_soma, :s_soma, :g_r_soma, :T_soma])

function input_soma!(input_current)
    function input_soma_function!(dv, v, edges, P, t; input=input_current)
        # dv[1] = dV, dv[2] = ds, dv[3] = dg_r
        V, s, g_r, T = v
        τ_m, τ_s = P
        dv[1] = if T < 0.
            (-V + - g_r * V + input(t)) / τ_m # dV
        elseif T ≥ 0.
            0. # No integration during refractory period
        end
        dv[2] = -s / τ_s # ds
        dv[3] = -g_r / (10*τ_m) # dg_r
        dv[4] = -smoothstep(T, 0.) # dT

        # dV += I_in
        for e in edges
            dv[1] += (e[1] / τ_m)
        end
    end
    input_soma! = ODEVertex(; f=input_soma_function!, dim=4, sym=[:V_soma, :s_soma, :g_r_soma, :T_soma])
end

function soma_output_function!(e, v_source, v_dest, K, t)
    e[1] = K[1] * v_source[2] #* (K[2] - v_dest[1]) # w * s_source(t) * (E - V_dest(t))
    e[2] = 0. # No flow in other direction
end
soma_output! = StaticEdge(; f=soma_output_function!, dim=2, coupling=:fiducial, sym=[:s_syn, :empty])

function synapse_function!(dv, v, edges, P, t)
    # dv[1] = dV, dv[3] = dT
    V, T = v
    τ_s = only(P)

    I = 0.
    for e in edges
        I += (e[1] / τ_s) # NB: already multiplied by R in the coupling function.
    end

    dv[1] = if T < 0. #ifelse(T < 0, (-V + I) / τ_s, 1e-3) # dV
        (-V + I) / τ_s # Only integrate charge when the timer is nonactive (T < 0.)
    elseif T ≥ 0.
        0. # No integration when timer is on
    end
    dv[2] = -smoothstep(T, 0.) # dT
end
synapse! = ODEVertex(; f=synapse_function!, dim=2, sym=[:V_syn, :T_syn])

function electric_coupling_function!(e, v_source, v_dest, P, t)
    R = only(P)
    e[1] = R * (v_source[1] - v_dest[1]) # K * (V_soma - V_syn)
    e[2] = R * (v_dest[1] - v_source[1]) # Equal, opposite flow in other direction
end
electric_coupling! = StaticEdge(; f=electric_coupling_function!, dim=2, coupling=:fiducial, sym=[:I_electric, :I_electric])

# Callbacks for spiking events, plateaus:

function plateau_cb(network, name::Symbol; θ = 1.0, T = 10.)

    # Define a margin between threshold and region where callback gets triggered
    # To prevent continuous triggering due to floating point roundoffs etc.
    ϵ = 1e-5 

    synapse_symbols = grouped_symbols(name, network)

    timer_activation = let
        callbacks = []
        for (V_sym, T_sym) in synapse_symbols
            V_idx = only(idx_of_symbol(network, V_sym))
            T_idx = only(idx_of_symbol(network, T_sym))

            function condition(u, t, integrator)
                u[V_idx] > θ + ϵ
            end
            function affect!(integrator)
                integrator.u[V_idx] = θ
                integrator.u[T_idx] = T
            end
            push!(callbacks, DiscreteCallback(condition, affect!))
        end
        CallbackSet(callbacks...)
    end
    plateau_holding = let
        callbacks = []
        for (V_sym, T_sym) in synapse_symbols
            V_idx = only(idx_of_symbol(network, V_sym))
            T_idx = only(idx_of_symbol(network, T_sym))

            function condition(u, t, integrator) 
                u[T_idx] > 0.
            end
            function affect!(integrator)
                integrator.u[V_idx] = θ
            end
            push!(callbacks, DiscreteCallback(condition, affect!))
        end
        CallbackSet(callbacks...)
    end
    # NB: The order in which these are put matters!
    # Discrete callbacks are called in the order they are provided
    # We need the plateau callback to be called before the timer
    return CallbackSet(plateau_holding, timer_activation)
end

function spiking_cb(network, name; θ = 0.99, T_refractory = 5.)
	callbacks = []
    for soma_variables in grouped_symbols(name, network)
        V_symbol = filter_variables(soma_variables, :V)
        g_r_symbol = filter_variables(soma_variables, :g_r)
        s_symbol = filter_variables(soma_variables, :s_)
        T_symbol = filter_variables(soma_variables, :T)

        V_index = only(idx_of_symbol(network, V_symbol))
        g_r_index = only(idx_of_symbol(network, g_r_symbol))
        s_index = only(idx_of_symbol(network, s_symbol))
        T_index = only(idx_of_symbol(network, T_symbol))
    
        condition(u, t, integrator) = u[V_index] - θ
        function affect!(integrator) 
            integrator.u[V_index] = 0.
            integrator.u[g_r_index] = 100.
            integrator.u[s_index] = 1.
            integrator.u[T_index] = T_refractory
        end
        push!(callbacks, ContinuousCallback(condition, affect!))
    end
    CallbackSet(callbacks...)
end

function test_spiking_cb(network, name; θ = 0.99, T_refractory = 5.)

    function spike_condition(out, u, t, integrator)
        for (i, soma_variables) in enumerate(grouped_symbols(name, network))
            V_symbol = filter_variables(soma_variables, :V)
            V_index = only(idx_of_symbol(network, V_symbol))
            
            out[i] = u[V_index] - θ
        end
    end

    function spike_affect!(integrator, idx)
        soma_variables = grouped_symbols(name, network)[idx]

        V_symbol = filter_variables(soma_variables, :V)
        g_r_symbol = filter_variables(soma_variables, :g_r)
        s_symbol = filter_variables(soma_variables, :s_)
        T_symbol = filter_variables(soma_variables, :T)

        V_index = only(idx_of_symbol(network, V_symbol))
        g_r_index = only(idx_of_symbol(network, g_r_symbol))
        s_index = only(idx_of_symbol(network, s_symbol))
        T_index = only(idx_of_symbol(network, T_symbol))

        integrator.u[V_index] = 0.
        integrator.u[g_r_index] = 100.
        integrator.u[s_index] = 1.
        integrator.u[T_index] = T_refractory
    end

    cb_length = length(grouped_symbols(name, network))

    return VectorContinuousCallback(spike_condition, spike_affect!, cb_length)
end

function test_plateau_cb(network, name; θ = 1.0, T = 10.)
    ϵ = 1e-5

    synapse_symbols = grouped_symbols(name, network)
    V_symbols = filter_variables.(synapse_symbols, :V)
    T_symbols = filter_variables.(synapse_symbols, :T)

    V_indices = vcat(idx_of_symbol.(network, V_symbols)...)
    T_indices = vcat(idx_of_symbol.(network, T_symbols)...)

    function plateau_timer_condition(out, u, t, integrator)
        for (i, synapse_symbols) in enumerate(synapse_symbols)
            V_symbol = filter_variables(synapse_symbols, :V)
            V_index = only(idx_of_symbol(network, V_symbol))
            
            out[i] = u[V_index] > θ + ϵ
        end
    end
    function plateau_timer_condition2(u, t, integrator)
        any(u[V_indices] .> θ + ϵ)
    end

    function plateau_timer_affect!(integrator, idx)
        synapse_variables = synapse_symbols[idx]

        V_sym = filter_variables(synapse_variables, :V)
        T_sym = filter_variables(synapse_variables, :T)

        V_idx = only(idx_of_symbol(network, V_sym))
        T_idx = only(idx_of_symbol(network, T_sym))

        integrator.u[V_idx] = θ
        integrator.u[T_idx] = T
    end
    function plateau_timer_affect2!(integrator)
        plateau_V_meta_idxs = findall(integrator.u[V_indices] .> (θ + ϵ))
        plateau_T_meta_idxs = T_indices[plateau_V_meta_idxs]

        integrator.u[plateau_T_meta_idxs] .= T
    end


    function plateau_hold_condition(out, u, t, integrator)
        for (i, synapse_symbols) in enumerate(grouped_symbols(name, network))
            T_sym = filter_variables(synapse_symbols, :T)
            T_idx = only(idx_of_symbol(network, T_sym))
            
            out[i] = u[T_idx] > 0.
        end
    end
    function plateau_hold_condition2(u, t, integrator)
        any(u[T_indices] .> 0.)
    end


    function plateau_hold_affect!(integrator)
        plateau_T_meta_idxs = findall(integrator.u[T_indices] .> 0.)
        plateau_V_meta_idxs = V_indices[plateau_T_meta_idxs]

        integrator.u[plateau_V_meta_idxs] .= θ
    end

    cb_length = length(synapse_symbols)

    timer_cb = DiscreteCallback(plateau_timer_condition2, plateau_timer_affect2!)
    hold_cb = DiscreteCallback(plateau_hold_condition2, plateau_hold_affect!)

    return CallbackSet(timer_cb, hold_cb)
end
