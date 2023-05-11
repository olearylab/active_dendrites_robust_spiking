# Plotting a variable from symbol:
"""
    plotsym(sol, network, symbol)

Plot the graph of `symbol` in `sol`.
"""
function plotsym(sol, network, symbol)
    @assert length(idx_of_symbol(network, symbol)) == 1 "That symbol was not found in the network."
    ts = sol.t
    vals = transpose(sol[idx_of_symbol(network, symbol), :])
    plot(ts, vals, c=:black)
end
