using Markowitz
using LinearAlgebra
import Plots

assets = [ "Bonds - US Government"
           "Bonds - US Corporate"
           "Bonds - International"
           "Bonds - High Yield"
           "Bonds - Bank Loans"
           "Bonds - Emerging USD Debt"
           "Bonds - Emerging Local Debt"
           "Alternative - Emerging FX"
           "Alternative - Commodities"
           "Alternative - REITs"
           "Stocks - US Large"
           "Stocks - US Small"
           "Stocks - International"
           "Stocks - Emerging" ]

V = [ 133.0  39.0  36.0  -17.0  -28.0   31.0   -5.0   -6.0  -34.0  -10.0  -46.0  -68.0  -52.0  -63.0
       39.0  16.0  17.0   10.0    0.0   22.0   13.0    7.0    1.0   19.0    3.0    0.0    5.0    9.0
       36.0  17.0  38.0   19.0    3.0   34.0   45.0   29.0   37.0   48.0   20.0   21.0   42.0   50.0
      -17.0  10.0  19.0   98.0   69.0   65.0   73.0   43.0   68.0  151.0   99.0  130.0  121.0  165.0
      -28.0   0.0   3.0   69.0   84.0   34.0   42.0   25.0   55.0   99.0   64.0   93.0   87.0  119.0
       31.0  22.0  34.0   65.0   34.0   93.0   83.0   48.0   57.0  115.0   76.0   85.0   99.0  145.0
       -5.0  13.0  45.0   73.0   42.0   83.0  154.0   87.0  111.0  168.0  113.0  143.0  165.0  225.0
       -6.0   7.0  29.0   43.0   25.0   48.0   87.0   57.0   75.0   95.0   71.0   89.0  104.0  142.0
      -34.0   1.0  37.0   68.0   55.0   57.0  111.0   75.0  286.0  117.0  101.0  125.0  164.0  239.0
      -10.0  19.0  48.0  151.0   99.0  115.0  168.0   95.0  117.0  491.0  231.0  327.0  259.0  322.0
      -46.0   3.0  20.0   99.0   64.0   76.0  113.0   71.0  101.0  231.0  214.0  256.0  217.0  269.0
      -68.0   0.0  21.0  130.0   93.0   85.0  143.0   89.0  125.0  327.0  256.0  380.0  265.0  342.0
      -52.0   5.0  42.0  121.0   87.0   99.0  165.0  104.0  164.0  259.0  217.0  265.0  297.0  359.0
      -63.0   9.0  50.0  165.0  119.0  145.0  225.0  142.0  239.0  322.0  269.0  342.0  359.0  556.0 ]

E = [ 0.1 0.7 0.8 2.3 2.2 1.9 5.6 5.6 2.2 1.3 0.7 -0.1 4.1 7.2 ]

class = vec([ :FI :FI :FI :FI :FI :FI :FI :ALT :ALT :ALT :EQ :EQ :EQ :EQ ])
colors = map(c-> c==:FI ? "blue" : (c==:ALT ? "green" : "red"),class)
alphas = 0.5 .+ 0.5*[ (1:7)/7 ; (1:3)/3 ; (1:4)/4 ]
    
function plot_frontier()
    m = f.problem
    unused = vec(sum(f.weights,1).==0)
    xrange = [0, max(maximum(sqrt.(diag(m.V))),maximum(f.vol))]*1.1
    yrange = [minimum(m.E),max(maximum(m.E),maximum(f.ret))]+(maximum(m.E)-minimum(m.E))*[-1,1]*0.1
    dense = smooth(f)
    p1 = Plots.plot(dense[:,1], dense[:,2], legend=false, label="", color="orange", width=2)
    for i = 1:length(m.E)
        Plots.scatter!([sqrt.(diag(m.V)[i])], [m.E[i]], label=m.names[i],
                       markersize=4-1.5*unused[i], markercolor=colors[i], markeralpha=alphas[i],
                       markerstrokecolor=colors[i], markerstrokealpha=alphas[i])
    end
    Plots.scatter!(f.vol, f.ret, legend=false, label="", color="orange", markerstrokecolor="orange", markersize=1)
    Plots.xlabel!("Volatility")
    Plots.ylabel!("Return")
    Plots.xlims!(Tuple(xrange))
    Plots.ylims!(Tuple(yrange))
    weights = f.weights[:,.!unused]
    indices = [1]
    for i in 2:size(weights,1)
        w0 = weights[i-1,:]
        w1 = weights[i,:]
        if sum((sign.(w0).*sign.(w1)).==-1)>0
            for j in 1:length(w0)
                if sign.(w0[j]) != sign.(w1[j])
                    indices = [indices ; i-1+abs(w0[j]/(w1[j]-w0[j]))]
                end
            end
        end
        indices = [indices ; i]
    end
    extweights = zeros(0,size(weights,2))
    extret = []
    for idx in sort(indices)
        if rem(idx,1) == 0
            extweights = [ extweights ; weights[[Int(idx)],:] ]
            extret = [ extret ; f.ret[Int(idx)] ]
        else
            w0 = weights[[Int(floor(idx))],:]
            w1 = weights[[Int(ceil(idx))],:]
            r0 = f.ret[[Int(floor(idx))]]
            r1 = f.ret[[Int(ceil(idx))]]
            extweights = [ extweights ; w0+rem(idx,1)*(w1-w0) ]
            extret = [ extret ; r0+rem(idx,1)*(r1-r0) ]
        end
    end
    wrange = [ minimum(sum((extweights.<0) .* (extweights),2)), maximum(sum((extweights.>0) .* (extweights),2)) ]
    wrange = [ wrange[1], wrange[1] + 2.5*(wrange[2]-wrange[1])]
    ylong = [ zeros(size(extweights,1),1) cumsum((extweights.>0) .* (extweights),2) ]
    yshort = [ zeros(size(extweights,1),1) cumsum((extweights.<0) .* (extweights),2) ]
    p2 = Plots.plot(legend=:right, background_color_legend="#ffffff00", legendfont=Plots.font(5), grid=nothing)
    Plots.xticks!([0,0.5,1])
    for i = 1:(size(ylong,2)-1)
        Plots.plot!(Plots.Shape([ vec(ylong[:,i]); reverse(vec(ylong[:,i+1])) ], [ extret; reverse(extret)]),
                    color=colors[.!unused][i], alpha=alphas[.!unused][i], label=m.names[.!unused][i])
    end
    for i = 1:(size(yshort,2)-1)
        Plots.plot!(Plots.Shape([ vec(yshort[:,i]); reverse(vec(yshort[:,i+1])) ], [ extret; reverse(extret)]),
                    color=colors[.!unused][i],alpha=alphas[.!unused][i], label="")
    end
    Plots.xlims!(Tuple(wrange))
    Plots.ylims!(Tuple(yrange))
    Plots.plot(p1, p2, size=(800, 400))
end


m = markowitz(E, V, names=assets)
unit_sum(m) # total weight = 100%
f=frontier(m)
plot_frontier()
optimal(f) # volatility, return and weights for the minimum variance portofolio
optimal(f,4) # volatility, return and weights for the optimal portofolio with return = 4


m=markowitz(E, V, names=assets, lower=0.05, upper=0.15) # min 5%, max 15% per position
unit_sum(m)
f=frontier(m)
plot_frontier()


m=markowitz(E, V, names=assets, # asset bounds by class: stocks -10/30, bonds 0/20, alt. 0/10
            lower = -0.1 * (class .== :EQ),
            upper = 0.3 * (class .== :EQ) + 0.2 * (class .== :FI) + 0.1 * (class .== :ALT))
unit_sum(m)
add_constraint(m, 1 * (class .== :EQ), '>', 0.3) # net equity exposure between 30% and 60%
add_constraint(m, 1 * (class .== :EQ), '<', 0.6)
add_constraint(m, [1 1 0 0 0 0 0 0 0 0 0 0 0 0], '=', 0.25) # US govt + Investment Grade = 25%
f=frontier(m)
plot_frontier()
