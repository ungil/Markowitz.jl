module Markowitz
using LinearAlgebra

export markowitz, unit_sum, add_constraint, frontier, smooth, optimal

function simplex(mu,A,b,L,U,epsilon=1e-10)
    n = size(mu,1)
    m = size(A,1)
    Ai = zeros(m,m)
    L = [ L; zeros(m) ]
    U = [ U; Inf*ones(m) ]
    z = zeros(n+m)
    x = zeros(n+m)
    state = Vector{Symbol}(n+m)
    nartificial = m
    adjRate = Vector{Float64}(n+m)
    phase = 1
    OUT = 1:n
    state[OUT] = :LOW
    x[OUT] = L[1:n]
    z[OUT] = 0
    IN = n + (1:m)
    tmp = b - A * L[1:n]
    for i in 1:m
        Ai[i,i] = (tmp[i]<0 ? -1 : 1)
    end
    A = hcat(A, Ai)
    state[IN] = :IN
    x[IN] = abs.(tmp)
    z[IN] = -1
    loop = 0
    while loop < 10*(n+m)
        loop = loop + 1
        price = -z[IN]' * Ai[1:length(IN),:]
        profit = ( z[OUT]' + price * A[:,OUT]) .* map(x-> x==:HIGH ? -1 : 1, state[OUT])'
        if maximum(profit) < epsilon
            if phase == 1
                error("unfeasible or degenerate")
            else
                muChanged = false
                for i in 1:length(OUT)
                    if profit[i] > -epsilon
                        muChanged = true
                        warn("tweaking mu[$i] to ensure the solution is unique")
                        mu[OUT[i]] = mu[OUT[i]] + epsilon * (state[OUT[i]] == :LOW ? -1 : 1)
                    end
                end
                return(Dict(:IN => IN, :OUT => OUT, :state => state[1:n],
                            :x => x[1:n], :Ai => Ai, :mu => muChanged ? mu : nothing))
            end
        end
        enters = OUT[indmax(profit)]
        directionIn = (state[enters]==:HIGH ? (:DOWN) : (:UP))
        adjRate = -Ai * A[:,enters] * (directionIn==:DOWN ? -1 : 1)
        exits = enters
        exitsIdx = NaN
        directionOut = directionIn
        if ( enters > n || U[enters] == Inf )
            theta = Inf
        else
            theta = U[enters]-L[enters]
        end
        for i in 1:m
            if adjRate[i] < -epsilon
                tmp = ( L[IN[i]] - x[IN[i]] ) /adjRate[i]
                if tmp < theta
                    theta = tmp
                    exits = IN[i]
                    exitsIdx = i
                    directionOut = :DOWN
                end
            elseif (adjRate[i] > epsilon) && (U[IN[i]] < Inf)
                tmp = ( U[IN[i]] - x[IN[i]] ) / adjRate[i]
                if tmp < theta
                    theta = tmp
                    exits = IN[i]
                    exitsIdx = i
                    directionOut = :UP
                end
            end
        end
        if theta == Inf
            error("simplex unbounded")
        end
        x[IN] = x[IN] + theta * adjRate
        x[enters] = x[enters] + theta * (directionIn==:UP ? 1 : -1)
        OUT = OUT[OUT.!=enters]
        IN = sort([ enters; IN ])
        state[enters] = :IN
        IN = IN[IN.!=exits]
        if exits <= n
            OUT = sort([ exits; OUT ])
        end
        state[exits] = (directionOut==:UP ? (:HIGH) : (:LOW))
        if !isnan(exitsIdx)
            for i in 1:m
                if i != exitsIdx
                    Ai[i,:] = Ai[i,:] - adjRate[i] / adjRate[exitsIdx] * Ai[exitsIdx,:]
                end
            end
            Ai[exitsIdx,:] = Ai[exitsIdx,:] / adjRate[exitsIdx] * (directionIn==:UP ? -1 : 1)
            entersIdx = find(IN .== enters)[1]
            if entersIdx > exitsIdx
                tmp = Ai[exitsIdx,:]
                for i in exitsIdx:(entersIdx-1)
                    Ai[i,:] = Ai[i+1,:]
                end
                Ai[entersIdx,:]= tmp
            elseif entersIdx < exitsIdx
                tmp = Ai[exitsIdx,:]
                for i in exitsIdx:-1:(entersIdx+1)
                    Ai[i,:] = Ai[i-1,:]
                end
                Ai[entersIdx,:] = tmp
            end
        end
        if phase == 1 && exits > n
            nartificial = nartificial-1
            if nartificial == 0
                phase = 2
                z = mu
            end
        end
    end
end

struct Change
    inout::Symbol
    direction::Symbol
    variable::Int
    lambda::Float64
end

function cla(C,mu,A,b,L,U,epsilon=1e-10,lambdaTarget=1e-6)
    @assert det(C)>=0
    scale = maximum(C)
    C = C/scale
    n = size(mu,1)
    m = size(A,1)
    alpha = zeros(n+m)
    beta = zeros(n+m)
    xi = zeros(n+m)
    bbar = zeros(n+m)
    extAi = zeros(n,m)
    M = [C     A'    ;
         A zeros(m,m)]
    lambdaE = Inf
    lambdaEold = NaN
    weights = Array{Float64}(0,n)
    ret = Vector{Float64}(0)
    vol = Vector{Float64}(0)
    lambda = Vector{Float64}(0)
    simp = simplex(mu,A,b,L,U,epsilon)
    if simp[:mu] != nothing; mu = simp[:mu]; end
    x = simp[:x]
    IN = simp[:IN]
    OUT = simp[:OUT]
    state = simp[:state]
    Ai = simp[:Ai]
    alpha[[ state.!=:IN; falses(m) ]] = x[state.!=:IN]
    IN=[IN; n+(1:m)]
    for in in IN
        bbar[in] = (in <= n ? 0 : b[in-n]) - dot(M[in,OUT], x[OUT])
    end
    extAi[IN[1:m],1:m] = Ai[1:m,1:m]
    Mi = [ zeros(n,n)           extAi            ;
           extAi'     -Ai'*C[IN[1:m],IN[1:m]]*Ai ]
    loop = 0
    while (loop < 10*(n+m))
        loop = loop + 1
        changes = []
        for in in IN
            alpha[in] = dot(Mi[in,IN], bbar[IN])
            beta[in] = dot(Mi[in,IN[IN.<=n]], mu[IN[IN.<=n]])
            if in <= n
                if beta[in] > epsilon
                    push!(changes,Change(:OUT,:DOWN,in,(L[in]-alpha[in])/beta[in]))
                elseif (U[in] < Inf) && (beta[in] < -epsilon)
                    push!(changes,Change(:OUT,:UP,in,(U[in]-alpha[in])/beta[in]))
                end
            end
        end
        for out in OUT
            gamma = dot(M[out,:], alpha)
            delta = -mu[out] + dot(M[out,:], beta)
            if (state[out] == :LOW) && (delta > epsilon)
                push!(changes,Change(:IN,:UP,out,-gamma/delta))
            elseif (state[out] == :HIGH) && (delta < -epsilon)
                push!(changes,Change(:IN,:DOWN,out,-gamma/delta))
            end
        end
        changes = sort(changes,by=x->-x.lambda)
        if length(changes) > 1 && changes[1].lambda == changes[2].lambda
            warn("more than one variable change could be done, will continue but note that the degenerate case is not properly handled")
        end
        lambdaEold = lambdaE
        if length(changes) > 0
            change = changes[1]
            changes = []
            lambdaE = max(0,change.lambda)
        else
            lambdaE = 0
        end
        x[IN[1:(length(IN)-m)]] = alpha[IN[1:(length(IN)-m)]] + lambdaE * beta[IN[1:(length(IN)-m)]]
        dEdLambdaE = beta[IN[1:(length(IN)-m)]]' * mu[IN[1:(length(IN)-m)]]
        if dEdLambdaE < epsilon
            E = x[1:n]' * mu[1:n]
            V = x[1:n]' * C * x[1:n]
        else
            a2 = 1/dEdLambdaE
            a1 = 2*(lambdaEold-a2*E)
            a0 = V-a1*E-a2*E*E
            E = E + (lambdaE-lambdaEold)*dEdLambdaE
            V = a0 + a1*E + a2*E*E
        end
        weights = vcat(weights,x')
        push!(lambda,lambdaE)
        push!(ret,E)
        push!(vol,sqrt(V*scale))
        if lambdaE < lambdaTarget
            return([weights,ret,vol,lambda])
        end
        if change.inout == :OUT
            out = change.variable
            alpha[out] = x[out]
            beta[out] = 0
            IN = IN[IN.!=out]
            if out <= n
                OUT = sort([ out; OUT ])
            end
            state[out] = (change.direction==:UP ? (:HIGH) : (:LOW))
            Mi = Mi - Mi[:,[out]] * Mi[[out],:] / Mi[out,out]
            bbar[IN] = bbar[IN] - M[IN,out] * x[out]
        else
            in = change.variable
            xi[IN] = Mi[IN,IN] * M[IN,[in]]
            xij = M[in,in] - dot(M[in,IN], xi[IN])
            Mi[IN,IN] = Mi[IN,IN] + xi[IN] * xi[IN]' / xij
            Mi[IN,in] = -xi[IN]/xij
            Mi[in,IN] = -xi[IN]/xij
            Mi[in,in] = 1/xij
            bbar[IN] = bbar[IN] + M[IN,in] * x[in]
            OUT = OUT[OUT.!=in]
            IN = sort([ in; IN ])
            state[in] = :IN
            bbar[in] = -dot(M[in,OUT], x[OUT])
        end
    end
end

mutable struct Form0
    V::Matrix{Float64}
    E::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    names::Vector{String}
    lhs::Matrix{Float64}
    op::Vector{Char}
    rhs::Vector{Float64}
end

function Base.display(m::Form0)
    if length(m.names) > 0
        println("Assets")
        Base.display(m.names)
    end
    println("Covariance matrix")
    Base.display(m.V)
    println("Expected returns")
    Base.display(m.E')
    println("Asset bounds")
    Base.display(vcat(m.lower', m.upper'))
    println("Additional constraints")
    Base.display(hcat(m.lhs, m.op, m.rhs))
end

function markowitz(E::Union{Vector,Matrix},V::Matrix;lower=0,upper=Inf,names=[])
    n = length(E)
    @assert size(V) == (n,n)
    if length(lower) == 1
        lower = lower*ones(n)
    end
    if length(upper) == 1
        upper = upper*ones(n)
    end
    Form0(V,vec(E),vec(lower),vec(upper),vec(names),
          Matrix{Float64}(0,n),Vector{Symbol}(0),Vector{Float64}(0))
end

function add_constraint(m::Form0,lhs::Vector,op::Char,rhs::Number)
    add_constraint(m,reshape(lhs,1,length(lhs)),op,rhs)
end

function add_constraint(m::Form0,lhs::Matrix,op::Char,rhs::Number)
    m.lhs = vcat(m.lhs, lhs)
    push!(m.op, op)
    push!(m.rhs, rhs)
    m
end

function unit_sum(m::Form0)
    add_constraint(m, ones(m.E), '=', 1)
end

struct Form2
    V::Matrix{Float64}
    E::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    A::Matrix{Float64}
    b::Vector{Float64}
end

function Base.display(m::Form2)
    println("Covariance matrix")
    Base.display(m.V)
    println("Expected returns")
    Base.display(m.E')
    println("Asset bounds")
    Base.display(vcat(m.lower', m.upper'))
    println("Equality constraints")
    Base.display(hcat(m.A, NaN*ones(length(m.b)), m.b))
end

function form2(m::Form0)
    n = length(m.E)
    nc = length(m.rhs)
    slack = sum(m.op.!='=')
    E = zeros(n+slack)
    E[1:n] = m.E
    V = zeros(n+slack,n+slack)
    V[1:n,1:n] = m.V
    lower = zeros(n+slack)
    lower[1:n] = m.lower
    upper = Inf*ones(n+slack)
    upper[1:n] = m.upper
    A = zeros(nc,n+slack)
    col = n+1
    for row in 1:nc
        A[row,1:n] = m.lhs[row,:]
        if m.op[row] != '='
            if m.op[row] == '<'
                A[row,col] = 1
            elseif m.op[row] == '>'
                A[row,col] = -1
            else
                error("invalid constraint type ",m.op[row])
            end
            col = col+1
        end
    end
    b = m.rhs
    Form2(V,E,lower,upper,A,b)
end

struct Frontier
    weights::Matrix{Float64}
    ret::Vector{Float64}
    vol::Vector{Float64}
    lambda::Vector{Float64}
    problem::Union{Form0,Form2}
end

function simplex(m::Form2)
    simplex(m.E,m.A,m.b,m.lower,m.upper)
end    

function frontier(m::Form2,form0=nothing)
    weights,ret,vol,lambda = cla(m.V,m.E,m.A,m.b,m.lower,m.upper)
    if form0 == nothing
        Frontier(weights,ret,vol,lambda,m)
    else
        Frontier(weights[:,1:length(form0.E)],ret,vol,lambda,form0)
    end
end    

function simplex(m::Form0)
    simplex(form2(m))
end    

function frontier(m::Form0)
    frontier(form2(m),m)
end    

function smooth(f,N=100)
    frontier = [ f.vol[1] f.ret[1] ]
    lengths = sqrt.((diff(f.ret)/(f.ret[1]-f.ret[end])).^2+(diff(f.vol)/(f.vol[1]-f.vol[end])).^2)
    nsegments = ceil.(N*lengths/sum(lengths))
    for i in 1:length(nsegments)
        w0 = f.weights[i,:]
        w1 = f.weights[i+1,:]
        if nsegments[i] > 0
            for p in (1:nsegments[i])/nsegments[i]
                w = w0 + p*(w1-w0)
                frontier = vcat(frontier, [ sqrt(w'*f.problem.V*w) f.problem.E'*w ])
            end
        end
    end
    frontier
end

struct Portfolio
    vol::Float64
    ret::Float64
    weights::Vector{Float64}
end

function optimal(f::Frontier,ret=nothing)
    if ret == nothing
        Portfolio(f.vol[end],f.ret[end],f.weights[end,:])
    else
        idx = findfirst(x-> x<ret, f.ret)
        if idx < 2
            error("the target return is not inside the frontier! min: ", f.ret[end], " max: ", f.ret[1])
        else
            p = (f.ret[idx-1]-ret)/(f.ret[idx-1]-f.ret[idx])
            w = f.weights[idx-1,:] + p*(f.weights[idx,:]-f.weights[idx-1,:])
            Portfolio(sqrt(w'*f.problem.V*w), f.problem.E'*w, w)
        end
    end
end

end
