using JuMP
import HiGHS
struct VRP
    customers::Number
    vehicles::Number
    capacity::Number
    demand::Vector
    depot_location::Vector
    positions::Vector
    distance_m::Matrix
    depot_distance::Vector
    # Node -> Position in node_pos
    node_pos::Dict{Number,Number}
    # demand, customer node number 
    node_demand::Vector{Tuple{Number,Number}}
end
mutable struct Route
    seq::Vector{Number}
    seqlen::Number
    load::Number
end

mutable struct Solution
    routes::Vector{Route}
    objective::Number
end

abstract type Option end
# Return a value
struct Some{T} <: Option
    value::T
end
# None - silent positive
struct None <: Option end
#Bad - failure
struct Bad <: Option end
# Write the above python function into julia
function euc_dist(a::Vector, b::Vector)
    sqrt(sum((a - b) .^ 2))
end

function create_distance_matrix(positions::Vector)
    c = length(positions)
    dist_m = zeros(c, c)
    for i = 1:c
        for j = 1:c
            dist_m[i, j] = euc_dist(positions[i], positions[j])
        end
    end
    dist_m
end
function depot_distance(positions::Vector, depot_pos::Vector)
    c = length(positions)
    depot_m = zeros(c)
    for i = 1:c
        depot_m[i] = euc_dist(positions[i], depot_pos)
    end
    depot_m
end
function getSortedDemand(demand::Vector)
    demand_index = map(i -> (i, demand[i]), 1:length(demand))
    sort!(demand_index)
    node_dict = Dict()
    for (index, elem) in enumerate(demand_index)
        node_dict[elem[1]] = index
    end
    return node_dict, demand_index
end

function read_input(fl::String)
    try
        open(fl, "r") do io
            lines = readlines(io)
            lines = [split(i) for i in lines]
            customers, vehicles, capacity = [parse(Int64, i) for i in lines[1]]
            depot_location = [floor(parse(Float64, i)) for i in lines[2][2:3]]
            demand = []
            positions = []
            for i = 3:length(lines)
                if length(lines[i]) == 0
                    break
                end
                push!(demand, parse(Float64, lines[i][1]))
                push!(positions, [floor(parse(Float64, lines[i][j])) for j = 2:3])
            end
            node_pos, node_demand = getSortedDemand(demand)
            VRP(customers - 1, vehicles, capacity, demand, depot_location,
                positions, create_distance_matrix(positions), depot_distance(positions, depot_location),
                node_pos, node_demand)
        end
    catch
        error("could not read file!")
    end
end

function slvr(fl::String)
    vrp = read_input(fl)
    model = Model(HiGHS.Optimizer)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)
    optimize!(model)
    lpvals = zeros(vrp.vehicles, vrp.customers)
    for i = 1:vrp.vehicles
        for j = 1:vrp.customers
            lpvals[i, j] = value(mtrx[i, j])
        end
    end
    lpvals
end
function read_test_vrp(fl :: String,vars::VRP)
    objective = 0
    route_mtx = slvr(fl)
    routes = []
    for i = 1:vars.vehicles
        ls = []
        for j = 1:vars.customers
            if route_mtx[i, j] == 1
                push!(ls, j)
            end
        end
        push!(routes, ls)
    end
    route_obj = []
    for route in routes
        load = 0
        for cust in route
            load += vars.demand[cust]
        end
        ln = length(route)
        push!(route_obj, Route(resize!(route, vars.customers), ln, load))
    end
    sol = Solution(route_obj, objective)
    sol.objective = recalc_obj_val(sol, vars)
    sol
end

function route_distance(route::Route, vars::VRP)
    distances = vars.depot_distance[route.seq[1]]
    for i = 1:route.seqlen
        if i == route.seqlen
            distances += vars.depot_distance[route.seq[i]]
        else
            distances += vars.distance_m[route.seq[i], route.seq[i+1]]
        end
    end
    distances
end
function recalc_obj_val(sol::Solution, vars::VRP)
    sum(map(r -> route_distance(r, vars), sol.routes))
end

#  1. take route[1] to route[i-1] and add them in order to new_route
#  2. take route[i] to route[k] and add them in reverse order to new_route
#  3. take route[k+1] to end and add them in order to new_route
# new_route
function opt2Swap(route::Route, i::Number, k::Number)
    route.seq[i:k] = reverse!(route.seq[i:k])
end


function complete2optSwap(route::Route, vars::VRP)
    improve_made = true
    while improve_made
        improve_made = false
        best_dist = route_distance(route, vars)
        for i = 1:route.seqlen
            for k = 1:route.seqlen
                opt2Swap(route, i, k)
                new_dist = route_distance(route, vars)
                if (new_dist < best_dist)
                    improve_made = true
                else
                    opt2Swap(route, i, k)
                end
            end
        end
    end
end

function sol2Opt(sol::Solution, vars::VRP)
    for route in sol.routes
        complete2optSwap(route, vars)
    end
    sol.objective  = recalc_obj_val(sol, vars)
end


function vis_output(sol::Solution, vars::VRP)
    lines = [string(sol.objective) * " 0"]
    for route in sol.routes
        out = "0"
        for i = 1:route.seqlen
            out *= " "
            out *= string(route.seq[i])
        end
        out *= " 0"
        push!(lines, out)
    end
    open("vistest.vrp", "w") do io
        for line in lines
            write(io, line * "\n")
        end
    end
end
function findCloc(sol::Solution, customer::Number)
    for (i, route) in enumerate(sol.routes)
        for j = 1:route.seqlen
            if route.seq[j] == customer
                return (i, j)
            end
        end
    end
    error("Could not find customer!")
end
function swapNodes(sol::Solution,vars::VRP, frloc::Number, fposloc::Number, srloc::Number, sposloc::Number)
    fdemand = vars.demand[sol.routes[frloc].seq[fposloc]]
    sdemand = vars.demand[sol.routes[srloc].seq[sposloc]]
    sol.routes[frloc].load = sol.routes[frloc].load - fdemand + sdemand
    sol.routes[srloc].load = sol.routes[srloc].load - sdemand + fdemand
    tmp = sol.routes[frloc].seq[fposloc]
    sol.routes[frloc].seq[fposloc] = sol.routes[srloc].seq[sposloc]
    sol.routes[srloc].seq[sposloc] = tmp
end
function calc_route_load(route::Route, vars::VRP)
    d = 0
    for i = 1:route.seqlen
        d += vars.demand[route.seq[i]]
    end
    if d > vars.capacity
        error("Demand exceeded")
    end
    d
end
function checkSwapNodes(sol::Solution, vars::VRP)
    first, second = rand(1:vars.customers, 2)
    if first == second
        return false
    end
    frloc, fposloc = findCloc(sol, first)
    srloc, sposloc = findCloc(sol, second)
    fdemand = vars.demand[first]
    sdemand = vars.demand[second]
    if frloc == srloc
        return false
    end
    if (sol.routes[frloc].load - fdemand + sdemand) > vars.capacity || (sol.routes[srloc].load - sdemand + fdemand) > vars.capacity
        return false
    end
    swapNodes(sol,vars, frloc, fposloc, srloc, sposloc)
    newdistance = recalc_obj_val(sol,vars)
    swapNodes(sol, vars,frloc, fposloc, srloc, sposloc)
    # @assert sol.routes[frloc].seq[fposloc] == first
    # @assert sol.routes[srloc].seq[sposloc] == second
    return newdistance,frloc,fposloc,srloc,sposloc
end

function localSearch(sol::Solution, vars::VRP)
    b = 0
    w = 0
    t =  1e+6
    sol.objective = recalc_obj_val(sol,vars)
    T = 10
    n=1
    mx = vars.customers * (vars.customers-1)/2
    best_sol = deepcopy(sol)
    prob_func = (ns) -> exp((best_sol.objective-ns)/T)
    for i = 1:t
        n+=1
        if n > mx
            n = 1
            T*=0.95
        end
        res = checkSwapNodes(sol, vars)
        if res isa Bool
            continue
        else
            nscore = res[1]
            f1,f2,s1,s2 = res[2:end]
            if nscore < best_sol.objective
                swapNodes(sol,vars,f1,f2,s1,s2)
                sol.objective = nscore
                best_sol=sol
                b+=1
            elseif rand() <= prob_func(nscore)
                best_sol = deepcopy(sol)
                sol.objective = nscore
                swapNodes(sol,vars,f1,f2,s1,s2)
                w+=1
            end
        end
    end
    sol = best_sol
    println("Made " * string(b) * " better swaps " * string(w) * " worse swaps out of " * string(t) * " swaps")
end
function mn()
    fl = "input/135_7_1.vrp"
    vars = read_input(fl)
    sol = read_test_vrp(fl, vars)
    sol2Opt(sol, vars)
    println("Distance before local search " * string(sol.objective))
    localSearch(sol, vars)
    sol2Opt(sol, vars)
    sol.objective = recalc_obj_val(sol, vars)
    println("Distance after local search " * string(sol.objective))

    vis_output(sol, vars)
end