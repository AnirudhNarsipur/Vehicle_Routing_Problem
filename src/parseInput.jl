include("./datastructs.jl")
include("./utils.jl")

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
function read_test_vrp(fl::String, vars::VRP)
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

