include("./datastructs.jl")
include("./utils.jl")

"""
Return a  vector of distances from nodes to depot
"""
function depot_distance(positions::Matrix{Float64}, depot_pos::Vector)

    c = size(positions)[1]
    depot_m = zeros(c)
    for i = 1:c
        depot_m[i] = euc_dist(positions[i,:], depot_pos)
    end
    depot_m
end
function read_input(fl::String)
    try
        open(fl, "r") do io
            lines = readlines(io)
            lines = [split(i) for i in lines]
            customers, vehicles, capacity = [parse(Int64, i) for i in lines[1]]
            customers-=1
            depot_location = [parse(Float64, i) for i in lines[2][2:3]]
            demand = []
            positions = []
            positions = zeros(customers,2)
            for i = 3:length(lines)
                if length(lines[i]) == 0
                    break
                end
                push!(demand, parse(Float64, lines[i][1]))
                positions[i-2,:] = [parse(Float64, lines[i][j]) for j = 2:3]
            end
            dist_m = create_distance_matrix(positions)
            VRP(customers, vehicles, capacity, demand, depot_location,
                positions, dist_m, depot_distance(positions, depot_location))
        end
    catch
        error("could not read file!")
    end
end
"""
Get Initial feasible solution to VRP from LP solver
Returns matrix of VxC binary values
"""
function lpSolver(vrp::VRP)

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    fix(x[1,1],1,force=true)
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
"""
Given VRP struct and LP result converts to Solution datastructure
"""
function solverToSol(vars::VRP, route_mtx)::Solution
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
    cust_route_dict = Dict()
    for i=1:length(routes)
        route = routes[i]
        load = 0
        for j=1:length(route)
            cust = route[j]
            load += vars.demand[cust]
            cust_route_dict[cust] = (i,j)
        end
        ln = length(route)
        push!(route_obj, Route(resize!(route, vars.customers), ln, load))
    end
    sol = Solution(route_obj, 0,cust_route_dict)
    sol.objective = recalc_obj_val(sol, vars)
    sol2Opt(sol,vars)
    # solutionCheck(sol,vars)
    sol
end
"""
Return a initial feasible solution
"""
function getInitialSol(vars::VRP)
    route_mtx = lpSolver(vars)
    solverToSol(vars, route_mtx)
end
function get_output(fl::String, time, sol::Solution)
    kv_json = (k, v) -> join(['"', k, """": """, '"', v, '"'])
    solo = []
    for route in sol.routes
        push!(solo, "0")
        for i = 1:route.seqlen
            push!(solo, string(route.seq[i]))
        end
        push!(solo, "0")
    end
    out = join([
        "{",
        kv_json("Instance", basename(fl)), ",",
        kv_json("Time", round(time, digits=2)), ",",
        kv_json("Result", string(round(sol.objective, digits=2))), ",",
        kv_json("Solution", join(solo, " ")),
        "}"
    ])
    println(out)
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