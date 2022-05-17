#Explorations that didn't result in anything
#Other code should be here but is somewhere in prev commits
function route_center(route::Route, vars::VRP)
    cx, cy = 0, 0
    for i = 1:route.seqlen
        cx += vars.positions[route[i]][1]
        cy += vars.positions[route[i]][2]
    end
    [cx / route.seqlen, cy / route.seqlen]
end
function average_center_distance(route::Route, center::Vector, vars::VRP)
    dist = 0
    for i = 1:route.seqlen
        dist += euc_dist(vars.positions[route[i]], center)
    end
    dist / route.seqlen
end
function get_closest_routes(sol::Solution, vars::VRP)
    centers = map(i -> route_center(sol[i], vars), 1:vars.vehicles)
    centerDists = map(i -> (i, average_center_distance(sol[i], centers[i], vars)), 1:vars.vehicles)
    sort!(centerDists, by=i -> i[2])
    randRoute = centerDists[end][1]
    println("center is ", randRoute)
    routeDist = map(i -> (i, euc_dist(centers[randRoute], centers[i])), 1:vars.vehicles)
    sort!(routeDist, by=i -> i[2])
    closeRoutes = map(i -> i[1], routeDist)[1:4]
    println("closest routes are ", closeRoutes)
    remove_c = []
    for r in closeRoutes
        x = sol[r].seq[2:sol[r].seqlen]
        append!(remove_c, x)
    end
    remove_c
end

function kmean_lp(vars::VRP)
    assigs = Clustering.kmeans(transpose(vars.positions), vars.vehicles).assignments
    mtx = zeros(vars.vehicles, vars.customers)
    for (i, a) in enumerate(assigs)
        mtx[a, i] = 1
    end
    mtx
end
function comp_diff(existing_mat::Matrix, vrp::VRP)
    # kmean_mat  = kmean_lp(vrp)x

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)
    fix(x[1, 1], 1, force=true)
    min_expr = @expression(model, 0)
    for i in V
        for j in C
            if existing_mat[i, j] == 1
                min_expr += x[i, j]
            end
        end
    end
    @objective(model, Min, min_expr)
    optimize!(model)
    lpvals = zeros(Int8, vrp.vehicles, vrp.customers)
    changes = 0
    for j = 1:vrp.customers
        for i = 1:vrp.vehicles
            lpvals[i, j] = Int(round(value(mtrx[i, j])))
        end
    end
    nsol = solverToSol(vrp, lpvals)
    sol2Opt(nsol, vrp)
    println("before local search ", nsol.objective)
    lnsol = localSearch(nsol, vrp)
    println("final objective is ", lnsol.objective)
    lnsol
end
function n_rand_groups(vars::VRP, g::Number)
    cust = shuffle(1:vars.customers)
    take = Int(round(vars.customers / g))
    ls = []
    i = 1
    while i <= length(cust)
        if i + take <= length(cust)
            push!(ls, cust[i:i+take])
            i += (take + 1)
        else
            push!(ls, cust[i:end])
            break
        end
    end
    ls
end
function fulldestroy(sol::Solution, vars::VRP)

    # take = Int(round(0.55 * vars.customers))

    # group = sample(1:vars.customers,take,replace=false)
    # group = map(x -> x[1] ,vars.node_demand[1:take])     
    # group = get_closest_routes(sol, vars)

    # destroy_tsp(sol, vars, group)
    # println("new sol ", sol.objective)
    # vis_output(sol, vars)
    sol1mat = sol_to_mat(sol, vars)
    sol2mat = sol_to_mat(comp_diff(sol1mat, vars), vars)
    sol3mat = sol_to_mat(comp_diff(sol1mat + sol2mat, vars), vars)
    s = comp_diff(sol1mat + sol2mat + sol3mat, vars)
    vis_output(s, vars)
    s
end
function kmean_sol(fl::String)
    vars = read_input(fl)
    kmean_based(vars)
end
"""
K Mean Clustering with v (num vehicles) clusters with v random customers as 
centers
Does not respect load
"""
function nn_to_lpmatrix(vars::VRP)
    solMat = zeros(Bool, vars.vehicles, vars.customers)
    centers = sample(1:vars.customers, vars.vehicles, replace=false)
    foreach(t -> solMat[t[1], t[2]] = 1, enumerate(centers))
    cust = filter(i -> !(i in centers), 1:vars.customers)
    for c in cust
        argmin, min = 1, Inf
        for (i, center) in enumerate(centers)
            if vars.distance_m[c, center] < argmin
                argmin, min = i, vars.distance_m[c, center]
            end
        end
        solMat[argmin, c] = 1
    end
    solMat
end
function nn_based(sol::Solution, vrp::VRP)
    nnMat = nn_to_lpmatrix(vrp) .+ 1
    remove_c = sample(1:vrp.customers, Int(round(1 * vrp.customers)), replace=false)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)
    cust_routes = get_customer_routes(sol)
    remove_mat = zeros(vrp.vehicles, vrp.customers)
    obj_expr = @expression(model, 0)
    for c in remove_c
        curr_route = cust_routes[c]
        obj_expr += @expression(model, x[curr_route, c])
    end
    @objective(model, Min, sum(nnMat .- x))
    optimize!(model)
    lpvals = zeros(vrp.vehicles, vrp.customers)
    for i = 1:vrp.vehicles
        for j = 1:vrp.customers
            lpvals[i, j] = value(mtrx[i, j])
        end
    end
    nsol = solverToSol(vrp, lpvals)
    sol2Opt(nsol, vrp)
    lnsol = localSearch(nsol, vrp)
    sol2Opt(lnsol, vrp)
    lnsol
end


function runExp(sol :: Solution,vars::VRP)
    mean(
        [localSearch(sol,vars) for i=1:3]
    )
end

"""
Returns a lm fraction of the nearest neighbors of a random node
"""
function nn_cluster(vars::VRP, lm::Float64)
    numtake = Int(round(lm * vars.customers))
    cent = rand(1:vars.customers)
    cust_d = [(i, vars.distance_m[cent, i]) for i in 1:vars.customers]
    sort!(cust_d, by=x -> x[2])
    map(x -> x[1], cust_d[1:numtake])
end
function destroy_tsp(sol::Solution, vrp::VRP, remove_c)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)
    min_expr = @expression(model, 0)
    cust_routes = get_customer_routes(sol)
    for c in C
        curr_route = cust_routes[c]
        if !(c in remove_c)
            @constraint(model, x[curr_route, c] == 1)
        else
            min_expr += @expression(model, x[curr_route, c])
        end
    end
    @objective(model, Min, min_expr)
    optimize!(model)
    lpvals = zeros(Int8, vrp.vehicles, vrp.customers)
    for j = 1:vrp.customers
        for i = 1:vrp.vehicles
            lpvals[i, j] = Int(round(value(mtrx[i, j])))
        end
    end
    solverToSol(vrp, lpvals)
end
function kmean_based(vrp::VRP)
    kmean_mat = kmean_lp(vrp)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)

    @objective(model, Min, sum(kmean_mat - x))
    optimize!(model)
    lpvals = zeros(Int8, vrp.vehicles, vrp.customers)
    changes = 0
    for j = 1:vrp.customers
        for i = 1:vrp.vehicles
            lpvals[i, j] = Int(round(value(mtrx[i, j])))
        end
    end
    nsol = solverToSol(vrp, lpvals)
    sol2Opt(nsol, vrp)
    println("before local search ", nsol.objective)
    lnsol = localSearch(nsol, vrp)
    println("final objective is ", lnsol.objective)
    lnsol
end

function read_goog_vrp(fl::String, vars::VRP)
    try
        open(fl, "r") do io
            lines = readlines(io)
            lines = [split(i) for i in lines]
            objective = parse(Float64, lines[1][1])
            # println("i is ",lines[2][2:end-1]," length ",length(lines[2][2:end-1]))
            routes = [map(x -> parse(Int64, x), i[2:end-1]) for i in lines[2:end]]
            # println("routes are ",routes)
            route_obj = []
            for route in routes
                load = 0
                for cust in route
                    load += vars.demand[cust]
                end
                ln = length(route)
                push!(route_obj, Route(resize!(route, vars.customers), ln, load))
            end
            sol = Solution(route_obj, objective,Dict())
            sol2Opt(sol,vars)
            sol
        end
    catch
        error("could not read test vrp file")
    end
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
function sorted_distance(dist_m::Matrix)
    l = size(dist_m)[1]
    sort_m = Matrix{Tuple{Number,Number}}(undef, l, l)
    for i = 1:l
        tmp = map(j -> (dist_m[i, j], j), 1:l)
        sort_m[i,:] = sort(tmp,by=i->i[1])
    end
    sort_m
end
function getSecondNode(dist::Weights, vars::VRP)::Int64
    vars.node_demand[sample(1:vars.customers, dist)][1]
end
function get_dist(vars::VRP)::Vector{Weights}
    ls = []
    v = mean(vars.demand) / 2
    for i = 1:vars.customers
        pdfvals = [abs(i - j) <= v ? 1 : 0 for j in 1:vars.customers]
        pdfvals[i] = 0
        push!(ls, Weights(StatsBase.LinearAlgebra.normalize(pdfvals)))
    end
    ls
end