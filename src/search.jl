# module TRP
using StatsBase
using Random
include("./parseInput.jl")

function swapNodes(sol::Solution, vars::VRP, frloc::Number, fposloc::Number, srloc::Number, sposloc::Number)
    fdemand = vars.demand[sol.routes[frloc].seq[fposloc]]
    sdemand = vars.demand[sol.routes[srloc].seq[sposloc]]
    sol.routes[frloc].load = sol.routes[frloc].load - fdemand + sdemand
    sol.routes[srloc].load = sol.routes[srloc].load - sdemand + fdemand
    sol.nodeloc[sol.routes[frloc].seq[fposloc]] = (srloc, sposloc)
    sol.nodeloc[sol.routes[srloc].seq[sposloc]] = (frloc, fposloc)
    tmp = sol.routes[frloc].seq[fposloc]
    sol.routes[frloc].seq[fposloc] = sol.routes[srloc].seq[sposloc]
    sol.routes[srloc].seq[sposloc] = tmp
    nothing
end

function pointdistance(vars::VRP, sol::Solution, rloc::Number, rpos::Number)
    c = sol.routes[rloc].seq[rpos]
    if sol.routes[rloc].seqlen <= 1
        route_distance(sol.routes[rloc], vars)
    elseif rpos == 1
        vars.depot_distance[c] + vars.distance_m[c, sol.routes[rloc].seq[rpos+1]]
    elseif rpos == sol.routes[rloc].seqlen
        vars.distance_m[c, sol.routes[rloc].seq[rpos-1]] + vars.depot_distance[c]
    else
        vars.distance_m[c, sol.routes[rloc].seq[rpos-1]] + vars.distance_m[c, sol.routes[rloc].seq[rpos+1]]
    end
end

function getrandomNodeSwap(sol::Solution, vars::VRP)
    route_load_change = (rn, f, s) -> sol.routes[rn].load - f + s
    function randomNodeSwap(wghts::Weights)
        first, second = sample(1:vars.customers, wghts, 2, replace=false)
        frloc, fposloc = sol.nodeloc[first]
        srloc, sposloc = sol.nodeloc[second]
        if frloc == srloc
            return Inf, () -> nothing
        end
        if route_load_change(frloc, vars.demand[first], vars.demand[second]) > vars.capacity || route_load_change(srloc, vars.demand[second], vars.demand[first]) > vars.capacity
            return Inf, () -> nothing
        end
        oldp = pointdistance(vars, sol, frloc, fposloc) + pointdistance(vars, sol, srloc, sposloc)
        swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        nd = pointdistance(vars, sol, frloc, fposloc) + pointdistance(vars, sol, srloc, sposloc)
        swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        newobj = sol.objective - oldp + nd
        function commitChange()
            swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
            sol.objective = newobj
        end
        return newobj, commitChange
    end
    randomNodeSwap
end
function localSearch(sol::Solution, vars::VRP)::Solution
    b = 0
    t = 0
    numItr = 5e+4 * vars.customers
    sol.objective = recalc_obj_val(sol, vars)
    Temperature = abs(sol.objective / log(MathConstants.e, 0.97))
    numRuns = 1
    annealing_cycle = vars.customers * (vars.customers - 1) / 2
    best_sol = deepcopy(sol)
    prob_func = (ns) -> exp((sol.objective - ns) / Temperature)
    randomNodeSwap = getrandomNodeSwap(sol, vars)
    # solutionCheck(sol,vars)
    wghts = Weights([1 / vars.customers for _ in 1:vars.customers])
    route_load_change = (rn, f, s) -> sol.routes[rn].load - f + s
    delta  =0 
    for i = 1:numItr
        numRuns += 1
        if numRuns > annealing_cycle
            numRuns = 1
            Temperature *= 0.95
        end
        #start_time
        first, second = sample(1:vars.customers, wghts, 2, replace=false)
        frloc, fposloc = sol.nodeloc[first]
        srloc, sposloc = sol.nodeloc[second]
        if frloc == srloc
            continue
        end
        if route_load_change(frloc, vars.demand[first], vars.demand[second]) > vars.capacity || route_load_change(srloc, vars.demand[second], vars.demand[first]) > vars.capacity
            continue
        end
        oldp = pointdistance(vars, sol, frloc, fposloc) + pointdistance(vars, sol, srloc, sposloc)
        swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        nd = pointdistance(vars, sol, frloc, fposloc) + pointdistance(vars, sol, srloc, sposloc)
        swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        newobj = sol.objective - oldp + nd
        function commitChange()
            swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
            sol.objective = newobj
        end
        #END END
        nscore, changeFunc = newobj, commitChange
        if nscore < best_sol.objective
            delta+=abs(vars.demand[first]-vars.demand[second])
            changeFunc()
            best_sol = deepcopy(sol)
            b += 1
            t += 1
        elseif rand() <= prob_func(nscore)
            changeFunc()
            delta+=abs(vars.demand[first]-vars.demand[second])
            t += 1
        end
    end
    println("delta is ",delta/numItr)
    sol2Opt(best_sol, vars)
    # solutionCheck(best_sol,vars)
    best_sol
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
    changes = 0
    for j = 1:vrp.customers
        curr_route = cust_routes[j]
        for i = 1:vrp.vehicles
            lpvals[i, j] = Int(round(value(mtrx[i, j])))
            if j in remove_c && curr_route == i && lpvals[i, j] == 0
                changes += 1
            end
        end
    end
    println(" num changes are : ", changes, " out of ", length(remove_c))
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
function fulldestroy(fl::String)
    vars = read_input(fl)
    sol = getInitialSol(vars)
    sol = localSearch(sol, vars)
    println("starting sol ", sol.objective)
    # groups = n_rand_groups(vars,3)
    take = Int(round(0.95 * vars.customers))
    # group = sample(1:vars.customers,take,replace=false)
    # group = map(x -> x[1] ,vars.node_demand[1:take])
    group = get_closest_routes(sol, vars)
    vis_output(sol, vars)
    return nothing
    println("rebuilding with close routes")
    sol = destroy_tsp(sol, vars, group)

    println("new sol ", sol.objective)
    vis_output(sol, vars)
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
    closeRoutes = map(i -> i[1], routeDist)[1:3]
    println("closest routes are ", closeRoutes)
    remove_c = []
    for r in closeRoutes
        x = sol[r].seq[2:sol[r].seqlen]
        append!(remove_c, x)
    end
    remove_c
end



function mn(fl::String)
    start_time = Base.Libc.time()
    vars = read_input(fl)
    sol = getInitialSol(vars)
    sol2Opt(sol, vars)
    org_d = sol.objective
    numRuns = 1
    bestsol = sol
    for i = 1:numRuns
        tmp = deepcopy(sol)
        # tmp = deepcopy(destroy_tsp(tmp,vars))
        tmp = localSearch(tmp, vars)
        # println("objective is ",tmp.objective)
        tmp.objective = recalc_obj_val(tmp, vars)
        if tmp.objective < bestsol.objective
            bestsol = deepcopy(tmp)
        end
    end
    f_d = bestsol.objective
    end_time = Base.Libc.time()
    println("initial obj is ", org_d, " final is ", f_d, " improv is ", round(((org_d - f_d) / org_d) * 100, digits=2), "%")
    vis_output(bestsol, vars)
    get_output(fl, end_time - start_time, bestsol, vars)

end

# function __init__()
#     mn(ARGS[1])
# end
# end