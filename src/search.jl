# module TRP
using StatsBase
include("./parseInput.jl")
"""
1. take route[1] to route[i-1] and add them in order to new_route
2. take route[i] to route[k] and add them in reverse order to new_route
3. take route[k+1] to end and add them in order to new_route
"""
function opt2Swap(route::Route, i::Number, k::Number)
    route.seq[i:k] = reverse!(route.seq[i:k])
end

"""
Given route and vars makes route 2-opt
"""
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
"""
Makes all routes in sol 2opt
"""
function sol2Opt(sol::Solution, vars::VRP)
    for route in sol.routes
        complete2optSwap(route, vars)
    end
    sol.objective = recalc_obj_val(sol, vars)
end

function swapNodes(sol::Solution, vars::VRP, frloc::Number, fposloc::Number, srloc::Number, sposloc::Number)
    fdemand = vars.demand[sol.routes[frloc].seq[fposloc]]
    sdemand = vars.demand[sol.routes[srloc].seq[sposloc]]
    sol.routes[frloc].load = sol.routes[frloc].load - fdemand + sdemand
    sol.routes[srloc].load = sol.routes[srloc].load - sdemand + fdemand
    tmp = sol.routes[frloc].seq[fposloc]
    sol.routes[frloc].seq[fposloc] = sol.routes[srloc].seq[sposloc]
    sol.routes[srloc].seq[sposloc] = tmp
    nothing
end
function calc_route_load(route::Route, vars::VRP)::Number
    d = 0
    for i = 1:route.seqlen
        d += vars.demand[route.seq[i]]
    end
    d
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
function nearNodeSwap(sol::Solution, vars::VRP)
    route_load_change = (rn, f, s) -> sol.routes[rn].load - f + s
    fnode = rand(1:vars.vehicles)
    snode = vars.sorted_d[fnode, rand(2:9)][2]
    frloc, fposloc = findCloc(sol, fnode)
    srloc, sposloc = findCloc(sol, snode)
    if frloc == srloc
        return Inf, () -> nothing
    end
    first = sol.routes[frloc].seq[fposloc]
    second = sol.routes[srloc].seq[sposloc]
    fdemand = vars.demand[first]
    sdemand = vars.demand[second]
    if route_load_change(frloc, fdemand, sdemand) > vars.capacity || route_load_change(srloc, sdemand, fdemand) > vars.capacity
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
function randomNodeSwap(sol::Solution, vars::VRP)
    route_load_change = (rn, f, s) -> sol.routes[rn].load - f + s
    frloc, srloc = rand(1:vars.vehicles, 2)
    if frloc == srloc
        return Inf, () -> nothing
    end
    fposloc = rand(1:sol.routes[frloc].seqlen)
    sposloc = rand(1:sol.routes[srloc].seqlen)
    first = sol.routes[frloc].seq[fposloc]
    second = sol.routes[srloc].seq[sposloc]
    fdemand = vars.demand[first]
    sdemand = vars.demand[second]
    if route_load_change(frloc, fdemand, sdemand) > vars.capacity || route_load_change(srloc, sdemand, fdemand) > vars.capacity
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

function localSearch(sol::Solution, vars::VRP)

    numItr = 5e+4 * vars.customers
    sol.objective = recalc_obj_val(sol, vars)
    Temperature = abs(sol.objective / log(MathConstants.e, 0.97))
    numRuns = 1
    annealing_cycle = vars.customers * (vars.customers - 1) / 2
    best_sol = deepcopy(sol)
    prob_func = (ns) -> exp((sol.objective - ns) / Temperature)
    for _ = 1:numItr
        numRuns += 1
        if numRuns > annealing_cycle
            numRuns = 1
            Temperature *= 0.95
        end
        nscore, changeFunc = randomNodeSwap(sol, vars)
        if nscore < best_sol.objective
            changeFunc()
            best_sol = deepcopy(sol)
        elseif rand() <= prob_func(nscore)
            changeFunc()
        end
    end
    sol.objective = recalc_obj_val(sol, vars)
    newsol = deepcopy(best_sol)
    sol2Opt(newsol,vars)
    newsol
end
function nn_heur(vars::VRP)
    groups = []
    cust = [1:vars.customers...]
    while length(cust) != 0
        ind = rand(1:length(cust))
        c = popat!(cust, ind)
        curr_group = [c]
        curr_load = vars.demand[c]
        cust_d = [(i, vars.distance_m[c, i]) for i in cust]
        sort!(cust_d, by=x -> x[2])
        i = 1
        while i <= length(cust_d)
            if curr_load + vars.demand[cust_d[i][1]] > vars.capacity
                break
            end
            push!(curr_group, cust_d[i][1])
            filter!(x -> x != cust_d[i][1], cust)
            curr_load += vars.demand[cust_d[i][1]]
            i += 1
        end
        push!(groups, Route(curr_group, length(curr_group), curr_load))
    end
    Solution(groups, 0)
end
function nn_method(vars::VRP)
    sol = nn_heur(vars)
    sol.objective = recalc_obj_val(sol, vars)
    sol2Opt(sol, vars)
    vis_output(sol, vars)
    sol
end
function routeLengths(sol :: Solution)
    print("Route length are ")
    for route in sol.routes
        print(route.seqlen," ")
    end
    println("")
end

function destroy_tsp(sol::Solution, vrp::VRP)
    routeLengths(sol)
    remove_c = sample(1:vrp.customers, Int(round(0.5* vrp.customers)), replace=false)
    model = Model(HiGHS.Optimizer)
    # set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)
    cust_routes = get_customer_routes(sol)
    min_expr = @expression(model, 0)
    for c in C
        curr_route = cust_routes[c]
        if !(c in remove_c)
            @constraint(model, x[curr_route, c] == 1)
        else
            min_expr += @expression(model, x[curr_route, c])
        end
    end
    @objective(model,Min,min_expr)
    optimize!(model)
    lpvals = zeros(Int8,vrp.vehicles, vrp.customers)
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
    println(" num changes are : " ,changes," out of ",length(remove_c))
    nsol = solverToSol(vrp, lpvals)
    sol2Opt(nsol, vrp)
    lnsol = localSearch(nsol, vrp)
    sol2Opt(lnsol, vrp)
    routeLengths(lnsol)
    lnsol
end

function fulldestroy(fl :: String)
    vars,sol = initStuff(fl)
    println("starting sol ",sol.objective)
    ns = destroy_tsp(sol,vars)
    println("new sol ",ns.objective)
    vis_output(ns, vars)
end

function nn_to_lpmatrix(vars::VRP)
    sol = nn_heur(vars)
    solMat = zeros(Bool, vars.vehicles, vars.customers)
    for r = 1:min(vars.vehicles, length(sol.routes))
        for cindex = 1:sol.routes[r].seqlen
            solMat[r, sol.routes[r].seq[cindex]] = 1
        end
    end
    solMat
end
function nn_based(sol::Solution, vrp::VRP)
    nnMat = nn_to_lpmatrix(vrp)
    remove_c = sample(1:vrp.customers, Int(round(0.35 * vrp.customers)), replace=false)
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
    # bad_cust_r = cust_routes[remove_c[1]]
    # @constraint(model,x[bad_cust_r,remove_c[1]] == 0)
    @objective(model, Max, obj_expr)
    optimize!(model)
    lpvals = zeros(vrp.vehicles, vrp.customers)
    for i = 1:vrp.vehicles
        for j = 1:vrp.customers
            lpvals[i, j] = value(mtrx[i, j])
        end
    end
    println("available sols : ", result_count(model))
    nsol = solverToSol(vrp, lpvals)
    sol2Opt(nsol, vrp)
    lnsol = localSearch(nsol, vrp)
    sol2Opt(lnsol, vrp)
    lnsol
end
function initStuff(fl::String)
    vars = read_input(fl)
    sol = getInitialSol(vars)
    sol2Opt(sol, vars)
    numRuns = 1
    for i = 1:numRuns
        tmp = deepcopy(sol)
        localSearch(tmp, vars)
        sol2Opt(tmp, vars)
        tmp.objective = recalc_obj_val(tmp, vars)
        if tmp.objective < sol.objective
            sol = tmp
        end
    end
    vars, sol
end

function mn(fl::String)
    start_time = Base.Libc.time()
    vars = read_input(fl)
    sol = getInitialSol(vars)
    sol2Opt(sol, vars)
    org_d = sol.objective
    numRuns = 2
    for i = 1:numRuns
        tmp = deepcopy(sol)
        localSearch(tmp, vars)
        sol2Opt(tmp, vars)
        tmp.objective = recalc_obj_val(tmp, vars)
        if tmp.objective < sol.objective
            sol = tmp
        end
    end
    f_d = sol.objective
    end_time = Base.Libc.time()
    println("initial obj is ", org_d, " final is ", f_d, " improv is ", round(((org_d - f_d) / org_d) * 100, digits=2), "%")
    vis_output(sol, vars)
    get_output(fl, end_time - start_time, sol, vars)
end

# function __init__()
#     mn(ARGS[1])
# end
# end