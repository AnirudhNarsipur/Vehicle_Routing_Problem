# module TRP
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

"""
Given a solution and a customer
return (routenumber,number in route)
"""
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
function swapNodes(sol::Solution, vars::VRP, frloc::Number, fposloc::Number, srloc::Number, sposloc::Number)
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
function pointdistance(vars::VRP, sol::Solution, rloc::Number, rpos::Number)
    c = sol.routes[rloc].seq[rpos]
    if sol.routes[rloc].seqlen <= 1
        route_distance(sol.routes[rloc], vars)
    elseif rpos == 1
        vars.depot_distance[c] + vars.distance_m[c,sol.routes[rloc].seq[rpos+1]]
    elseif rpos == sol.routes[rloc].seqlen
        vars.distance_m[c,sol.routes[rloc].seq[rpos-1]] + vars.depot_distance[c]
    else
        vars.distance_m[c,sol.routes[rloc].seq[rpos-1]] + vars.distance_m[c,sol.routes[rloc].seq[rpos+1]]
    end
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
    oldp = pointdistance(vars,sol,frloc,fposloc) + pointdistance(vars,sol,srloc,sposloc)
    swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
    nd = pointdistance(vars,sol,frloc,fposloc) + pointdistance(vars,sol,srloc,sposloc)
    swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
    newobj = sol.objective - oldp + nd
    function commitChange()
        swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        sol.objective = newobj
    end
    return newobj, commitChange
end

function localSearch(sol::Solution, vars::VRP)
    b = 0
    w = 0
    t = 1e+6
    sol.objective = recalc_obj_val(sol, vars)
    T = 1000
    n = 1
    mx = vars.customers * (vars.customers - 1) / 2
    best_sol = deepcopy(sol)
    prob_func = (ns) -> exp((sol.objective - ns) / T)
    for _ = 1:t
        n += 1
        if n > mx
            n = 1
            T *= 0.95
        end
        nscore, changeFunc = randomNodeSwap(sol, vars)
        if nscore < best_sol.objective
            best_sol = sol
            changeFunc()
            b += 1
        elseif rand() <= prob_func(nscore)
            best_sol = deepcopy(sol)
            changeFunc()
            w += 1
        end
    end
    sol = best_sol
end

function nn_heur(vars :: VRP)
    groups = []
    cust = [1:vars.customers...]
    while length(cust) != 0
        ind = rand(1:length(cust))
        c = popat!(cust,ind)
        curr_group = [c]
        curr_load = vars.demand[c]
        cust_d = [(i,vars.distance_m[c,i]) for i in cust]
        sort!(cust_d,by = x->x[2])
        i = 1
        while i <= length(cust_d)
            if curr_load + vars.demand[cust_d[i][1]] > vars.capacity
                break
            end
            push!(curr_group,cust_d[i][1])
            filter!(x->x!=cust_d[i][1],cust)
            curr_load+=vars.demand[cust_d[i][1]]
            i+=1
        end
        push!(groups,Route(curr_group,length(curr_group),curr_load))
    end
    Solution(groups,0)
end
function nn_method(vars :: VRP)
    sol = nn_heur(vars)
    sol.objective = recalc_obj_val(sol,vars)
    sol2Opt(sol,vars)
    vis_output(sol,vars)
    sol
end
function mn(fl::String)
    start_time = Base.Libc.time()
    vars = read_input(fl)
    sol = getInitialSol(fl, vars)
    org_d = sol.objective
    sol2Opt(sol, vars)
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
    println("initial obj is ",org_d," final is ",f_d ," improv is ",round( ((org_d-f_d)/org_d)*100,digits=2),"%")
    vis_output(sol, vars)
    get_output(fl, end_time - start_time, sol, vars)
end

# function __init__()
#     mn(ARGS[1])
# end
# end