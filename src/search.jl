# module TRP
using StatsBase
include("./parseInput.jl")
using JuMP
import HiGHS
import Ipopt

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
function swapNodes(vars::VRP, fr::Route, fposloc::Number, sr::Route, sposloc::Number)
    fdemand = vars.demand[fr.seq[fposloc]]
    sdemand = vars.demand[sr.seq[sposloc]]
    fr.load = fr.load - fdemand + sdemand
    sr.load = sr.load - sdemand + fdemand
    tmp = fr.seq[fposloc]
    fr.seq[fposloc] = sr.seq[sposloc]
    sr.seq[sposloc] = tmp
    complete2optSwap(fr, vars)
    complete2optSwap(sr, vars)
    @assert calc_route_load(fr, vars) == fr.load
    @assert calc_route_load(sr, vars) == sr.load

end

function calc_route_load(route::Route, vars::VRP)
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
function randomNodeSwap(sol::Solution, vars::VRP)
    route_load_change = (rn, f, s) -> sol.routes[rn].load - f + s
    frloc, srloc = rand(1:vars.vehicles, 2)
    if frloc == srloc
        return true, Inf, () -> nothing
    end
    # TODO : ROUTES could be empty !
    fposloc = rand(1:sol.routes[frloc].seqlen)
    sposloc = rand(1:sol.routes[srloc].seqlen)
    first = sol.routes[frloc].seq[fposloc]
    second = sol.routes[srloc].seq[sposloc]
    fdemand = vars.demand[first]
    sdemand = vars.demand[second]
    PF = 0
    if route_load_change(frloc, fdemand, sdemand) > vars.capacity || route_load_change(srloc, sdemand, fdemand) > vars.capacity
        return true, Inf, () -> nothing
    end
    oldp = route_distance(sol.routes[frloc], vars) + route_distance(sol.routes[srloc], vars)
    frcp, srcp = deepcopy(sol.routes[frloc]), deepcopy(sol.routes[srloc])
    swapNodes(vars, frcp, fposloc, srcp, sposloc)
    # @assert frcp.load == route_load_change(frloc,fdemand,sdemand)
    # @assert srcp.load ==  route_load_change(srloc, sdemand, fdemand)
    nd = route_distance(frcp, vars) + route_distance(srcp, vars) + PF

    newobj = sol.objective - oldp + nd
    function commitChange()
        sol.routes[frloc] = frcp
        sol.routes[srloc] = srcp
        sol.objective = newobj
    end
    return (PF == 0), newobj, commitChange
end
function checkSol(sol::Solution, vars::VRP)
    for route in sol.routes
        if route.load > vars.capacity
            return false
        end
    end
    return true
end

function localSearch(sol::Solution, vars::VRP)

    b = 0
    w = 0
    t = 1e+4 * vars.customers

    sol.objective = recalc_obj_val(sol, vars)
    init = sol.objective
    temperature = 1000
    num_runs = 1
    mx_runs = vars.customers * (vars.customers - 1) / 2
    best_sol = deepcopy(sol)
    prob_func = (ns) -> exp((sol.objective - ns) / temperature)
    checkSol(best_sol, vars)
    println("running local search")
    for i = 1:numItr
        checkSol(best_sol, vars)
        num_runs += 1
        if num_runs > mx_runs
            num_runs = 1
            temperature *= 0.95
        end
        incapacity, nscore, changeFunc = randomNodeSwap(sol, vars)
        if nscore < best_sol.objective
            changeFunc()

            best_sol = deepcopy(sol)
            b += 1

        elseif rand() <= prob_func(nscore)
            changeFunc()
        end
    end

    sol.objective = recalc_obj_val(sol,vars)
    newsol = deepcopy(best_sol)
    newsol
end

# function localSearch(sol :: Solution, vars :: VRP)
#     t = 1e+4 * vars.customers
#     sol.objective = recalc_obj_val(sol, vars)
#     for _ = 1:t
#         nscore, changeFunc = randomNodeSwap(sol, vars)
#         if nscore < sol.objective
#             changeFunc()
#         end
#     end
#     nothing
# end
function nn_heur(vars :: VRP)

    groups = []
    cust = [1:vars.customers...]
    while length(cust) != 0
        ind = rand(1:length(cust))
        c = popat!(cust, ind)
        curr_group = [c]
        curr_load = vars.demand[c]
        cust_d = [(i, vars.distance_m[c, i]) for i in cust]
        sort!(cust_d, by=x -> x[2])
        for i = 1:length(cust_d)
            if curr_load + vars.demand[cust_d[i][1]] > vars.capacity
                break
            else
                push!(curr_group, cust_d[i][1])
                filter!(x -> x != cust_d[i][1], cust)
                curr_load += vars.demand[cust_d[i][1]]
            end
        end
        push!(groups, Route(curr_group, length(curr_group), curr_load))
    end
    Solution(groups, 0)
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
function nn_method(vars::VRP, l::Bool)
    sol = nn_heur(vars)
    sol.objective = recalc_obj_val(sol, vars)
    sol2Opt(sol, vars)
    if l
        localSearch(sol, vars)
        sol.objective = recalc_obj_val(sol, vars)
        sol2Opt(sol, vars)
    else
        sol.objective = recalc_obj_val(sol, vars)
    end
    sol
end

function get_customer_routes(sol :: Solution)
    c_r = Dict{Number,Number}()
    for i = 1:length(sol.routes)
        for j=1:sol.routes[i].seqlen
            c_r[sol.routes[i].seq[j]] =  i 
        end
    end
    c_r
end
function destroy_tsp(sol:: Solution, vrp :: VRP)
    remove_c = sample(1:vrp.customers,Int(round(0.75*vrp.customers)),replace=false)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    C = 1:vrp.customers
    V = 1:vrp.vehicles
    mtrx = @variable(model, x[V, C], binary = true)
    @constraint(model, [c in C], sum(x[:, c]) == 1)
    @constraint(model, [v in V], sum(x[v, :] .* vrp.demand) <= vrp.capacity)
    cust_routes = get_customer_routes(sol)
    remove_mat = zeros(vrp.vehicles,vrp.customers)
    min_expr = @expression(model,0)
    for c in C
        curr_route = cust_routes[c]
        if !(c in remove_c)
            
            @constraint(model,  x[curr_route,c] == 1)
        else
            min_expr += @expression(model,x[curr_route,c])
        end
    end
    # bad_cust_r = cust_routes[remove_c[1]]
    # @constraint(model,x[bad_cust_r,remove_c[1]] == 0)
    @objective(model,Min,min_expr)
    optimize!(model)
    lpvals = zeros(vrp.vehicles, vrp.customers)
    for i = 1:vrp.vehicles
        for j = 1:vrp.customers
            lpvals[i, j] = value(mtrx[i, j])
        end
    end
    println("available sols : ",result_count(model))
    nsol = solverToSol(vrp,lpvals)
    sol2Opt(nsol,vrp)
    lnsol = localSearch(nsol,vrp)
    sol2Opt(lnsol,vrp)
    lnsol
end
function initStuff(fl :: String) 
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
    vars,sol
end
    
function mn(fl::String)
    start_time = Base.Libc.time()
    vars = read_input(fl)
    sol = getInitialSol(vars)
    sol2Opt(sol, vars)
    org_d = sol.objective
    numRuns = 2
    for i = 1:numRuns
        tmp = nn_method(vars, false)
        if tmp.objective < sol.objective
            sol = tmp
        end
    end
    end_time = Base.Libc.time()
    sol2Opt(sol, vars)
    vis_output(sol, vars)
    get_output(fl, end_time - start_time, sol, vars)
end
# function __init__()
#     nn_mn(ARGS[1])
# end
# end