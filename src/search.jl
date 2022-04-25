module TRP
include("./parseInput.jl")

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
    sol.objective = recalc_obj_val(sol, vars)
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
    swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
    newdistance = recalc_obj_val(sol, vars)
    swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
    # @assert sol.routes[frloc].seq[fposloc] == first
    # @assert sol.routes[srloc].seq[sposloc] == second
    return newdistance, frloc, fposloc, srloc, sposloc
end

function localSearch(sol::Solution, vars::VRP)
    b = 0
    w = 0
    t = 1e+5
    sol.objective = recalc_obj_val(sol, vars)
    T = 100
    n = 1
    mx = vars.customers * (vars.customers - 1) / 2
    best_sol = deepcopy(sol)
    prob_func = (ns) -> exp((best_sol.objective - ns) / T)
    for _ = 1:t
        n += 1
        if n > mx
            n = 1
            T *= 0.95
        end
        res = checkSwapNodes(sol, vars)
        if res isa Bool
            continue
        else
            nscore = res[1]
            f1, f2, s1, s2 = res[2:end]
            if nscore < best_sol.objective
                swapNodes(sol, vars, f1, f2, s1, s2)
                sol.objective = nscore
                best_sol = sol
                b += 1
            elseif rand() <= prob_func(nscore)
                best_sol = deepcopy(sol)
                sol.objective = nscore
                swapNodes(sol, vars, f1, f2, s1, s2)
                w += 1
            end
        end
    end
    sol = best_sol
    # println("Made " * string(b) * " better swaps " * string(w) * " worse swaps out of " * string(t) * " swaps")
end
function mn(fl :: String)
    start_time = Base.Libc.time()
    vars = read_input(fl)
    sol = read_test_vrp(fl, vars)
    sol2Opt(sol, vars)
    # println("Distance before local search " * string(sol.objective))
    t = 3
    sols = Vector{Solution}(undef, t)
    for i = 1:t
        tmp = deepcopy(sol)
        localSearch(tmp, vars)
        sol2Opt(tmp, vars)
        tmp.objective = recalc_obj_val(tmp, vars)
        sols[i] = tmp
    end
    for i = 1:t
        # println("sols ", i, " objective is ", sols[i].objective)
        if sols[i].objective < sol.objective
            sol = sols[i]
        end
    end
    end_time = Base.Libc.time()
    # println("Distance after local search " * string(sol.objective))
    vis_output(sol, vars)
    get_output(fl,end_time-start_time,sol,vars)
end
function get_output(fl :: String,time,sol :: Solution,vars::VRP)
    kv_json = (k,v) -> join(['"',k,"""": """,'"',v,'"'])
    solo = []
    for route in sol.routes
        push!(solo,"0")
        for i=1:route.seqlen
            push!(solo,string(route.seq[i]))
        end
        push!(solo,"0")
    end
    out = join([
        "{",
        #TODO : Change to filename
        kv_json("Instance",fl),",",
        kv_json("Time",round(time,digits=2)),",",
        kv_json("Result",string(round(sol.objective,digits=2))) ,
        kv_json("Solution",join(solo," ")),
        "}"
    ])
    println(out)
end
function __init__()
    mn(ARGS[1])
end
end