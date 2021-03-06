module TRP
include("./parseInput.jl")
"""
Swaps Two customers in the Solution
"""
function swapNodes(sol::Solution, vars::VRP, frloc::Int64, fposloc::Int64, srloc::Int64, sposloc::Int64)
    fdemand = vars.demand[sol.routes[frloc].seq[fposloc]]
    sdemand = vars.demand[sol.routes[srloc].seq[sposloc]]
    sol.routes[frloc].load = sol.routes[frloc].load - fdemand + sdemand
    sol.routes[srloc].load = sol.routes[srloc].load - sdemand + fdemand
    tmp = sol.routes[frloc].seq[fposloc]
    sol.routes[frloc].seq[fposloc] = sol.routes[srloc].seq[sposloc]
    sol.routes[srloc].seq[sposloc] = tmp
    nothing
end
"""
Given a node : Y in the route A->X->Y->Z->B calculates distance X->Y->Z
"""
function pointdistance(vars::VRP, sol::Solution, rloc::Int64, rpos::Int64)::Float64
    c::Int64 = sol.routes[rloc].seq[rpos]
    if sol.routes[rloc].seqlen <= 1
        vars.depot_distance[c]
    elseif rpos == 1
        vars.depot_distance[c] + vars.distance_m[c, sol.routes[rloc].seq[rpos+1]]
    elseif rpos == sol.routes[rloc].seqlen
        vars.distance_m[c, sol.routes[rloc].seq[rpos-1]] + vars.depot_distance[c]
    else
        vars.distance_m[c, sol.routes[rloc].seq[rpos-1]] + vars.distance_m[c, sol.routes[rloc].seq[rpos+1]]
    end
end
"""
Core Local Search (Simulated Annealing)
"""
function localSearch(sol::Solution, vars::VRP)::Solution
    numItr::UInt64 = round((154579.1256 + 150.2166590 * vars.customers + 54.10538 * vars.customers^2))
    Temperature::Float64 = abs(sol.objective / log(MathConstants.e, 0.97)) * 0.001
    numRuns = 1
    annealing_cycle::Int64 = (vars.customers * (vars.customers - 1) / 2)
    prob_func = (ns) -> exp((sol.objective - ns) / Temperature)
    best_sol = deepcopy(sol)
    for _ = 1:numItr
        numRuns += 1
        if numRuns > annealing_cycle
            numRuns = 1
            Temperature *= 0.95
        end
        frloc::Int64, srloc::Int64 = sample(1:vars.vehicles, 2, replace=false)
        if sol[frloc].seqlen == 0 || sol[srloc].seqlen == 0
            continue
        end
        fposloc::Int64 = rand(1:sol[frloc].seqlen)
        sposloc::Int64 = rand(1:sol[srloc].seqlen)
        first::Int64 = sol[frloc][fposloc]
        second::Int64 = sol[srloc][sposloc]
        if (sol[frloc].load - vars.demand[first] + vars.demand[second]) > vars.capacity || (sol[srloc].load - vars.demand[second] + vars.demand[first]) > vars.capacity
            continue
        end
        oldp::Float64 = pointdistance(vars, sol, frloc, fposloc) + pointdistance(vars, sol, srloc, sposloc)
        swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        nd::Float64 = pointdistance(vars, sol, frloc, fposloc) + pointdistance(vars, sol, srloc, sposloc)
        nscore::Float64 = sol.objective - oldp + nd
        if nscore < best_sol.objective
            sol.objective = nscore
            best_sol = deepcopy(sol)
            continue
        elseif frloc != srloc && rand() <= prob_func(nscore)
            sol.objective = nscore
            continue
        else
            swapNodes(sol, vars, frloc, fposloc, srloc, sposloc)
        end
    end
    sol2Opt(best_sol, vars)
    best_sol
end
"""
Driver code
"""
function mn(fl::String)
    start_time = Base.Libc.time()
    vars = read_input(fl)
    sol = getInitialSol(vars)
    numRuns = 20
    bestsol = sol
    for _ = 1:numRuns
        tmp = deepcopy(sol)
        tmp = localSearch(tmp, vars)
        if tmp.objective < bestsol.objective
            bestsol = deepcopy(tmp)
        end
    end
    end_time = Base.Libc.time()
    get_output(fl, end_time - start_time, bestsol)
end
function __init__()
    if length(ARGS) == 0
        nothing
    else
        mn(ARGS[1])
    end
end
end


