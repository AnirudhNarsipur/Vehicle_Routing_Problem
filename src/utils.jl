function recalc_obj_val(sol::Solution, vars::VRP) :: Float64
    sum(map(r -> route_distance(r, vars), sol.routes))
end

function route_distance(route::Route, vars::VRP)
    if route.seqlen==0
        return 0
    end
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
"""
Return a dictionary of customer route assigments
"""
function get_customer_routes(sol :: Solution) :: Dict{Number,Number}
    c_r = Dict{Number,Number}()
    for i = 1:length(sol.routes)
        for j=1:sol.routes[i].seqlen
            c_r[sol.routes[i].seq[j]] =  i 
        end
    end
    c_r
end

"""
1. take route[1] to route[i-1] and add them in order to new_route
2. take route[i] to route[k] and add them in reverse order to new_route
3. take route[k+1] to end and add them in order to new_route
"""
function opt2Swap(route::Route, i::Number, k::Number)
    reverse!(view(route.seq, i:k))
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
                if (route_distance(route, vars) < best_dist)
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
    sol.nodeloc = customer_route_loc(sol)

end

"""
Calculate the total load of the current route
"""
function calc_route_load(route::Route, vars::VRP)::Number
    reduce((base,elem) -> base + vars.demand[route.seq[elem]],1:route.seqlen,init=0)
end

"""
Prints current route lengths
"""
function routeLengths(sol::Solution)
    print("Route length are ")
    for route in sol.routes
        print(route.seqlen, " ")
    end
    println("")
end

#### Nearest Neighbor Based

  


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
    Solution(groups, 0,Dict())
end
function nn_method(vars::VRP)
    sol = nn_heur(vars)
    sol.objective = recalc_obj_val(sol, vars)
    sol2Opt(sol, vars)
    vis_output(sol, vars)
    sol
end

function customer_route_loc(sol :: Solution)
    dct = Dict{Number,Tuple{Number,Number}}()
    for (i,route) in enumerate(sol.routes)
        for j in 1:route.seqlen
            dct[route.seq[j]] = (i,j)
        end
    end
    dct
end

function solutionCheck(sol :: Solution,vars :: VRP)
    #Loads are correct 
    for route in sol.routes
        @assert(route.load == calc_route_load(route,vars))
    end
    for c in 1:vars.customers
        @assert(findCloc(sol,c) == sol.nodeloc[c])
    end
end