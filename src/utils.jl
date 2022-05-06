function recalc_obj_val(sol::Solution, vars::VRP)
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