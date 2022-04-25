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