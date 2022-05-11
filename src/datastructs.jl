using JuMP
import HiGHS
struct VRP
    customers::Number
    vehicles::Number
    capacity::Number
    demand::Vector
    depot_location::Vector
    positions::Vector
    distance_m::Matrix
    depot_distance::Vector
    # Node -> Position in node_pos
    node_pos::Dict{Number,Number}
    #  customer node number  , demand
    node_demand::Vector{Tuple{Number,Number}}
    sorted_d :: Matrix
end
mutable struct Route
    seq::Vector{Number}
    seqlen::Number
    load::Number
end

mutable struct Solution
    routes::Vector{Route}
    objective::Number
    nodeloc :: Dict{Number,Tuple{Number,Number}}
end
Base.getindex(s :: Solution,i :: Int64) =  s.routes[i]
Base.getindex(r :: Route,i :: Int64) = if i <= r.seqlen r.seq[i] else missing end