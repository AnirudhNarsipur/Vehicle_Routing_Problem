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
    # demand, customer node number 
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