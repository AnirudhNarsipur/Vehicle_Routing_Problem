using JuMP
import HiGHS
using StatsBase
using Random
struct VRP
    customers::Int64
    vehicles::Int64
    capacity::Int64
    demand::Vector{Int64}
    depot_location::Vector{Float64}
    positions::Vector{Vector{Float64}}
    distance_m::Matrix{Float64}
    depot_distance::Vector{Float64}
    # Node -> Position in node_pos
    node_pos::Dict{Int64,Int64}
    #  customer node number  , demand
    node_demand::Vector{Tuple{Int64,Int64}}
    sorted_d :: Matrix
end
mutable struct Route
    seq::Vector{Int64}
    seqlen::Int64
    load::Int64
end

mutable struct Solution
    routes::Vector{Route}
    objective::Float64
    nodeloc :: Dict{Int64,Tuple{Int64,Int64}}
end
Base.getindex(s :: Solution,i :: Int64) =  s.routes[i]
Base.getindex(r :: Route,i :: Int64) = if i <= r.seqlen r.seq[i] else missing end