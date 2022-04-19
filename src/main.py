from dataclasses import dataclass
from math import dist
import os
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path
import numpy as np
from docplex.cp.config import context
from docplex.cp.config import *
from docplex.cp.model import *


def set_context():
    if not("CP_SOLVER_EXEC" in os.environ):
        solver_exec = Path("src/bin/cpoptimizer")
    else:
        solver_exec = Path(os.environ["CP_SOLVER_EXEC"])
        if not solver_exec.exists():
            solver_exec = Path("src/bin/cpoptimizer")
    # print(f"Using cp installation at {solver_exec}")
    context.solver.agent = "local"
    context.solver.local.execfile = str(solver_exec)


parser = ArgumentParser()
parser.add_argument("input_file", type=str)

@dataclass()
class VRP:
    customers : int
    vehicles : int
    capacity : int
    demand  : list[int]
    depot_location : tuple
    positions : np.ndarray
    distance_matrix : np.ndarray
@dataclass()
class Solution:
    objective : float
    opt : int
    routes : list

def get_distance(pos1,pos2):
    return np.sqrt( ((pos1[0] -pos2[0])**2) + ((pos1[1] -pos2[1])**2))
def create_distance_matrix(positions):
    dists = np.zeros((len(positions),len(positions)))
    for i in builtin_range(positions.shape[0]):
        for j in builtin_range(positions.shape[1]):
            if i!=j:
                dists[i,j] = get_distance(positions[i],positions[j])
    return dists
def read_input(fl : Path):
    with open(fl,"r") as f:
        lines = f.readlines()
    lines = [i.split() for i in lines]
    customers,vehicles,capacity = [int(float(i)) for i in lines[0]]
    depot_location = (int(float(lines[1][1])),int(float(lines[1][2])))
    demand = []
    positions = []
    for i in builtin_range(2,customers+1):
        demand.append(int(float(lines[i][0])))
        positions.append([int(float(lines[i][j])) for j in builtin_range(1,3)])
    positions = np.array([np.array(i) for i in positions])
    
    return VRP(customers-1,vehicles,capacity,demand,depot_location,positions,create_distance_matrix(positions))
def main(args):
    set_context()
    input_file = Path(args.input_file)
    vars =  read_input(input_file)
    solve(vars)

def solve(vars : VRP):
    model = CpoModel()
    customer_assigs = []
    for c in builtin_range(0,vars.customers):
        customer_assigs.append(
          [integer_var(1,vars.vehicles),integer_var(1,vars.customers)]
        )
    customer_assigs = np.array([np.array(i) for i in customer_assigs])
    #Ordering is unique
    for i in builtin_range(0,vars.customers):
        for j in builtin_range(0,vars.customers):
            if i!=j:
                model.add(if_then(customer_assigs[i,0] == customer_assigs[j,0],
            customer_assigs[i,1] != customer_assigs[j,1]))
    truck_capacities  =[[] for i in builtin_range(vars.vehicles)] # List of Lists 

    for c in builtin_range(0,len(customer_assigs)):
        for t in builtin_range(1,vars.vehicles+1):
           truck_capacities[t-1].append(times(count([customer_assigs[c,0]],t),vars.demand[c]))
    for vehicle_demand in truck_capacities:
        model.add(sum(vehicle_demand) <= vars.capacity)
    #Got to be consecutive
    # customers_per_truck = [0 for i in range(vars.vehicles)]
    for t in builtin_range(1,vars.vehicles+1):
        customers_per_truck = count(customer_assigs[:,0].ravel().tolist(),t)
        for c in builtin_range(vars.customers):
            model.add(if_then(
                customer_assigs[c,0] == t,
                customer_assigs[c,1] <= customers_per_truck
            ))

     # Objective   
    # Transition matrix
    trans_mat = []
    for i in range(vars.customers):
        trans_mat.append(binary_var_list(vars.customers))
    trans_mat = np.array([np.array(i) for i in trans_mat])
    for i in range(vars.customers):
        for j in range(vars.customers):
            if i==j:
                model.add(trans_mat[i,j] == 0)
            else:
                
                model.add(
                    if_then(
                        logical_and((customer_assigs[i,0] == customer_assigs[j,0]) , (plus(customer_assigs[i,1],1) == customer_assigs[j,1])),
                        trans_mat[i,j] == 1

                    )
                 )
                model.add(
                    if_then(
                       logical_not(logical_and(equal(customer_assigs[i,0],customer_assigs[j,0]) ,  equal(plus(customer_assigs[i,1],1) , customer_assigs[j,1]))),
                        trans_mat[i,j] == 0

                    )
                 )
    distExpr = []
    for i in builtin_range(vars.customers):
        for j in builtin_range(vars.customers):
            if i!=j:
                distExpr.append(times(trans_mat[i,j],vars.distance_matrix[i,j]))
    model.minimize(sum(distExpr))
    model.set_parameters(CpoParameters(SolutionLimit=3))
    solution : CpoSolveResult = model.solve()
    if not solution.is_solution():
        return 
    sol_routes = np.zeros((vars.vehicles,vars.customers))
    
    for c in range(vars.customers):
        t,n = solution[customer_assigs[c,0]],solution[customer_assigs[c,1]]
        sol_routes[t-1,n] = c + 1
    
    
    s  = get_stdout(vars,sol_routes.astype(int),solution.is_solution_optimal())
    write_sol(s,"test.vrp")
def write_sol(s : Solution,fl):
    with open(fl,"w") as f:
        f.write(f"{s.objective} {0} \n")
        for route in s.routes:
            f.write(f'{" ".join([str(i) for i in route])} \n')
    f.close()
        

def get_stdout(vars :VRP,sol_routes,opt):
    opt = 1 if opt else 0
    obj = 0
    all_routes = []
    for i in builtin_range(sol_routes.shape[0]):
        # print(f"i is {i}")/
        route = []
        for j in builtin_range(sol_routes.shape[1]):
            # print(f"j is {j}")
            route.append(sol_routes[i,j])
            if j==0:
                obj+=get_distance(vars.depot_location,vars.positions[sol_routes[i,j+1]-1])
            elif sol_routes[i,j+1] == 0:
                obj+=get_distance(vars.positions[sol_routes[i,j]-1],vars.depot_location)
                route.append(0)
                break
            else:
                obj+=vars.distance_matrix[sol_routes[i,j]-1,sol_routes[i,j+1]-1]
        all_routes.append(route)
    return Solution(obj,opt,all_routes)

    
if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
