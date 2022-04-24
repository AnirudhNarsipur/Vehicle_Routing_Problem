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
    '''
    Sets up CPLEX 
    '''
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
    '''
    Represents VRP internally
    '''
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
    solve(vars,input_file)

def solve(vars : VRP,inp):
    model = CpoModel()
    truck_assigs = []
    for i in builtin_range(vars.vehicles):
        truck_assigs.append([binary_var() for i in builtin_range(vars.customers)])
    truck_assigs = np.array([np.array(i) for i in truck_assigs])

    for i in builtin_range(vars.vehicles):
        model.add(scal_prod(truck_assigs[i].tolist() ,vars.demand) <= vars.capacity)
    for i in builtin_range(vars.customers):
        model.add(sum(truck_assigs[:,i].tolist()) == 1)
    model.set_parameters(CpoParameters(LogVerbosity="Terse"))
    solution : CpoSolveResult = model.solve()
    if not solution.is_solution():
        return 
    sol_routes = []
    for t in builtin_range(vars.vehicles):
        tmp =[]
        for c in builtin_range(vars.customers):
            if solution[truck_assigs[t,c]] == 1:
                tmp.append(c)
        sol_routes.append(tmp)
    sol = np.zeros((vars.vehicles,vars.customers))
    for t in builtin_range(vars.vehicles):
        for c in builtin_range(vars.customers):
            if c < len(sol_routes[t]):
                sol[t,c] = sol_routes[t][c] + 1
    s  = get_stdout(vars,sol.astype(int),solution.is_solution_optimal())
    print_sol(s,solution.get_solve_time(),inp)
    write_sol(s,"test.vrp")
def print_sol(s : Solution,time,inp :Path):
    dct = {"Instance" :os.path.basename(inp),"Time":time,"Result":s.objective,"Solution": " ".join([" ".join([str(j) for j in i]) for i in s.routes]) }
    print(json.dumps(dct))
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
        route = [0]
        for j in builtin_range(sol_routes.shape[1]):
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
