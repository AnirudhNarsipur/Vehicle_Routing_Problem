from dataclasses import dataclass
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

def read_input(fl : str):
    with open(fl,"r") as f:
        lines = f.readlines()
    lines = [i.split() for i in lines]
    customers,vehicles,capacity = [int(i) for i in lines[0]]
    depot_location = (int(lines[1][1]),int(lines[1][2]))
    demand = []
    positions = []
    for i in range(2,customers+1):
        demand.append(int(lines[i][0]))
        positions.append([int(lines[i][j]) for j in range(1,3)])
    positions = np.array([np.array(i) for i in positions])
    return VRP(customers-1,vehicles,capacity,demand,depot_location,positions)
def main(args):
    set_context()
    input_file = Path(args.input_file)
    filename = input_file.name
    vars =  read_input(filename)
    solve(vars)

def solve(vars):
    model = CpoModel()
    truck_routes = []
    for _ in range(0,vars.vehicles):
        truck_routes.append(
            integer_var_list(vars.customers+2,-1,vars.customers)
        )
    truck_routes = np.array(truck_routes)
    #always start at depot
    for i in range(0,len(truck_routes)):
        model.add(equals(
            truck_routes[i,0],0
        ))
    # each customer occurs exactly once
    for c in range(1,vars.customers+1):
        model.add(equal(count(truck_routes.ravel().tolist(),c), 1))
    # total served by each truck <= capacity
    for i in range(0,len(truck_routes)):
        acc = 0
        for j in range(0,truck_routes.shape[1]):
            model.add(if_then(truck_routes[i,j] > 0,sum(acc,vars.demand[truck_routes[i,j]])))
        model.add(acc <= vars.capacity)



if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
