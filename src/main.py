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

def read_input(fl : Path):
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
    vars =  read_input(input_file)
    solve(vars)

def solve(vars : VRP):
    model = CpoModel()
    customer_assigs = []
    for c in range(0,vars.customers):
        customer_assigs.append(
          [integer_var(1,vars.vehicles),integer_var(1,vars.customers)]
        )
    customer_assigs = np.array([np.array(i) for i in customer_assigs])
    #Ordering is unique
    for i in range(0,len(customer_assigs)):
        for j in range(0,len(customer_assigs)):
            if i!=j:
                model.add(if_then(customer_assigs[i,0] == customer_assigs[j,0],
            customer_assigs[i,1] != customer_assigs[j,1]))
    truck_capacities  =[[] for i in range(vars.vehicles)] # List of Lists 

    for c in range(0,len(customer_assigs)):
        for t in range(1,vars.vehicles+1):
           truck_capacities[t-1].append(times(count([customer_assigs[c,0]],t),vars.demand[c]))
    for vehicle_demand in truck_capacities:
        model.add(sum(vehicle_demand) <= vars.capacity)
    model.solve()



if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
