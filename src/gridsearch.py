import sched
import numpy as np
from docplex.cp.config import *
from docplex.cp.model import *
from scheduler import *
import itertools
var_evaluators = {
"domain_size":domain_size,
"domain_max":domain_max,
"domain_min":domain_min,
"var_impact":var_impact,
"var_local_impact":var_local_impact,
"var_success_rate":var_success_rate,
"impact_of_last_branch":impact_of_last_branch,
}
val_selectors = [select_largest,select_smallest,select_random_value] 
var_selectors = [select_largest,select_smallest,select_random_var]
val_evaluators = {
  "value" : value,
  "value_impact":value_impact,
  "value_success_rate" : value_success_rate
  }
need_param = ["value_index","explicit_value_eval","var_index","explicit_var_eval"]
files = ["input/7_17.sched","input/14_14.sched","input/28_40.sched"]
def hours_vars(model : Scheduler):
    return model.hours.ravel().tolist()
def shift_vars(model : Scheduler):
    return model.shifts.ravel().tolist()

def no_vars(model : Scheduler):
    return None
def training_phase(model : Scheduler):
    return model.hours[:, : model.config.n_shifts].ravel().tolist()  + model.shifts[:, : model.config.n_shifts].ravel().tolist()
    

def run_experiment(file : str):
    cfg = Config.load(file)
    timelimit = 1
    search_space = itertools.product(var_evaluators.keys(),val_evaluators.keys(),var_selectors,val_selectors)
    for search in search_space:
        print("Starting search ",search)
        var_k , val_k ,var_sel,val_sel  = search      
        sched = Scheduler(cfg)
        params = CpoParameters(SearchType="DepthFirst",TimeLimit=timelimit)
        sched.model.set_parameters(params)
        sched.model.set_search_phases([search_phase(varchooser=var_sel(var_evaluators[var_k]()),valuechooser=val_sel(val_evaluators[val_k]()))])
        sched.model.solve()
        print("For ",file," and params ",var_k," : ",val_k," in " ,timelimit," search status is ",sched.model.get_solve_status()," took ",sched.model.get_solve_time(), " seconds")
run_experiment(files[0])
        




