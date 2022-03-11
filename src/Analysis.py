#!/usr/bin/env python
# coding: utf-8

# In[1]:


from scheduler import *
from main import *
import os


# In[2]:


os.environ["CP_SOLVER_EXEC"] = "/Applications/CPLEX_Studio129/cpoptimizer/bin/x86-64_osx/cpoptimizer"
set_context()


# In[3]:


def calc_waste(fl):
    cfg = Config.load(fl)
    sched = Scheduler(cfg)
    solution = sched.solve().schedule
    numWork = np.zeros(cfg.min_shifts.shape)
    numDaily = np.zeros((cfg.n_days))
    minDaily = np.full((cfg.n_days),cfg.min_daily)
    for worker in solution:
        for index in range(len(worker)):
            shift = worker[index]
            hours = shift[1] - shift[0]
            numDaily[index]+=hours
            if 0 <= shift[1] <= 8:
                numWork[index][1]+=1
                continue
            elif shift[1] <= 16:
                numWork[index][2]+=1
                continue
            elif  shift[1] <= 24:
                numWork[index][3]+=1
                continue
    sched_score =  (np.mean(numWork[:,1:]/cfg.min_shifts[:,1:])) -1  
    h_score =  np.mean(numDaily/minDaily) - 1
    return (sched_score * 100,h_score*100)


# In[4]:


files = ["../input/" + i for i in os.listdir("../input/")]
for file in files:
    res = calc_waste(file)
    print(f"Shift Waste for {file} is {res[0]} and hours waste is {res[1]}")


# In[ ]:




