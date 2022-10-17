#basic algorithm to compute the maximu from a list
#in this demo I will use CQM
from dimod import dimod, Binaries, ConstrainedQuadraticModel, quicksum
from dwave.system import LeapHybridCQMSampler
import numpy as np

#give number of entries in your problem
num_i=5

#define an array of size num_i
w=np.array([9, -1, 3, 0, -2])

#define binary variables
bin_variables = list(Binaries(num_i))

#instantiate a cqm
cqm=ConstrainedQuadraticModel()

#define the objective
objective=-quicksum(w[i]*bin_variables[i] for i in range(num_i))

#set the objective
cqm.set_objective(objective)

#define our constraint
#pick just one and only one element
cqm.add_constraint(quicksum(bin_variables[i] for i in range(num_i))==1)

#submit to a solver
# CQM sampler
#instantiate a cqm solver (sampler)
cqm_sampler = LeapHybridCQMSampler()

#solve with the instance of the cqm solver
sampleset = cqm_sampler.sample_cqm(cqm, time_limit=5, label ='Minimum Demo')


#filter results. Only feasible solutions.
feasible_sampleset = sampleset.filter(lambda d: d.is_feasible)

first_feasible_sol = feasible_sampleset.first.sample

#print the resulting
print("values are:" )
print(w)
print(first_feasible_sol)
print(feasible_sampleset)
print(sampleset)

#results of my last run
# values are:
#[ 9 -1  3  0 -2]
#{0.0, 0.0, 0.0, , 'vf73de6e89f6b418fa1ec71f92d466870': 1.0}

  