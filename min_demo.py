#basic algorithm to compute the minimum from a list
#in this demo I will use CQM
from dimod import Binary, ConstrainedQuadraticModel, quicksum
from dwave.system import LeapHybridCQMSampler
import numpy as np

#give number of entries in your problem
num_i=5

#define an array of size num_i
w=np.array([9, -1, -2, 0, -2])

#define binary variables
bin_variables = [Binary(f'x_{i}') for i in range(num_i)]

#print what kinf of variables we are using
print(bin_variables[0])

#instantiate a cqm
cqm=ConstrainedQuadraticModel()

#define the objective
objective=quicksum(w[i]*bin_variables[i] for i in range(num_i))

print(objective)

#set the objective
cqm.set_objective(objective)

#define our constraint
#pick just one and only one element
cqm.add_constraint(quicksum(bin_variables[i] for i in range(num_i))==1)


#print the number of variables 
print("number of variables for this problem {}".format(len(cqm.variables)))

#submit to a solver
# CQM sampler
#instantiate a cqm solver (sampler)
cqm_sampler = LeapHybridCQMSampler()

#solve with the instance of the cqm solver
sampleset = cqm_sampler.sample_cqm(cqm, time_limit=5, label ='Minimum Demo')
print("all samplest aggregated")
print(sampleset.aggregate())

# filter results. Only feasible solutions.
print("feasible solutions aggregated")
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
fsamples_agr=feasible_sampleset.aggregate()
print(fsamples_agr)


# for sample in fsamples_agr.samples():   
#     print(sample)


#print(feasible_sampleset)

#

# #check if we have a feasible solution
# #print number of feasible solutions among all soltions
# if len(feasible_sampleset):
#     best=feasible_sampleset.first
#     print(best.sample.items())
#     print("{} feasible solutions of {}.".format(len(feasible_sampleset), len(sampleset)))



# #print the best solution
# selected_sols=[key for key, val in feasible_sampleset.Minimum if "x_" in key and val]
# print("{} minima found".format(len(selected_sols)))
# print("selected_sols")
# print(selected_sols)




# #get the indices of the best soltion(s)
# #def get_indices(name):
#     #return [int(digs) for digs in name.split('_') if digs.isdigit()]

# #

# #print all best solution(s)
# for sol in selected_sols:
#     in_sol=[key for key, val in best.sample.items() if "x_" in key]
#     result=[get_indices(item)[0] for item in in_sol]
#     print("sol with index {}  has value {}.".format(result, w[result]))





#print(feasible_sampleset)

#results of my last run
# values are:
#[ 9 -1  3  0 -2]
#{0.0, 0.0, 0.0, ,
#  'vf73de6e89f6b418fa1ec71f92d466870': 1.0}

  