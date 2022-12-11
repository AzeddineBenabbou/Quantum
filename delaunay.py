#basic algorithm to compute the delaunay triangulation of a set points in 
# 1D and 2D 
# I use CQM
from dimod import Binary, ConstrainedQuadraticModel, quicksum
from dwave.system import LeapHybridCQMSampler
import numpy as np

#define geometry of points
n=5
geom=np.array([9, -1, 1, 0, 2])

#define the edges (variables of current problem)
# edges ={{0, 1}, {1, 2}, {2, 3} ...{n-2, n-1} {0, 2}, {1, 3} ... {n-3, n-1} .....{0, n-2}, {1, n-1}, {0, n-1}}
# e1=[1 2]
nb_var=int((n*n-n)/2)
edges=np.array([0, 1, 1, 2, 2, 3, 3, 4, 0, 2, 1, 3, 2, 4, 0, 3, 1, 4, 0, 4 ])

#define lengths
length=[abs(geom[edges[2*i]]-geom[edges[2*i+1]]) for i in range(nb_var)]

print(length)

#define lifting map y=x^2
y=[pow(geom[i], 2) for i in range(geom.size)]

print(y)

#define weights
w=[length[i]*(y[edges[2*i]]+y[edges[2*i+1]]) for i in range(nb_var)]
print(w)

# define binary variables
bin_variables = [Binary(f'e_{i}') for i in range(nb_var)]

# print what kind of variables we are using
#print(bin_variables[0])

# instantiate a cqm
cqm=ConstrainedQuadraticModel()

# define the objective function
objective1=quicksum(w[i]*bin_variables[i] for i in range(nb_var))

print(objective1)

objective2=quicksum(-1*bin_variables[i] for i in range(nb_var))

print(objective2)


# set objectives
#cqm.set_objective(objective1)
cqm.set_objective(objective2)

# #define our constraint
# the number of edges = npts-1
cqm.add_constraint(quicksum(bin_variables[i] for i in range(nb_var))==n-1)
# sum(lenghts)=(9+1)=10
cqm.add_constraint(quicksum(length[i]*bin_variables[i] for i in range(nb_var))==10)


# #print the number of variables 
# print("number of variables for this problem {}".format(len(cqm.variables)))

# submit to a solver
#CQM sampler
# instantiate a cqm solver (sampler)
cqm_sampler = LeapHybridCQMSampler()

# #solve with the instance of the cqm solver
sampleset = cqm_sampler.sample_cqm(cqm, time_limit=5, label ='Delaunay 1d Demo with max nb_edges')
print("all samplest aggregated")
print(sampleset.aggregate())

# filter results. Only feasible solutions.
print("feasible solutions aggregated")
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
fsamples_agr=feasible_sampleset.aggregate()
print(fsamples_agr)




  
