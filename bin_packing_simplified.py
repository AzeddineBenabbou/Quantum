from dimod import ConstrainedQuadraticModel, Binaries, ExactCQMSolver, quicksum, Binary
from dwave.system import LeapHybridCQMSampler
import numpy as np

#example ffrom dwave documentation: https://docs.ocean.dwavesys.com/en/stable/docs_dimod/reference/models.html?highlight=cqm#module-dimod.constrained.constrained
#the obective if minimize the number of used bins while not exceding the capacity of each bin
#weights of packs
weights = [.9, .7, .2, .1]
#capacity of each bin
capacity = 2
#variables are for bins: if yj=1 means that bin j is used (note that we suppose number of bins = number of packs)
y = [Binary(f'y_{j}') for j in range(len(weights))]
#xij are variables that indicate if a pack i is put in bin j
x = [[Binary(f'x_{i}_{j}') for j in range(len(weights))]
     for i in range(len(weights))]

#create empty cqm
cqm=ConstrainedQuadraticModel()

#objective
cqm.set_objective(sum(y))

#constraints

#constraint1: each pack must be put at exactly one bin
for i in range(len(weights)):
    cqm.add_constraint(sum(x[i]) == 1, label=f'item_placing_{i}')

#constraint2: weights of packs put in a bin cannot exceed bin's capacity
for j in range(len(weights)):
    cqm.add_constraint(
        sum(weights[i] * x[i][j] for i in range(len(weights))) - y[j] * capacity <= 0,
        label=f'capacity_bin_{j}')

#variant with additional constraint from me
#constraint3: packs with big different weights cannot be in the same bin (crieteria w[i]-w[k]<0.3)
#cqm.add_constraint(quicksum(x[i][j]*x[k][j] for j in range(len(weights)) for i in range(len(weights)-1) for k in range(i+1, len(weights)) if (abs(weights[i]-weights[k])>=0.3))==0, label=f'packs with small wieghts diff in the same bin')

#solve
sampleset = ExactCQMSolver().sample_cqm(cqm)
all_samples=sampleset.samples
#print(all_samples)
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
fsamples_agr=feasible_sampleset.aggregate()
print('---------------------------------------------FEASIBLE-------------------------------------------------------------')
s=feasible_sampleset.first
print(s)

