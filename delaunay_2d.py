#basic algorithm to compute the delaunay triangulation of a set points in 2d 
# I use CQM
from dimod import Binary, ConstrainedQuadraticModel, quicksum
from dwave.system import LeapHybridCQMSampler
import numpy as np

#define geometry of points
xy=np.array([0., 0., 1., 0., 0.9, 1., 0., 1.])
n=xy.size/2

#define the edges (variables of current problem)
# edges ={{0, 1}, {1, 2}, {2, 3} ...{n-2, n-1} {0, 2}, {1, 3} ... {n-3, n-1} .....{0, n-2}, {1, n-1}, {0, n-1}}
# e1=[1 2]
nb_var=int((n*n-n)/2)
#define edges for the current problem
edges=np.array([0, 1, 1, 2, 2, 3, 0, 2, 1, 3, 0, 3])

# define lifting map z=x^2+y^2
z=[pow(xy[2*i], 2)+pow(xy[2*i+1], 2) for i in range(xy.size//2)]
print('z=', z)

#defin connectivity matrix nb_var x nb_var
#define area matrix
#define table of triangles
m=np.zeros((nb_var, nb_var))
a=np.zeros((nb_var, nb_var))
vol=np.zeros((nb_var, nb_var))
tri=np.zeros((nb_var, nb_var, 3))
for i in range(nb_var):
    v0=edges[2*i]
    v1=edges[2*i+1]
    for j in range(nb_var):
        if (j!=i):
            v2=edges[2*j]
            v3=edges[2*j+1]
            if ((v2==v0) | (v2==v1) | (v3==v0) | (v3==v1)):
                m[i][j]=1
                min_v=min(min(v0, v1), min(v2, v3))
                max_v=max(max(v0, v1), max(v2, v3))
                mid_v=v0
                for v in [v0, v1, v2, v3]:
                    if (v!=min_v):
                        if (v!=max_v):
                            mid_v=v
                            break
                        
                tri[i][j]=[min_v, mid_v, max_v]

                # print(' ', v0, ' ', v1, ' ', v2, ' ', v3)
                # print('min_v=', min_v, ' mid__v=', mid_v, 'max_v=', max_v)

                x0=xy[2*min_v]
                y0=xy[2*min_v+1]
                x1=xy[2*mid_v]
                y1=xy[2*mid_v+1]
                x2=xy[2*max_v]
                y2=xy[2*max_v+1]

                #area=x0(y1-y2)-x1(y0-y2)+x2(y0-y1) (i ignore 2)
                area=abs(x0*(y1-y2)-x1*(y0-y2)+x2*(y0-y1))
                a[i][j]=area
                #volume of prism under the lifting map = sum(z0+z1+z2)xarea (i ingnre 3)
                vol[i][j]=area*(z[min_v]+z[mid_v]+z[max_v])
               

print('m=')
print(m)

# print('tri=')
# print(tri)

print('areas=')
print(a)

print('volumes=')
print(vol)

#total area for this example
total_area=1.9






# define binary variables
bin_variables = [Binary(f'e_{i}') for i in range(nb_var)]

# print what kind of variables we are using
# print(bin_variables[0])

# instantiate a cqm
cqm=ConstrainedQuadraticModel()

# define the objective function
objective1=quicksum(m[i][j]*vol[i][j]*bin_variables[i]*bin_variables[j] for i in range(nb_var) for j in range(i+1, nb_var))

print(objective1)

#objective2=quicksum(-1*bin_variables[i] for i in range(nb_var))


# set objectives
cqm.set_objective(objective1)


# add constraints
#constraint 1= existing edges
cqm.add_constraint(bin_variables[0]+bin_variables[1]+bin_variables[2]+bin_variables[5]==4, 'existing edge constraint')
#constraint2 = sum of incident areas of incident triangles to internal edges = total area (for all edges in the solution)
cqm.add_constraint(bin_variables[3]*quicksum(a[3][i] for i in range(nb_var))+ bin_variables[4]*quicksum(a[4][i] for i in range(nb_var))==2.*total_area, 'total area')

#constraints on border edges
# cqm.add_constraint(quicksum(m[0][j]*bin_variables[0]*bin_variables[j] for j in range(nb_var))==2, 'border edge0 one triangle')
# cqm.add_constraint(quicksum(m[1][j]*bin_variables[1]*bin_variables[j] for j in range(nb_var))==2, 'border edge1 one triangle')
# cqm.add_constraint(quicksum(m[2][j]*bin_variables[2]*bin_variables[j] for j in range(nb_var))==2, 'border edge2 one triangle')
# cqm.add_constraint(quicksum(m[5][j]*bin_variables[5]*bin_variables[j] for j in range(nb_var))==2, 'border edge5 one triangle')
#constraints on internal edges
#cqm.add_constraint(quicksum(m[3][j]*bin_variables[3]*bin_variables[j] for j in range(nb_var))%4==0, 'border edge3 two triangles')
#cqm.add_constraint(quicksum(m[4][j]*bin_variables[4]*bin_variables[j] for j in range(nb_var))%4==0, 'border edge4 two triangles')

# print the number of variables 
# print("number of variables for this problem {}".format(len(cqm.variables)))

# submit to a solver
# CQM sampler
# instantiate a cqm solver (sampler)
cqm_sampler = LeapHybridCQMSampler()

# # #solve with the instance of the cqm solver
sampleset = cqm_sampler.sample_cqm(cqm, time_limit=5, label ='Delaunay 2d with objective and constraints')
print("all samplest aggregated")
print(sampleset.aggregate())

# filter results. Only feasible solutions.
print("feasible solutions aggregated")
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
fsamples_agr=feasible_sampleset.aggregate()
print(fsamples_agr)




  

