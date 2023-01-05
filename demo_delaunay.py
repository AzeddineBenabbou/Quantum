#basic algorithm to compute the delaunay triangulation of a set points in 2d 
# I use CQM
from dimod import Binary, ConstrainedQuadraticModel, quicksum
from dwave.system import LeapHybridCQMSampler
import numpy as np

#define geometry of points
xy=np.array([0., 0., 1., 0., 0.9, 1., 0., 1., 0.6, 0.7, 0.4, 0.3])
border_edges=[0, 1, 1, 2, 2, 3, 3, 0]
# n is number of points
n=int(xy.size/2)
#define lifting map
z=[pow(xy[2*i], 2)+pow(xy[2*i+1], 2) for i in range(n)]
#define triangles
nb_var=int((n*(n-1)*(n-2))/6)
#dictionnary for internal vertices
vertex_to_triangle={i:set() for i in range(n)}
var_to_triangle=np.zeros(3*nb_var, np.int8)
#define necessary info
var_index=int(0)
for i in range(n-2):
    for j in range(i+1, n-1):
        for k in range (j+1, n):
            vertex_to_triangle[i].add(var_index)
            vertex_to_triangle[j].add(var_index)
            vertex_to_triangle[k].add(var_index)
            var_to_triangle[3*var_index]=i
            var_to_triangle[3*var_index+1]=j
            var_to_triangle[3*var_index+2]=k
            var_index=var_index+1


# define binary variables
bin_variables = [Binary(f't_{i}') for i in range(nb_var)]

#define volumes
vol=np.zeros(nb_var)
for i in range(nb_var):
    v0=var_to_triangle[3*i]
    v1=var_to_triangle[3*i+1]
    v2=var_to_triangle[3*i+2]
    #area
    x0=xy[2*v0]
    y0=xy[2*v0+1]
    x1=xy[2*v1]
    y1=xy[2*v1+1]
    x2=xy[2*v2]
    y2=xy[2*v2+1]
    #area=x0(y1-y2)-x1(y0-y2)+x2(y0-y1) (i ignore 2)
    area=abs(x0*(y1-y2)-x1*(y0-y2)+x2*(y0-y1))
    #volume of prism under the lifting map = sum(z0+z1+z2)xarea (i ingnre 3)
    vol[i]=area*(z[v0]+z[v1]+z[v2])
    #vol[i]=-area*((x0+x1+x2)*(x0+x1+x2)+(y0+y1+y2)*(y0+y1+y2))
    print(v0, ' ', v1, ' ', v2, ' ', vol[i])

# print('----------------------------------------------------------------------------------------')
# for i in range(int(len(border_edges)/2)):
#     v0=border_edges[2*i]
#     v1=border_edges[2*i+1]
#     s0=vertex_to_triangle[v0]
#     s1=vertex_to_triangle[v1]
#     s=s0 & s1
#     print('edge ', v0, ' ', v1)
#     for j in s:
#         print('var ', j, ' ', var_to_triangle[3*j], ' ', var_to_triangle[3*j+1],' ', var_to_triangle[3*j+2])

# print('---------------------------------------------------------------------------------------------------')    
# border_list=[[0,1],[1,0], [1,2], [2,1], [0, 3], [3,0]]
# bin_variables[2]=1
# bin_variables[11]=1
# bin_variables[13]=0
# bin_variables[15]=0
# for i in range(n-1):
#     s0=vertex_to_triangle[i]
#     for j in range(i+1, n):
#         edge=[i,j]
#         print(edge)
#         if (border_list.count(edge)==1):
#             print(' exists ')
#             continue


#         s1=vertex_to_triangle[j]
#         s=list(s0 & s1)
#         print('s ')
#         for jj in s:
#             print('var ', jj, ' ', var_to_triangle[3*jj], ' ', var_to_triangle[3*jj+1],' ', var_to_triangle[3*jj+2])

#         print(quicksum(bin_variables[k] for k in s)-2*quicksum(bin_variables[s[k]]*bin_variables[s[l]]   for k in range(len(s)-1) for l in range(k+1, len(s))))
        



        #cqm.add_constraint(quicksum(bin_variables[k] for k in s)-2*quicksum(bin_variables[s[k]]*bin_variables[s[l]]   for k in range(len(s)-1) for l in range(k+1, len(s)))==0, f'internal_edge_{i}_{j}')


# instantiate a cqm
cqm=ConstrainedQuadraticModel()

# define the objective function
objective=quicksum(vol[i]*bin_variables[i] for i in range(nb_var))

# set objectives
cqm.set_objective(objective)

# add constraints
#constraint on border edges
for i in range(int(len(border_edges)/2)):
    v0=border_edges[2*i]
    v1=border_edges[2*i+1]
    s0=vertex_to_triangle[v0]
    s1=vertex_to_triangle[v1]
    s=s0 & s1
    cqm.add_constraint(quicksum(bin_variables[j] for j in s)==1, f'border_edge_{i}')

#constraint on internal edges
#make a list of border edges
border_list=[[0,1],[1,0], [1,2], [2,1], [0, 3], [3,0]]
for i in range(n-1):
    s0=vertex_to_triangle[i]
    for j in range(i+1, n):
        edge=[i,j]
        if (border_list.count(edge)==1):
            print(edge, ' exists ')
            continue

        s1=vertex_to_triangle[j]
        s=list(s0 & s1)
        cqm.add_constraint(quicksum(bin_variables[k] for k in s)-2*quicksum(bin_variables[s[k]]*bin_variables[s[l]]   for k in range(len(s)-1) for l in range(k+1, len(s)))==0, f'internal_edge_{i}_{j}')

#constraint of incident triangles of internal vertices
for i in range(4, n):
    s=vertex_to_triangle[i]
    cqm.add_constraint(quicksum(bin_variables[j] for j in s)>=3, f'internal_vertex_{i}')




# # submit to a solver
# # CQM sampler
# # instantiate a cqm solver (sampler)
cqm_sampler = LeapHybridCQMSampler()

# # # #solve with the instance of the cqm solver
sampleset = cqm_sampler.sample_cqm(cqm, time_limit=25, label ='Demo Delaunay 2d 6 pts')
print("all samplest aggregated")
print(sampleset.aggregate())

# # filter results. Only feasible solutions.
print("feasible solutions aggregated")
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
fsamples_agr=feasible_sampleset.aggregate()
print(fsamples_agr)

# vertex_to_triangle[0].add(0)
# A = {1, 2, 3, 4, 5}
# B = {4, 5, 6, 7, 8}

# # use & operator
# # Output: {4, 5}
# C=A & B
# for e in C:
#     print(e)
# print(C)







# points=np.array([[0., 0.], [1, 0.], [0.9,1.], [0.,1.], [0.6, 0.7], [0.4, 0.3]])
# #tri=Delaunay(points[0:4], furthest_site=False, incremental=True, qhull_options=None)# unstable!
# tri=Delaunay(points)

# print('dalunay triangles')
# print(tri.simplices)
# plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
# plt.plot(points[:,0], points[:,1], 'o')
# plt.show()



# #define the edges (variables of current problem)
# # edges ={{0, 1}, {1, 2}, {2, 3} ...{n-2, n-1} {0, 2}, {1, 3} ... {n-3, n-1} .....{0, n-2}, {1, n-1}, {0, n-1}}
# # e1=[1 2]

# #define edges for the current problem
# edges=np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 0, 2, 1, 3, 2, 4, 3, 5, 0, 3, 1, 4, 2, 5, 0, 4, 1, 5, 0, 5])

# print('edges')
# for i in range(nb_var):
#     print(edges[2*i], edges[2*i+1])

# # define lifting map z=x^2+y^2
# z=[pow(xy[2*i], 2)+pow(xy[2*i+1], 2) for i in range(xy.size//2)]
# print('z=', z)

# #defin connectivity matrix nb_var x nb_var
# #define area matrixclear
# #define table of triangles
# m=np.zeros((nb_var, nb_var))
# a=np.zeros((nb_var, nb_var))
# vol=np.zeros((nb_var, nb_var))
# tri=np.zeros((nb_var, nb_var, 3))
# for i in range(nb_var):
#     v0=edges[2*i]
#     v1=edges[2*i+1]
#     for j in range(nb_var):
#         if (j!=i):
#             v2=edges[2*j]
#             v3=edges[2*j+1]
#             if ((v2==v0) | (v2==v1) | (v3==v0) | (v3==v1)):
#                 m[i][j]=1
#                 min_v=min(min(v0, v1), min(v2, v3))
#                 max_v=max(max(v0, v1), max(v2, v3))
#                 mid_v=v0
#                 for v in [v0, v1, v2, v3]:
#                     if (v!=min_v):
#                         if (v!=max_v):
#                             mid_v=v
#                             break
                        
#                 tri[i][j]=[min_v, mid_v, max_v]

#                 # print(' ', v0, ' ', v1, ' ', v2, ' ', v3)
#                 # print('min_v=', min_v, ' mid__v=', mid_v, 'max_v=', max_v)

#                 x0=xy[2*min_v]
#                 y0=xy[2*min_v+1]
#                 x1=xy[2*mid_v]
#                 y1=xy[2*mid_v+1]
#                 x2=xy[2*max_v]
#                 y2=xy[2*max_v+1]

#                 #area=x0(y1-y2)-x1(y0-y2)+x2(y0-y1) (i ignore 2)
#                 area=abs(x0*(y1-y2)-x1*(y0-y2)+x2*(y0-y1))
#                 a[i][j]=area
#                 #volume of prism under the lifting map = sum(z0+z1+z2)xarea (i ingnre 3)
#                 vol[i][j]=area*(z[min_v]+z[mid_v]+z[max_v])
               

# print('m=')
# print(m)

# # print('tri=')
# # print(tri)

# print('areas=')
# print(a)

# print('volumes=')
# print(vol)

# #total area for this example
# total_area=1.9






# define binary variables
# bin_variables = [Binary(f'{i}') for i in range(nb_var)]

# # # print what kind of variables we are using
# # # print(bin_variables[0])

# # # instantiate a cqm
# cqm=ConstrainedQuadraticModel()

# # # define the objective function
# objective1=quicksum(m[i][j]*vol[i][j]*bin_variables[i]*bin_variables[j] for i in range(nb_var) for j in range(i+1, nb_var))

# # print(objective1)

# # #objective2=quicksum(-1*bin_variables[i] for i in range(nb_var))


# # # set objectives
# cqm.set_objective(objective1)


# # # add constraints
# # #constraint 1= existing edges
# unknown_edges=np.ones(nb_var)
# unknown_edges[0]=unknown_edges[1]=unknown_edges[2]=unknown_edges[9]=0
# print("unknown edges")
# print(unknown_edges)
# cqm.add_constraint(quicksum((1-unknown_edges[i])*bin_variables[i] for i in range(nb_var))==4, 'existing edge constraint')

# # #constraint2 = sum of areas of incident triangles to internal edges = total area (for all edges in the solution)



# #cqm.add_constraint(quicksum(unknown_edges[i]*bin_variables[i]*quicksum(a[i][j] for j in range(nb_var)) for i in range(nb_var))==2.*total_area, 'total area')

# cqm.add_constraint(quicksum(unknown_edges[i]*bin_variables[i] for i in range(nb_var))==7, 'total edges')



# # #constraints on border edges
# # # cqm.add_constraint(quicksum(m[0][j]*bin_variables[0]*bin_variables[j] for j in range(nb_var))==2, 'border edge0 one triangle')
# # # cqm.add_constraint(quicksum(m[1][j]*bin_variables[1]*bin_variables[j] for j in range(nb_var))==2, 'border edge1 one triangle')
# # # cqm.add_constraint(quicksum(m[2][j]*bin_variables[2]*bin_variables[j] for j in range(nb_var))==2, 'border edge2 one triangle')
# # # cqm.add_constraint(quicksum(m[5][j]*bin_variables[5]*bin_variables[j] for j in range(nb_var))==2, 'border edge5 one triangle')
# # #constraints on internal edges
# # #cqm.add_constraint(quicksum(m[3][j]*bin_variables[3]*bin_variables[j] for j in range(nb_var))%4==0, 'border edge3 two triangles')
# # #cqm.add_constraint(quicksum(m[4][j]*bin_variables[4]*bin_variables[j] for j in range(nb_var))%4==0, 'border edge4 two triangles')

# # # print the number of variables 
# # # print("number of variables for this problem {}".format(len(cqm.variables)))

# # # submit to a solver
# # # CQM sampler
# # # instantiate a cqm solver (sampler)
# cqm_sampler = LeapHybridCQMSampler()

# # # # #solve with the instance of the cqm solver
# sampleset = cqm_sampler.sample_cqm(cqm, time_limit=5, label ='Demo Delaunay 2d 6 pts')
# print("all samplest aggregated")
# print(sampleset.aggregate())

# # # filter results. Only feasible solutions.
# print("feasible solutions aggregated")
# feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
# fsamples_agr=feasible_sampleset.aggregate()
# print(fsamples_agr)




  


