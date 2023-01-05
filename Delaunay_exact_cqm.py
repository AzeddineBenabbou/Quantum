# use of exact solver for debugging
from dimod import ConstrainedQuadraticModel, Binaries, ExactCQMSolver, quicksum, Binary
from dwave.system import LeapHybridCQMSampler
import numpy as np

# instiate example from dwave documentation
# cqm = ConstrainedQuadraticModel()
# x, y, z = Binaries(['x', 'y', 'z'])
# cqm.set_objective(x*y + 2*y*z)
# cqm.add_constraint(x*y == 1, label='constraint_1')
# sampleset = ExactCQMSolver().sample_cqm(cqm)
# print(sampleset)

# define geometry of points
xy = np.array([0., 0., 1., 0., 0.9, 1., 0., 1., 0.6, 0.7, 0.25, 0.25, 0.5, 0.5, 0.25, 0.75])
border_edges = [0, 1, 1, 2, 2, 3, 3, 0]
# n is number of points
n = int(xy.size/2)
# define lifting map
z = [pow(xy[2*i], 2)+pow(xy[2*i+1], 2) for i in range(n)]
# define triangles
nb_var = int((n*(n-1)*(n-2))/6)
# dictionnary for internal vertices
vertex_to_triangle = {i: set() for i in range(n)}
var_to_triangle = np.zeros(3*nb_var, np.int8)
# define necessary info
var_index = int(0)
for i in range(n-2):
    for j in range(i+1, n-1):
        for k in range(j+1, n):
            vertex_to_triangle[i].add(var_index)
            vertex_to_triangle[j].add(var_index)
            vertex_to_triangle[k].add(var_index)
            var_to_triangle[3*var_index] = i
            var_to_triangle[3*var_index+1] = j
            var_to_triangle[3*var_index+2] = k
            var_index = var_index+1


# define binary variables
# labels = ['t00', 't01', 't02', 't03', 't04', 't05', 't06', 't07', 't08',
#           't09', 't10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19']
labels=[]
for i in range(nb_var):
    if (i <10):
        labels.append(f't0{i}')
    else:
        labels.append(f't{i}')

print(labels)


bin_variables = [Binary(labels[i]) for i in range(nb_var)]

print(bin_variables)

# define volumes
vol = np.zeros(nb_var)
for i in range(nb_var):
    v0 = var_to_triangle[3*i]
    v1 = var_to_triangle[3*i+1]
    v2 = var_to_triangle[3*i+2]
    # area
    x0 = xy[2*v0]
    y0 = xy[2*v0+1]
    x1 = xy[2*v1]
    y1 = xy[2*v1+1]
    x2 = xy[2*v2]
    y2 = xy[2*v2+1]
    # area=x0(y1-y2)-x1(y0-y2)+x2(y0-y1) (i ignore 2)
    area = abs(x0*(y1-y2)-x1*(y0-y2)+x2*(y0-y1))
    # volume of prism under the lifting map = sum(z0+z1+z2)xarea (i ingnre 3)
    vol[i] = area*(z[v0]+z[v1]+z[v2])
    # vol[i]=-area*((x0+x1+x2)*(x0+x1+x2)+(y0+y1+y2)*(y0+y1+y2))
    print(v0, ' ', v1, ' ', v2, ' ', vol[i])

# define volumes
border_list = [[0, 1], [1, 0], [1, 2], [2, 1], [0, 3], [3, 0], [2, 3], [3, 2]]
longest_edge = np.zeros(nb_var)
radii = np.zeros(nb_var)
for i in range(nb_var):
    v0 = var_to_triangle[3*i]
    v1 = var_to_triangle[3*i+1]
    v2 = var_to_triangle[3*i+2]
    # area
    x0 = xy[2*v0]
    y0 = xy[2*v0+1]
    x1 = xy[2*v1]
    y1 = xy[2*v1+1]
    x2 = xy[2*v2]
    y2 = xy[2*v2+1]
    # area=x0(y1-y2)-x1(y0-y2)+x2(y0-y1) (i ignore 2)
    area = 0.5*abs(x0*(y1-y2)-x1*(y0-y2)+x2*(y0-y1))

    # a, b, c are lenghts of edges 01, 12 and 20 respectively
    a = ((x1-x0)**2+(y1-y0)**2)**0.5
    b = ((x2-x1)**2+(y2-y1)**2)**0.5
    c = ((x2-x0)**2+(y2-y0)**2)**0.5
    # radius abc/4*area
    radii[i] = (a*b*c)/(4*area)
    #
    edge = [v0, v1]
    longest = 0
    if (border_list.count(edge) == 0):
        longest = a

    edge = [v1, v2]
    if (border_list.count(edge) == 0 and b > longest):
        longest = b

    edge = [v2, v0]
    if (border_list.count(edge) == 0 and c > longest):
        longest = c

    longest_edge[i] = longest
    # vol[i]=-area*((x0+x1+x2)*(x0+x1+x2)+(y0+y1+y2)*(y0+y1+y2))
    print(v0, ' ', v1, ' ', v2, ' ', radii[i], ' ', longest)

#


# instantiate a cqm
cqm = ConstrainedQuadraticModel()

# compute triangles for border edges
_ones = set()
_zeros = set()
# add constraints
# constraint on border edges
#this method is not sufficient
# for i in range(int(len(border_edges)/2)):
#     v0 = border_edges[2*i]
#     v1 = border_edges[2*i+1]
#     s0 = vertex_to_triangle[v0]
#     s1 = vertex_to_triangle[v1]
#     s = s0 & s1
#     #cqm.add_constraint(quicksum(bin_variables[j] for j in s)==1, f'border_edge_{v0}_{v1}')
#     t_min = -1
#     vol_min = 0
#     for t in s:
#         vol_t = vol[t]
#         if (t_min == -1):
#             t_min = t
#             vol_min = vol_t
#         elif (vol_min > vol_t):
#             vol_min = vol_t
#             t_min = t
#         # set 1
#         bin_variables[t_min] = 1

#     for t in s:
#         if (t != t_min):
#             # set 0 for others
#             bin_variables[t] = 0

# constraint on border edges
# with this method we use the circumsphere criteria
for i in range(int(len(border_edges)/2)):
    a = border_edges[2*i]
    b = border_edges[2*i+1]
    s0 = vertex_to_triangle[a]
    s1 = vertex_to_triangle[b]
    s = list(s0 & s1)
    ax = xy[2*a]
    ay = xy[2*a+1]
    bx = xy[2*b]
    by = xy[2*b+1]
    del_t=-1
    print("border edge ", a, ' ', b)
    for k in range(len(s)):
        t=s[k]
        v0 = var_to_triangle[3*t]
        v1 = var_to_triangle[3*t+1]
        c = var_to_triangle[3*t+2]
        #check for c
        if (c==a or c==b):
            c=v1
        if (c==a or c==b):
            c=v0
        
        print('3rd vertex ', c)
        cx = xy[2*c]
        cy = xy[2*c+1]
        #compute the center and radius of the circumsphere
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        centerx = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
        centery = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
        radius_sq=(centerx-ax)*(centerx-ax)+(centery-ay)*(centery-ay)

        for l in range(len(s)):
            other_t = s[l]
            if (other_t==t):
                continue
            v0 = var_to_triangle[3*other_t]
            v1 = var_to_triangle[3*other_t+1]
            c = var_to_triangle[3*other_t+2]
            #check for c
            if (c==a or c==b):
                c=v1
            if (c==a or c==b):
                c=v0
            cx = xy[2*c]
            cy = xy[2*c+1]
            dist=(centerx-cx)*(centerx-cx)+(centery-cy)*(centery-cy)
            if (dist<radius_sq):
                t=-1
                break
        
        if (t!=-1):
            #this is a Delaunay triangle
            del_t=t
            print('delaunay found here ')
            break
    
    #
    bin_variables[del_t] = 1
    for t in s:
        if (t==del_t):
            continue
        bin_variables[t]=0







# update border edges
for i in range(n-1):
    s0 = vertex_to_triangle[i]
    for j in range(i+1, n):
        edge = [i, j]
        if (border_list.count(edge) == 1):
            print(edge, ' exists ')
            continue

        s1 = vertex_to_triangle[j]
        s = list(s0 & s1)
        sum = 0
        for k in s:
            if (bin_variables[k] == 1):
                sum = sum+1

        if (sum == 2):
            print('sum==2 case for ', edge)
            for k in s:
                if (bin_variables[k] != 1):
                    bin_variables[k] = 0

# update border edges
objective4 = 0
for i in range(n-1):
    s0 = vertex_to_triangle[i]
    for j in range(i+1, n):
        edge = [i, j]
        if (border_list.count(edge) == 1):
            print(edge, ' exists ')
            continue

        s1 = vertex_to_triangle[j]
        s = list(s0 & s1)
        sum_1 = 0
        sum_0 = 0
        for k in s:
            if bin_variables[k] == 1:
                sum_1 = sum_1+1
            elif bin_variables[k] == 0:
                sum_0 = sum_0+1
            
        
        doit=False
        if sum_1 == 0 and sum_0 < len(s):
            print('190 constraint for edge ', edge)
            #cqm.add_constraint(quicksum(bin_variables[l] for l in s) == 2 or quicksum(bin_variables[l] for l in s) == 0, f'internal_edge_{i}_{j}')
            cqm.add_constraint(quicksum(bin_variables[k] for k in s)-2*quicksum(bin_variables[s[k]]*bin_variables[s[l]]   for k in range(len(s)-1) for l in range(k+1, len(s)))==0, f'internal_edge_{i}_{j}')
            print(quicksum(bin_variables[l] for l in s))
            doit=True
        if sum_1 == 1:
            print('193 constraint for edge ', edge)
            print(quicksum(bin_variables[l] for l in s if bin_variables[l]!=0 and bin_variables[l]!=1))
            cqm.add_constraint(quicksum(bin_variables[l] for l in s if bin_variables[l]!=0 and bin_variables[l]!=1)==1, f'internal_edge_{i}_{j}')
            doit=True
        
        if (doit==True):
            print('objective for ', edge, ' ', quicksum(vol[k] for k in s))
            objective4 = objective4+quicksum(vol[k] for k in s)*quicksum(bin_variables[k] for k in s)
        

#print(objective4)

for k in range(nb_var):
    value=-1
    if (bin_variables[k]==1):
        value=1
    elif (bin_variables[k]==0):
        value=0
    print (var_to_triangle[3*k], ' ', var_to_triangle[3*k+1], ' ',  var_to_triangle[3*k+2], ' ', value)


# # define the objective function
objective1=quicksum(vol[i]*bin_variables[i] for i in range(nb_var))

# #
# objective2=-quicksum(bin_variables[i] for i in range(nb_var))/nb_var

# #
# objective3=quicksum(radii[i]*bin_variables[i] for i in range(nb_var))

# #
# objective4 = 0
# # constraint on internal edges
# for i in range(n-1):
#     s0 = vertex_to_triangle[i]
#     for j in range(i+1, n):
#         edge = [i, j]
#         if (border_list.count(edge) == 1):
#             print(edge, ' exists ')
#             continue

#         s1 = vertex_to_triangle[j]
#         s = list(s0 & s1)

#         #cqm.add_constraint(quicksum(bin_variables[k] for k in s)-2*quicksum(bin_variables[s[k]]*bin_variables[s[l]]   for k in range(len(s)-1) for l in range(k+1, len(s)))==0, f'internal_edge_{i}_{j}')
#         print('for inernal edge ', edge, ' s is ', s)
#         cqm.add_constraint(quicksum(bin_variables[k] for k in s) == 0 or quicksum(
#             bin_variables[k] for k in s) == 2, f'internal_edge_{i}_{j}')
#         objective4 = objective4+quicksum(vol[k]*bin_variables[k] for k in s)


#constraint of incident triangles of internal vertices
for i in range(4, n):
    s=vertex_to_triangle[i]
    cqm.add_constraint(quicksum(bin_variables[j] for j in s)>=3, f'internal_vertex_{i}')


# # set objectives
#print(objective1)
cqm.set_objective(objective1)

# # solve

#sampleset = ExactCQMSolver().sample_cqm(cqm)
cqm_sampler = LeapHybridCQMSampler()
sampleset=cqm_sampler.sample_cqm(cqm, label ='Delaunay_general_case')
all_samples=sampleset.samples
#print(all_samples)
feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
fsamples_agr=feasible_sampleset.aggregate()
print('---------------------------------------------FEASIBLE-------------------------------------------------------------')
delaunay_trgl=[3,8,11,15,16,19]
print('feasible_sampleset')
for s in feasible_sampleset:
    #if (quicksum(s[labels[i]] for i in delaunay_trgl)==6):
    print(s)

print(feasible_sampleset.first)
