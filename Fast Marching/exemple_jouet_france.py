
# In[2]:


import numpy as np
import matplotlib.pyplot as plt
n = 100
neigh = np.array([[1, -1, 0, 0], [0, 0,  1, -1]])
boundary = lambda x: np.mod(x,n)
ind2sub1 = lambda k: [int( (k-np.fmod(k,n))/n ), np.fmod(k,n)]
sub2ind1 = lambda u: int( u[0]*n + u[1] )
Neigh = lambda k,i: sub2ind1(boundary(ind2sub1(k) + neigh[:,i]))


# Check that these functions are indeed bijections.

# In[7]:


print( ind2sub1( sub2ind1([13, 27]) ) )
print( sub2ind1( ind2sub1(134) ) )


# Dikstra Algorithm
# -----------------
# The Dijkstra algorithm compute the geodesic distance on a graph.
# We use here a graph whose nodes are the pixels, and whose edges defines the
# usual 4-connectity relationship.
# 
# 
# In the following, we use the notation $i \sim j$ to indicate that an
# index $j$ is a neighbor of $i$ on the graph defined by the discrete
# grid.
# 
# 
# The metric $W(x)$. We use here a constant metric.



# In[13]:

n = 40
neigh = np.array([[1, -1, 0, 0], [0, 0,  1, -1]])
boundary = lambda x: np.mod(x,n)
ind2sub1 = lambda k: [int( (k-np.fmod(k,n))/n ), np.fmod(k,n)]
sub2ind1 = lambda u: int( u[0]*n + u[1] )
Neigh = lambda k,i: sub2ind1(boundary(ind2sub1(k) + neigh[:,i]))
extract   = lambda x,I: x[I]
extract1d = lambda x,I: extract(x.flatten(),I)


# In[20]:


def init(x0 = [n/2, n/2], W = np.ones( (n,n) ) ):
    I = [sub2ind1(x0)]
    D = np.ones( (n,n) ) + np.inf
    u = ind2sub1(I)
    D[u[0],u[1]] = 0
    S = np.zeros( (n,n) )
    S[u[0],u[1]] = 1
    
    return [W, D, S, I]

def classic_step(W, D, S, I):
    j = np.argsort(extract1d(D, I))
    if len(j) == 0:
        end = True
        return D,S,I, [], end
    j = j[0]
    i = I.pop(j)
    u = ind2sub1(i);
    S[u[0], u[1]] = -1
    J = []
    for k in np.arange(0, 4):
        j = Neigh(i, k)
        if extract1d(S, j) != -1:
            # add to the list of point to update
            J.append(j)
            if extract1d(S, j) == 0:
                # add to the front
                u = ind2sub1(j)
                S[u[0], u[1]] = 1
                I.append(j)
    return [np.copy(D), np.copy(S), I, J, False]


# In[21]:


def update_D_norme1(W, D, J):
    #Update distance by looking at all available places
    for j in J:
        DNeigh = lambda D, k: extract1d(D, Neigh(j, k))
        dx = min(DNeigh(D,0), DNeigh(D,1))
        dy = min(DNeigh(D,2), DNeigh(D,3))
        u = ind2sub1(j)
        w = extract1d(W,j);
        D[u[0],u[1]] = min(dx + w, dy + w)
    return D


# In[22]:


def dijkstra(update_norme, step, x0 = [n/2, n/2], W = np.ones( (n,n) ) ):
    [W, D, S, I] = init(x0, W)
    iteration = 0
    J = []
    end = False
    while True:
        [D, S, I, J, end] = step(W, D, S, I)
        if end:
            break
        else:
            D = update_norme(W, D, J)

        # iteration update
        if iteration % 400 == 0:
            plt.figure()
            plt.imshow(D)
            plt.set_cmap('jet')
            plt.show()
        iteration += 1
        start = False
    plt.figure(1)
    plt.imshow(D)
    plt.show()
    return D


def update_D_norme2(W, D, J):
    for j in J:
        DNeigh = lambda D, k: extract1d(D, Neigh(j, k))
        dx = min(DNeigh(D, 0), DNeigh(D, 1))
        dy = min(DNeigh(D, 2), DNeigh(D, 3))
        u = ind2sub1(j)
        w = extract1d(W, j)
        Delta = 2 * w ** 2 - (dx - dy) ** 2
        if (Delta >= 0):
            D[u[0], u[1]] = (dx + dy + np.sqrt(Delta)) / 2
        else:
            D[u[0], u[1]] = min(dx + w, dy + w)
    return D

def boundaries_control(x):
    inside_boundary = ((x[0])//n == 0) and ((x[1])//n == 0)
    return np.mod(x,n), inside_boundary


def newNeigh(pixel_initial,i):
    neighbor, inside_boundary = boundaries_control(ind2sub1(pixel_initial) + neigh[:,i])
    return sub2ind1(neighbor), inside_boundary


def newDNeigh(D,k, j):
    n, b = newNeigh(j,k)
    if b==False:
        return np.inf
    else:
        return extract1d(D,n)


def update_D_norme2_with_boundaries(W, D, J):
    for j in J:
        n1 = newDNeigh(D, 0, j)
        n2 = newDNeigh(D, 1, j)
        n3 = newDNeigh(D, 2, j)
        n4 = newDNeigh(D, 3, j)
        dx = min(n1, n2)
        dy = min(n3, n4)
        u = ind2sub1(j)
        w = extract1d(W, j)
        Delta = 2 * w ** 2 - (dx - dy) ** 2
        if (Delta >= 0):
            D[u[0], u[1]] = (dx + dy + np.sqrt(Delta)) / 2
        else:
            D[u[0], u[1]] = min(dx + w, dy + w)
    return D


def step_with_boundaries(W, D, S, I):
    j = np.argsort(extract1d(D, I))
    if len(j) == 0:
        end = True
        return D,S,I, [], end
    j = j[0]
    i = I.pop(j)
    u = ind2sub1(i);
    S[u[0], u[1]] = -1
    J = []
    for k in np.arange(0, 4):
        j, inside_boundary = newNeigh(i, k)
        if inside_boundary:
            if extract1d(S, j) != -1:
                # add to the list of point to update
                J.append(j)
                if extract1d(S, j) == 0:
                    # add to the front
                    u = ind2sub1(j)
                    S[u[0], u[1]] = 1
                    I.append(j)
    return [np.copy(D), np.copy(S), I, J, False]


def Geval(X, G):
    x, y = X
    xinf = np.int(x)
    xsup = xinf + 1
    yinf = np.int(y)
    ysup = yinf + 1
    
    return (G[:, xinf, yinf] * (xsup - x) * (ysup - y) +
            G[:, xinf, ysup] * (xsup - x) * (y - yinf) +
            G[:, xsup, yinf] * (x - xinf) * (ysup - y) +
            G[:, xsup, ysup] * (x - xinf) * (y - yinf) )


# In[58]:

import os
print("Current working directory: ", os.getcwd())
from PIL import Image
image_france = Image.open("carte_france_relief.jpg")
carte_france = np.array(image_france)


# In[73]:

# Paris
x0 = np.array([55,25])
# Nice
x1 = np.array([87, 78])
# Montpellier
x2 = np.array([62,80])
# Toulouse
x3 = np.array([40,80])
# Brest
x4 = np.array([1, 28])
# Strasbourg
x5 = np.array([87,25])
# Bordeaux
x6 = np.array([25,62])
# Lille
x7 = np.array([55,5])

cities = ["Nice", "Montpellier", "Toulouse", "Brest", "Strasbourg", "Bordeaux", "Lille"]


# In[74]:


carte_france.shape
n = 100
W = carte_france[:,:,0] + np.random.normal(0,1)**2
plt.imshow(W)
plt.set_cmap('gray')
plt.plot(x0[0], x0[1], '.r', markersize=20, label = "Paris")
plt.plot(x1[0], x1[1], '.g', markersize=20, label = "Nice")
plt.plot(x2[0], x2[1], '.y', markersize=20, label="Montpellier")
plt.plot(x3[0], x3[1], '.m', markersize=20, label="Toulouse")
plt.plot(x4[0], x4[1], '.c', markersize=20, label="Brest")
plt.plot(x5[0], x5[1], '.', markersize=20, label="Strasbourg")
plt.plot(x6[0], x6[1], '.', markersize=20, label="Bordeaux")
plt.plot(x7[0], x7[1], '.', markersize=20, label="Lille")
plt.legend()
plt.show()
X = np.array([x1,x2,x3,x4,x5,x6,x7])
# W = np.concatenate((255*np.ones((7,W.shape[1])), W), axis=0)
D = np.full(W.shape, np.infty)
S = np.zeros(W.shape)



# In[83]:


D_france = dijkstra(update_D_norme2_with_boundaries, step_with_boundaries, [x0[1],x0[0]], W)


# In[84]:


G0 = np.gradient(D_france.T)
d = np.sqrt(np.sum(np.array(G0)**2, axis=0))
U = np.zeros((2,n,n))
U[0,:,:] = d+10
U[1,:,:] = d+10
G0 = np.array(G0)
G = G0 / U
# plt.imshow(D_france[:20,:20])
# plt.show()


# In[113]:


# x0 = np.array([50,20])
# x1 = np.array([75, 40])
nb_ville = 0
h = plt.figure()
h = plt.plot(x0[0], x0[1], '.', markersize=20, label="Paris")
for x1 in X:
    city = cities[nb_ville]
    nb_ville += 1
    gamma = [x1]
    tau = 0.1
    itermax = 100000
    niter = 0
    while np.linalg.norm(gamma[-1] - x0) / np.linalg.norm(x0) > 0.05 and niter < itermax:
        g = Geval(gamma[-1], G)
        gamma.append( gamma[-1] - tau*g )
        niter +=1
        # if niter%1000==0:
            # print(round(niter/itermax*100), "%")
    
    W_plus = np.ones((W.shape[0], W.shape[1]+50)) * 255
    W_plus[:W.shape[0], :W.shape[1]] = np.copy(W)
    plt.imshow(W_plus) 
    plt.set_cmap('gray')
    h = plt.plot(np.array(gamma)[:,0], np.array(gamma)[:,1], '.b', linewidth=2)
    h = plt.plot(x1[0], x1[1], '.', markersize=20, label=city)
h = plt.legend()
h = plt.plot()




