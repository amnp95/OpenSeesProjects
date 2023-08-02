#%%
import numpy as np
import matplotlib.pyplot as plt


# %%
G = np.loadtxt('G.txt', delimiter=",")
K = np.loadtxt('K.txt', delimiter=",")
M = np.loadtxt('M.txt', delimiter=",")
C = np.loadtxt('C.txt', delimiter=",")
# %%
# plt.imshow(K)       
# np.allclose(K, K.T, atol=1e-8)   
keff = K + G + M  + C
# plt.spy(keff, markersize=0.5)
# plt.title('Keff pattern')
# plt.spy(G, color="red",markersize=0.5)
# plt.spy(M, color="green",markersize=0.5)
# plt.spy(M, color="blue",markersize=0.5)
plt.spy(G, color="green",markersize=0.5)
plt.spy(K, color="red",markersize=0.5)

print(np.allclose(C, C.T, atol=1e-4))
print(np.allclose(M, M.T, atol=1e-4))
print(np.allclose(G, G.T, atol=1e-4))
print(np.allclose(K, K.T, atol=1e-4))

# %%
G.shape
# %%
# check symmetry
np.allclose(G, G.T, atol=1e-8)
# %%
# store G as a sparse symmetric matrix
from scipy.sparse import csr_matrix
G = csr_matrix(G)
# show the sparsity pattern
plt.spy(G, markersize=0.5)

# %%
# compare the memory usage of the dense and sparse matrix
import sys
print(sys.getsizeof(keff))
print(sys.getsizeof(csr_matrix(keff)))

# %%
