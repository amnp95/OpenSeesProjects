# %%
import numpy as np
import matplotlib.pyplot as plt


def read(filename):
    lines = open(filename).readlines()
    # for each line ignore the latest string
    # and convert the rest to float
    return np.array([list(map(float, line.split()[:])) for line in lines], dtype=np.float64)


M = read("M.txt")
C = read("C.txt")
K = read("K.txt")
G = read("G.txt")
# check if M is symmetric
if np.allclose(M, M.T):
    print("M is symmetric")


if np.allclose(C, C.T):
    print("C is symmetric")

if np.allclose(K, K.T):
    print("K is symmetric")

if np.allclose(G, G.T):
    print("G is symmetric")



# # plt.spy(M)
# fig, ax = plt.subplots(1,1, figsize=(8, 8),sharex=True, sharey=True)
# ax.spy(K, markersize=1.5,color = "Blue")
# ax.spy(G, markersize=1.5,color = "Red")
# ax.spy(C, markersize=1.5,color = "Green")
# ax.spy(M, markersize=1.0,color = "Black")

a = K[0,0]/ G[0,0]
print(a)
if np.allclose(a*G, K):
    print("K = a*G")


if np.allclose(np.where(C != 0, 1, 0), np.where(K != 0, 1, 0)):
    print("same sparse pattern")
# %%
