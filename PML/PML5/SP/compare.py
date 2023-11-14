# %%
import numpy as np
import matplotlib.pyplot as plt



Mmat = np.loadtxt('MMatlab.txt')
Cmat = np.loadtxt('CMatlab.txt')
Kmat = np.loadtxt('KMatlab.txt')

# %%
Mop = np.loadtxt('MOpenSees.txt', delimiter=",")
Cop = np.loadtxt('COpenSees.txt', delimiter=",")
Kop = np.loadtxt('KOpenSees.txt', delimiter=",")
# %%
Kop.shape
# %%
np.allclose(Mmat, Mop)

# %%
np.allclose(Cmat, Cop)
# %%
if (not (np.allclose(Kmat, Kop))):
    # show the difference
    print(np.where(Kmat != Kop))



# %%
def find_different_indices_with_tolerance(matrix1, matrix2, tolerance):
    # Create a boolean matrix indicating where the elements are not close within the given tolerance
    diff_matrix = ~np.isclose(matrix1, matrix2, atol=tolerance)
    
    # Find the indices where the elements are different
    different_indices = np.transpose(np.nonzero(diff_matrix))
    
    return different_indices

# Example matrices
matrix1 = np.array([[1.1, 2.2, 3.3],
                    [4.4, 5.5, 6.6]])

matrix2 = np.array([[1.12, 2.21, 3.32],
                    [4.43, 5.53, 6.64]])

# Tolerance for comparison
tolerance = 0.01

# Find indices where elements are different within the tolerance
different_indices = find_different_indices_with_tolerance(Kop, Kmat, tolerance)

print("Indices with different elements within the tolerance:")
print(different_indices)
# %%
