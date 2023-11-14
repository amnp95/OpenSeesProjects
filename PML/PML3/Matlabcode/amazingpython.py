# %%
# read a text file strings from file K.txt
# %%
import numpy as np
K = np.loadtxt("K.txt", dtype=str, delimiter=',')
K.shape
f = open("syntax.txt", "w")
for i in range(K.shape[0]):
    for j in range(K.shape[1]-1):
        KString = K[i][j].item()
        KString = KString.replace(" ", "")
        KString = KString.replace("b^2", "b*b")
        KString = KString.replace("a^2", "a*a")
        KString = KString.replace("L_PML_x^2", "L_PML_x*L_PML_x")
        KString = KString.replace("L_PML_y^2", "L_PML_y*L_PML_y")
        KString = KString.replace("nu^2", "nu*nu")
        KString = KString.replace("yj^2", "yj*yj")
        KString = KString.replace("xi^2", "xi*xi")
        KString = KString.replace("*", " * ")
        KString = KString.replace("/", " / ")
        KString = KString.replace("+", " + ")
        KString = KString.replace("-", " - ")
        f.write(f"\tK[{i*11 + j}] = " + KString + ";\n")
        j+=1
f.close()
# %%
import numpy as np
M = np.loadtxt("M.txt ", dtype=str, delimiter=',')
f = open("syntax.txt", "w")
for i in range(M.shape[0]):
    for j in range(M.shape[1]-1):
        # erase all the spaces
        # M[i][j] = M[i][j].replace(" ", "")
        # change numpy.str index to python str
        MString = M[i][j].item()
        MString = MString.replace(" ", "")
        MString = MString.replace("b^2", "b*b")
        MString = MString.replace("a^2", "a*a")
        MString = MString.replace("L_PML_x^2", "L_PML_x*L_PML_x")
        MString = MString.replace("L_PML_y^2", "L_PML_y*L_PML_y")
        MString = MString.replace("nu^2", "nu*nu")
        MString = MString.replace("yj^2", "yj*yj")
        MString = MString.replace("xi^2", "xi*xi")
        MString = MString.replace("*", " * ")
        MString = MString.replace("/", " / ")
        MString = MString.replace("+", " + ")
        MString = MString.replace("-", " - ")
        f.write(f"\tM[{i*11 + j}] = " + MString + ";\n")
        j+=1
f.close()
# %%
import numpy as np
C = np.loadtxt("C.txt", dtype=str, delimiter=',')
C.shape
f = open("syntax.txt", "w")
for i in range(C.shape[0]):
    for j in range(C.shape[1]-1):
        CString = C[i][j].item()
        CString = CString.replace(" ", "")
        CString = CString.replace("b^2", "b*b")
        CString = CString.replace("a^2", "a*a")
        CString = CString.replace("a^3", "a*a*a")
        CString = CString.replace("b^3", "b*b*b")
        CString = CString.replace("L_PML_x^2", "L_PML_x*L_PML_x")
        CString = CString.replace("L_PML_y^2", "L_PML_y*L_PML_y")
        CString = CString.replace("nu^2", "nu*nu")
        CString = CString.replace("yj^2", "yj*yj")
        CString = CString.replace("xi^2", "xi*xi")
        CString = CString.replace("*", " * ")
        CString = CString.replace("/", " / ")
        CString = CString.replace("+", " + ")
        CString = CString.replace("-", " - ")
        f.write(f"\tC[{i*11 + j}] = " + CString + ";\n")
        j+=1
f.close()
# %%
