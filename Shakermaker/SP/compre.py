# %%
import numpy as np
import matplotlib.pyplot as plt
import os 

# list all file with name "G**.mat"
files = os.listdir("./PML")
files = [f for f in files if f.endswith('.mat') and f.startswith('M')]
files.sort()
# %%
tags = []
for filename in files:
    data = np.loadtxt("./PML/" + filename)
    # delete the first char and the last for '.mat' and convert to int
    tag = int(filename[1:-4])
    # chek if the G matrix is not all zero with tolerance 1e-9
    if np.all(np.abs(data) < 1e-9):
        # print('G matrix is all zero')
        tags.append(tag) 
len(tags)  
len(files)  
#  save tags to file
np.savetxt('tags.txt', tags, fmt='%d')

# delete all the files in the ./PML directory
# %%
files = os.listdir("./PML")
for filename in files:
    os.remove("./PML/" + filename)
# %% 

