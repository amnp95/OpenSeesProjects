# %%
import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('NodeDisp0.out')
plt.plot(data[:, 0], data[:, 1])
plt.plot(data[:, 0], data[:, 2])
plt.plot(data[:, 0], data[:, 3])
plt.xlabel('Time (s)')
plt.ylabel('Y Displacement (m)')
plt.legend(labels=['node 1', 'node 2', 'node 3'])
plt.ylim(-0.0002, 0.0003)
plt.xlim(right=1.1)
plt.title('with PML')
plt.grid(linestyle='--')
plt.show()

# # %%
# import numpy as np
# data = np.loadtxt('force.dat')
# plt.plot(np.arange(0, data.shape[0]*0.001,0.001), data)
# plt.grid(linestyle='--')
# plt.xlabel('Time (s)')
# plt.ylabel('Force (N)')
# plt.title('Load Time History')
# plt.show()

# %%
# fig, axes = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(12, 4))

# data = np.loadtxt('NodeDisp0.out')
# axes[0].plot(data[:, 0], data[:, 1], label='PML2D_3')
# axes[1].plot(data[:, 0], data[:, 2], label='PML2D_3')
# axes[2].plot(data[:, 0], data[:, 3], label='PML2D_3')


# data = np.loadtxt('/home/amin/Projects/Github/OpenSeesProjects/PML/Verification/2D/2Dfield/SingleCore/NodeDisp0.out')
# axes[0].plot(data[:, 0], data[:, 1], label='PML2D')
# axes[1].plot(data[:, 0], data[:, 2], label='PML2D')
# axes[2].plot(data[:, 0], data[:, 3], label='PML2D')



# axes[0].set_xlabel('Time (s)')
# axes[0].set_ylabel('Y Displacement (m)')
# axes[0].set_ylim(-0.0002, 0.0003)
# axes[0].grid(linestyle='--')
# axes[0].legend()
# axes[0].set_title('Node 1')

# axes[1].set_xlabel('Time (s)')
# axes[1].set_ylabel('Y Displacement (m)')
# axes[1].grid(linestyle='--')
# axes[1].set_title('Node 2')

# axes[2].set_xlabel('Time (s)')
# axes[2].set_ylabel('Y Displacement (m)')
# axes[2].grid(linestyle='--')
# axes[2].set_title('Node 3')




# %%
