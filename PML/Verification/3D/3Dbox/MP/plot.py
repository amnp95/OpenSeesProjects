# %%
import numpy as np
import matplotlib.pyplot as plt
# %%
# data = np.loadtxt('NodeDispwithout.out')
# plt.plot(data[:, 0], data[:, 1])
# # plt.plot(data[:, 0], data[:, 2])
# # plt.plot(data[:, 0], data[:, 3])
# # plt.plot(data[:, 0], data[:, 4])
# # plt.plot(data[:, 0], data[:, 5])
# # plt.plot(data[:, 0], data[:, 6])
# # plt.plot(data[:, 0], data[:, 7])
# # plt.plot(data[:, 0], data[:, 8])

# plt.xlabel('Time (s)')
# plt.ylabel('Z Displacement (m)')
# plt.legend(labels=['1', '2', '3', '4', '5', '6', '7', '8'])
# plt.ylim(-0.002, 0.002)
# plt.grid(linestyle='--')
# plt.title('without PML')

# %%
data = np.loadtxt('NodeDisp0.out')
plt.plot(data[:, 0], data[:, 1])
# plt.plot(data[:, 0], data[:, 2])
# plt.plot(data[:, 0], data[:, 3])
# plt.plot(data[:, 0], data[:, 4])
# plt.plot(data[:, 0], data[:, 5])
# plt.plot(data[:, 0], data[:, 6])
# plt.plot(data[:, 0], data[:, 7])
# plt.plot(data[:, 0], data[:, 8])

plt.xlabel('Time (s)')
plt.ylabel('Z Displacement (m)')
plt.legend(labels=['node 1', 'node 2', 'node 3', 'node 4', 'node 5', 'node 6', 'node 7', 'node 8'])
# plt.ylim(-0.0021, 0.001)
plt.xlim(right=1.1)
plt.title('with PML')
plt.grid(linestyle='--')
plt.show()

# %%
