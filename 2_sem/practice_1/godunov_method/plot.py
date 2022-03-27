import numpy as np
import matplotlib.pyplot as plt

t = np.load('t.npy')
x = np.load('x.npy')
solution = np.load('solution.npy')

i = 1500
# print(t[i])


fig, ax = plt.subplots(2, 2)

ax[0, 0].plot(x, solution[i, :, 0], linestyle='none', marker='.')
ax[0, 0].grid()
ax[0, 0].set_title('rho')

ax[0, 1].plot(x, solution[i, :, 1] / solution[i, :, 0], linestyle='none', marker='.')
ax[0, 1].grid()
ax[0, 1].set_title('u')

ax[1, 0].plot(x, solution[i, :, 2] / solution[i, :, 0], linestyle='none', marker='.')
ax[1, 0].grid()
ax[1, 0].set_title('e')

ax[1, 1].plot(x, solution[i, :, 3], linestyle='none', marker='.')
ax[1, 1].grid()
ax[1, 1].set_title('p')

plt.show()
