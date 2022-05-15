import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

result = np.loadtxt('temperature.txt')
x = np.linspace(0, 1, result.shape[1] - 1)

fig, ax = plt.subplots()
line, = ax.plot(x, result[0, 1:])
ax.grid()
ax.set_xlabel('x, м')
ax.set_ylabel('T, К')


def animate(i):
    line.set_ydata(result[i, 1:])
    return line,


anim = animation.FuncAnimation(fig, animate, interval=50, repeat=False, frames=result.shape[0])

anim.save('temperature.gif')

# plt.plot(x, result[-1, 1:])
plt.show()
