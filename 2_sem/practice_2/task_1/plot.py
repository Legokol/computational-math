import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

result = np.loadtxt('pressure.txt')
x = np.linspace(0, 500, result.shape[1] - 1)

fig, ax = plt.subplots()
line, = ax.plot(x, result[0, 1:])
ax.grid()
ax.set_xlabel('x, м')
ax.set_ylabel('p, Па')


def animate(i):
    line.set_ydata(result[i, 1:])
    return line,


anim = animation.FuncAnimation(fig, animate, interval=25, repeat=False, frames=result.shape[0])

anim.save('pressure.gif')


# plt.plot(x, result[-1, 1:])
plt.show()
