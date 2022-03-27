import numpy as np
from scipy.optimize import fsolve

# Параметры задачи
gamma = 5 / 3
L = 10
T = 0.02

# Начальные условия
pL = 10e5
pR = 1e5
rhoL = 13
rhoR = 1.3
vL = 0
vR = 0


# Вектор потоков
def f(u):
    return np.array([u[1], u[1] ** 2 / u[0] + u[3], u[1] / u[0] * (u[2] + u[3])])


# Функции для аппроксимации потоков
def a(rho):
    return 2 / (gamma + 1) / rho


def b(p):
    return (gamma - 1) / (gamma + 1) * p


def f_p(p, w):
    if p >= w[2]:
        return (p - w[2]) * np.sqrt(a(w[0]) / (p + b(w[2])))
    else:
        c = np.sqrt(gamma * w[2] / w[0])
        return 2 * c / (gamma - 1) * (np.power(p / w[2], (gamma - 1) / 2 / gamma) - 1)


# Аппроксимация потока
def calc_w(w_left, w_right):
    p0 = 0.5 * (w_left[2] + w_right[2])
    p = fsolve(lambda p: f_p(p, w_left) + f_p(p, w_right) + w_right[1] - w_left[1], p0)[0]
    u = 0.5 * (w_left[1] + w_right[1]) + 0.5 * (f_p(p, w_right) - f_p(p, w_left))
    if u >= 0:
        c_l = np.sqrt(gamma * w_left[2] / w_left[0])
        if p > w_left[2]:
            s_l = w_left[1] - c_l * np.sqrt((gamma + 1) / 2 / gamma * p / w_left[2] + (gamma - 1) / 2 / gamma)
            if s_l >= 0:
                return w_left
            else:
                rho = w_left[0] * (p / w_left[2] + (gamma - 1) / (gamma + 1)) / (
                        (gamma - 1) / (gamma + 1) * p / w_left[2] + 1)
                return np.array([rho, u, p])
        else:
            s_hl = w_left[1] - c_l
            s_tl = u - c_l * np.power(p / w_left[2], (gamma - 1) / 2 / gamma)
            if s_hl >= 0:
                return w_left
            elif s_tl <= 0:
                rho = w_left[0] * np.power(p / w_left[2], 1 / gamma)
                return np.array([rho, u, p])
            else:
                rho = w_left[0] * np.power(2 / (gamma + 1) + (gamma - 1) / (gamma + 1) / c_l * w_left[1],
                                           2 / (gamma - 1))
                u = 2 / (gamma + 1) * (c_l + (gamma - 1) / 2 * w_left[1])
                p = w_left[2] * np.power(2 / (gamma + 1) + (gamma - 1) / (gamma + 1) / c_l * w_left[1],
                                         2 * gamma / (gamma - 1))
                return np.array([rho, u, p])
    else:
        c_r = np.sqrt(gamma * w_right[2] / w_right[0])
        if p > w_right[2]:
            s_r = w_right[1] + c_r * np.sqrt(
                (gamma + 1) / 2 / gamma * p / w_right[2] + (gamma - 1) / 2 / gamma)
            if s_r <= 0:
                return w_right
            else:
                rho = w_right[0] * (p / w_right[2] + (gamma - 1) / (gamma + 1)) / (
                        (gamma - 1) / (gamma + 1) * p / w_right[2] + 1)
                return np.array([rho, u, p])
        else:
            s_hr = w_right[1] + c_r
            s_tr = u + c_r * np.power(p / w_right[2], (gamma - 1) / 2 / gamma)
            if s_hr <= 0:
                return w_right
            elif s_tr >= 0:
                rho = w_right[0] * np.power(p / w_right[2], 1 / gamma)
                return np.array([rho, u, p])
            else:
                rho = w_right[0] * np.power(2 / (gamma + 1) - (gamma - 1) / (gamma + 1) / c_r * w_right[1],
                                            2 / (gamma - 1))
                u = 2 / (gamma + 1) * (-c_r + (gamma - 1) / 2 * w_right[1])
                p = w_right[2] * np.power(2 / (gamma + 1) - (gamma - 1) / (gamma + 1) / c_r * w_right[1],
                                          2 * gamma / (gamma - 1))
                return np.array([rho, u, p])


# Переход от переменных u к w и обратно
def u_to_w(u):
    return np.array([u[0], u[1] / u[0], u[3]])


def w_to_u(w):
    return np.array([w[0], w[1] * w[0], (gamma - 1) * w[2], w[2]])


# Шаги по пространству и времени
h = 0.2
tau = 1e-5

solution = np.empty([1, int(2 * L / h), 4])

x = np.linspace(-L + 0.5 * h, L - 0.5 * h, int(2 * L / h))

for j in range(solution.shape[1]):
    if -L + (j + 0.5) * h < 0:
        solution[0][j] = [rhoL, rhoL * vL, pL / (gamma - 1), pL]
    else:
        solution[0][j] = [rhoR, rhoR * vR, pR / (gamma - 1), pR]

t = tau
time = np.array([t])
while t <= T:
    timeLayer = np.empty(solution[0].shape)
    w_current = u_to_w(solution[-1][0])

    f_left = f(w_to_u(calc_w(w_current, w_current)))
    for j in range(timeLayer.shape[0] - 1):
        w_right = u_to_w(solution[-1][j + 1])
        f_right = f(w_to_u(calc_w(w_current, w_right)))
        timeLayer[j, 0:3] = solution[-1, j, 0:3] + tau / h * (f_left - f_right)
        timeLayer[j, 3] = timeLayer[j, 2] * (gamma - 1)

        f_left = f_right
        w_current = w_right

    f_right = f(w_to_u(calc_w(w_current, w_right)))
    timeLayer[-1, 0:3] = solution[-1, -1, 0:3] + tau / h * (f_left - f_right)
    timeLayer[-1, 3] = timeLayer[-1, 2] * (gamma - 1)
    solution = np.append(solution, [timeLayer], axis=0)

    t += tau
    time = np.append(time, [t])

# np.save('t.npy', time)
# np.save('x.npy', x)
# np.save('solution.npy', solution)
