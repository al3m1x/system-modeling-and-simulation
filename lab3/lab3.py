from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).parent

MU = 0.01  # gęstość liniowa struny [kg/m]
TENSION = 100.0  # naprężenie struny [N]
L = 1.0  # długość struny [m]
M = 100  # liczba węzłów w przestrzeni
N = 2000  # liczba poziomów czasowych
B_LIST = (0.0, 0.5)  # współczynniki tłumienia do symulacji (0 = brak tłumienia)

S = lambda x: np.sin(3 * np.pi * x + np.pi / 4)  # kształt początkowy y(x, 0)
G = lambda x: x - x**2  # prędkość początkowa dy/dt(x, 0)
LEFT = lambda t: np.sqrt(2) / 2  # lewy koniec: wychylenie y(0, t) (warunek Dirichleta)
RIGHT = lambda t: -3 * np.pi * np.sqrt(2) / 2  # prawy koniec: nachylenie dy/dx(L, t) (warunek Neumanna)


def simulate(b):
    beta = b / (2 * MU)  # tempo wygasania drgań [1/s]; z porównania b/mu = 2*beta
    x = np.linspace(0, L, M)
    dx = x[1] - x[0]  # = L / (M - 1), bo M to liczba węzłów, nie przedziałów
    v = np.sqrt(TENSION / MU)  # prędkość fali w strunie
    dt = dx / v  # warunek stabilności CFL
    t = np.arange(N) * dt

    p = (v * dt / dx) ** 2
    q = 1 + beta * dt
    u = 1 - beta * dt

    y = np.zeros((N, M))
    y[0] = S(x)

    # z wyprowadzenia
    y[1, 1:-1] = p / 2 * (y[0, 2:] - 2 * y[0, 1:-1] + y[0, :-2]) + y[0, 1:-1] + u * dt * G(x[1:-1])
    y[1, 0] = LEFT(t[1])  # warunek Dirichleta: wartość zadana wprost
    y[1, -1] = p * (y[0, -2] - y[0, -1] + dx * RIGHT(t[0])) + y[0, -1] + u * dt * G(x[-1])

    for k in range(1, N - 1):
        # wzór ogólny (wykład 5): nowy punkt z jego dwóch poprzednich chwil i sąsiadów
        y[k + 1, 1:-1] = (p / q * (y[k, 2:] - 2 * y[k, 1:-1] + y[k, :-2])
                          + 2 / q * y[k, 1:-1] - u / q * y[k - 1, 1:-1])
        y[k + 1, 0] = LEFT(t[k + 1])
        # prawy koniec (Neumann)
        y[k + 1, -1] = 2 * p / q * (y[k, -2] - y[k, -1] + dx * RIGHT(t[k])) + 2 / q * y[k, -1] - u / q * y[k - 1, -1]

    return x, t, y

# stan równowagi (nic się już nie zmienia w czasie): y_t = y_tt = 0 => y_xx = 0,
# czyli y_E = c*x + d; c z nachylenia na prawym końcu, d z wychylenia na lewym (y_E(0) = d)
c = RIGHT(0)
d = LEFT(0)
print(f"equilibrium solution: y_E(x) = {c:g}*x {d:+g}")

for b in B_LIST:
    x, t, y = simulate(b)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection="3d")
    X, T = np.meshgrid(x, t)
    ax.plot_surface(X, T, y, cmap="viridis", rstride=5, cstride=1, linewidth=0)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("t [s]")
    ax.set_zlabel("y [m]")
    ax.set_title(f"String vibrations, b = {b}")
    plt.savefig(OUT / f"fig_surface_b{b}.png", dpi=200)

    if b > 0:
        # porównanie ostatniej iteracji ze stanem równowagi
        plt.figure()
        plt.plot(x, y[-1], linewidth=3, label="last iteration")
        plt.plot(x, c * x + d, "r", linewidth=1, label="equilibrium")
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.title(f"Equilibrium solution, b = {b}")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.savefig(OUT / f"fig_equilibrium_b{b}.png", dpi=200)
        print(f"b = {b}: max |last iteration - equilibrium| = {np.abs(y[-1] - (c * x + d)).max():.4f}")

plt.show()
