from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).parent
G = 9.81 # gravitational acceleration [m/s^2]
M = 1.0 # mass of the body [kg]
F_COEF = (1, -1, 0, 0, 0, 0) # f(x,y) = a*x^2 + b*y^2 + c*x*y + p*x + q*y + r
START = (-10, 0.01, 200) # initial position xyz
V_START = (10, 0, 0) # initial velocity v
K_LIST = (1.0, 0.9) # energy losses
N_BOUNCES = 5
EPS = 1e-9 # roots smaller than this are the starting point itself (t = 0)


def surface(x, y):
    a, b, c, p, q, r = F_COEF
    return a * x**2 + b * y**2 + c * x * y + p * x + q * y + r # universal for any quadratic surface


def surface_gradient(x, y):
    # partial derivatives df/dx and df/dy
    a, b, c, p, q, r = F_COEF
    return 2 * a * x + c * y + p, 2 * b * y + c * x + q


def impact_time(pos, vel):
    # z(t) = f(x(t), y(t)) is a quadratic equation A*t^2 + B*t + C = 0 for f of degree <= 2
    a, b, c, p, q, r = F_COEF
    x0, y0, z0 = pos
    vx, vy, vz = vel
    A = -G / 2 - (a * vx**2 + b * vy**2 + c * vx * vy)
    B = vz - (2 * a * x0 * vx + 2 * b * y0 * vy + c * (x0 * vy + y0 * vx) + p * vx + q * vy)
    C = z0 - surface(x0, y0)

    roots = np.roots([A, B, C])
    roots = roots[np.abs(roots.imag) < 1e-9].real # reject complex solutions
    roots = roots[roots > EPS] # remove negative solutions and zeros
    if len(roots) == 0:
        return None
    return roots.min()


def simulate(k):
    # returns a table of impacts (x, y, z, dt, Ekin, Epot, Etotal) and the paths between them
    pos = np.array(START, dtype=float)
    vel = np.array(V_START, dtype=float)
    rows, paths = [], []

    for _ in range(N_BOUNCES):
        t_i = impact_time(pos, vel) # check the time of possible impact
        if t_i is None:
            print("no collision")
            break

        # build trajecrory
        t = np.linspace(0, t_i, 100)
        # equation of motion: x(t) = x0 + vx*t, y(t) = y0 + vy*t, z(t) = z0 + vz*t - 0.5*G*t^2
        paths.append(np.array([pos[0] + vel[0] * t,
                               pos[1] + vel[1] * t,
                               pos[2] + vel[2] * t - 0.5 * G * t**2]))
        pos = paths[-1][:, -1]
        v_before = np.array([vel[0], vel[1], vel[2] - G * t_i])

        # normal vector to the surface at the impact
        fx, fy = surface_gradient(pos[0], pos[1])
        n = np.array([-fx, -fy, 1.0])
        n /= np.linalg.norm(n)

        e_kin = 0.5 * M * v_before @ v_before
        e_pot = M * G * pos[2]
        rows.append([*pos, t_i, e_kin, e_pot, e_kin + e_pot])

        vel = np.sqrt(k) * (v_before - 2 * (v_before @ n) * n) # reflected velocity

    return np.array(rows), paths


for k in K_LIST:
    table, paths = simulate(k)

    print(f"\nk = {k}") # table
    print(f"{'x_I':>10}{'y_I':>10}{'z_I':>10}{'dt':>10}{'E_kin':>12}{'E_pot':>12}{'E_total':>12}")
    for row in table:
        print("".join(f"{v:>{w}.4f}" for v, w in zip(row, (10, 10, 10, 10, 12, 12, 12))))

    # surface and trajectories
    lim = 1.2 * max(np.abs(path[:2]).max() for path in paths)
    X, Y = np.meshgrid(np.linspace(-lim, lim, 40), np.linspace(-lim, lim, 40))

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection="3d")
    ax.plot_wireframe(X, Y, surface(X, Y), color="steelblue", linewidth=0.4, alpha=0.6)
    for path in paths:
        ax.plot(*path, color="red", linewidth=2)
    ax.scatter(*table[:, :3].T, color="black", s=15)  # impact points
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(f"First {len(paths)} bounces, k = {k}")
    plt.savefig(OUT / f"fig_bounces_k{k}.png", dpi=200)

plt.show()
