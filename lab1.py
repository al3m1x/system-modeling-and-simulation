import numpy as np
import matplotlib.pyplot as plt

# parameters - change only this section
G = 10.0  # gravitational acceleration [m/s^2]
V0_RANGE = (1, 15)  # initial speed range [m/s]
ALPHA_RANGE = (1, 89)  # launch angle range [deg]
N_TRIALS = 100_000
N_PLOTTED = 5  # number of valid trajectories to draw

# regions as lists of segments (x1, y1, x2, y2)
ALLOWED = {
    "R1": [(2, 1, 2, 4)],
    "R3": [(4, 1, 7, 1)],
}
FORBIDDEN = {
    "R2": [(2.8, 1.5, 4, 1.5), (4, 1.5, 4, 3), (4, 3, 2.8, 3), (2.8, 3, 2.8, 1.5)],  # rectangle
}


def trajectory(v0, alpha, n_points=300):
    # positions from launch until the body falls back to y = 0
    a = np.radians(alpha)
    t = np.linspace(0, 2 * v0 * np.sin(a) / G, n_points)
    return v0 * np.cos(a) * t, v0 * np.sin(a) * t - 0.5 * G * t**2


def hit_at_time(t, v0, alpha, seg):
    # check if the body is within the flight and on the segment at time t
    x1, y1, x2, y2 = seg
    a = np.radians(alpha)
    vx, vy = v0 * np.cos(a), v0 * np.sin(a)
    x, y = vx * t, vy * t - 0.5 * G * t**2
    s = ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)) / ((x2 - x1)**2 + (y2 - y1)**2)
    return (t >= 0) & (t <= 2 * vy / G) & (s >= 0) & (s <= 1)


def crosses_segment(v0, alpha, seg):
    # substitute the segment's line into x(t), y(t) -> quadratic in t
    x1, y1, x2, y2 = seg
    a = np.radians(alpha)
    vx, vy = v0 * np.cos(a), v0 * np.sin(a)
    dx, dy = x2 - x1, y2 - y1
    A = -dx * G / 2
    B = -dy * vx + dx * vy
    C = dy * x1 - dx * y1

    if abs(A) < 1e-12:  # horizontal segment -> linear equation
        with np.errstate(divide="ignore", invalid="ignore"):
            return hit_at_time(-C / B, v0, alpha, seg)
    delta = B**2 - 4 * A * C
    root = np.sqrt(np.maximum(delta, 0))
    t1, t2 = (-B - root) / (2 * A), (-B + root) / (2 * A)
    return (delta >= 0) & (hit_at_time(t1, v0, alpha, seg) | hit_at_time(t2, v0, alpha, seg))


def crosses_region(v0, alpha, segments):
    return np.any([crosses_segment(v0, alpha, s) for s in segments], axis=0)


# 1. Monte Carlo search for valid initial values
v0 = np.random.uniform(*V0_RANGE, N_TRIALS)
alpha = np.random.uniform(*ALPHA_RANGE, N_TRIALS)

valid = np.ones(N_TRIALS, dtype=bool)
for segments in ALLOWED.values():
    valid &= crosses_region(v0, alpha, segments)
for segments in FORBIDDEN.values():
    valid &= ~crosses_region(v0, alpha, segments)

plt.figure()
plt.plot(alpha[valid], v0[valid], "k.", markersize=1)
plt.plot([ALPHA_RANGE[0], ALPHA_RANGE[1], ALPHA_RANGE[1], ALPHA_RANGE[0], ALPHA_RANGE[0]],
         [V0_RANGE[0], V0_RANGE[0], V0_RANGE[1], V0_RANGE[1], V0_RANGE[0]], "b-", linewidth=2)
plt.xlabel("alpha [deg]")
plt.ylabel("v0 [m/s]")
plt.title("Initial values fulfilling the conditions")
plt.grid(alpha=0.3)
plt.savefig("fig_initial_values.png", dpi=200)

# 2. regions and valid trajectories
plt.figure()
for regions, color in ((ALLOWED, "green"), (FORBIDDEN, "red")):
    for name, segments in regions.items():
        for x1, y1, x2, y2 in segments:
            plt.plot([x1, x2], [y1, y2], color=color, linewidth=3)
        plt.text(np.mean([s[0] for s in segments]) + 0.1, np.mean([s[1] for s in segments]) + 0.1,
                 name, color=color, fontsize=12)

for i in np.random.choice(np.flatnonzero(valid), min(N_PLOTTED, valid.sum()), replace=False):
    plt.plot(*trajectory(v0[i], alpha[i]), label=f"v0={v0[i]:.2f}, alpha={alpha[i]:.1f}")

all_segments = [s for r in (ALLOWED, FORBIDDEN) for segs in r.values() for s in segs]
plt.xlim(0, max(max(s[0], s[2]) for s in all_segments))
plt.ylim(0, 1.2 * max(max(s[1], s[3]) for s in all_segments))
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Allowed (green) and forbidden (red) regions")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(alpha=0.3)
plt.savefig("fig_trajectories.png", dpi=200, bbox_inches="tight")

# 3. probability of fulfilling the conditions in a single try
print(f"P = {valid.sum()} / {N_TRIALS} = {valid.mean():.5f}")
plt.show()
