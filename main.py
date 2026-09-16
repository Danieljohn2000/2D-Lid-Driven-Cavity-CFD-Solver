import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. GRID
# ============================================================

nx = 51
ny = 51

L = 1.0

dx = L / (nx - 1)
dy = L / (ny - 1)


# ============================================================
# 2. PHYSICAL PARAMETERS
# ============================================================

U_lid = 1.0

Re = 100

nu = U_lid * L / Re

dt = 0.001

num_steps = 2000


# ============================================================
# 3. CREATE FLOW VARIABLES
# ============================================================

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))

omega = np.zeros((ny, nx))
psi = np.zeros((ny, nx))


# ============================================================
# 4. MAIN TIME LOOP
# ============================================================

for step in range(num_steps):

    # --------------------------------------------------------
    # Solve streamfunction Poisson equation
    # --------------------------------------------------------

    for iteration in range(100):

        psi_new = psi.copy()

        psi_new[1:-1, 1:-1] = 0.25 * (
            psi[1:-1, 2:]
            + psi[1:-1, :-2]
            + psi[2:, 1:-1]
            + psi[:-2, 1:-1]
            + dx**2 * omega[1:-1, 1:-1]
        )

        psi = psi_new

    # --------------------------------------------------------
    # Calculate velocity from streamfunction
    # --------------------------------------------------------

    u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2 * dy)

    v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, :-2]) / (2 * dx)

    # --------------------------------------------------------
    # Boundary conditions for velocity
    # --------------------------------------------------------

    # Bottom wall
    u[0, :] = 0
    v[0, :] = 0

    # Left wall
    u[:, 0] = 0
    v[:, 0] = 0

    # Right wall
    u[:, -1] = 0
    v[:, -1] = 0

    # Moving top wall
    u[-1, :] = U_lid
    v[-1, :] = 0

    # --------------------------------------------------------
    # Update vorticity inside the fluid
    # --------------------------------------------------------

    d_omega_dx = (omega[1:-1, 2:] - omega[1:-1, :-2]) / (2 * dx)

    d_omega_dy = (omega[2:, 1:-1] - omega[:-2, 1:-1]) / (2 * dy)

    d2_omega_dx2 = (omega[1:-1, 2:] - 2 * omega[1:-1, 1:-1] + omega[1:-1, :-2]) / dx**2

    d2_omega_dy2 = (omega[2:, 1:-1] - 2 * omega[1:-1, 1:-1] + omega[:-2, 1:-1]) / dy**2

    # --------------------------------------------------------
    # Navier-Stokes vorticity equation
    # --------------------------------------------------------

    omega_new = omega.copy()

    omega_new[1:-1, 1:-1] = (
        omega[1:-1, 1:-1]
        - dt * (u[1:-1, 1:-1] * d_omega_dx + v[1:-1, 1:-1] * d_omega_dy)
        + nu * dt * (d2_omega_dx2 + d2_omega_dy2)
    )

    # --------------------------------------------------------
    # Vorticity boundary conditions
    # --------------------------------------------------------

    # Bottom wall
    omega_new[0, 1:-1] = -2 * psi[1, 1:-1] / dy**2

    # Top moving wall
    omega_new[-1, 1:-1] = -2 * psi[-2, 1:-1] / dy**2 - 2 * U_lid / dy

    # Left wall
    omega_new[1:-1, 0] = -2 * psi[1:-1, 1] / dx**2

    # Right wall
    omega_new[1:-1, -1] = -2 * psi[1:-1, -2] / dx**2

    omega = omega_new

    # --------------------------------------------------------
    # Display progress
    # --------------------------------------------------------

    if step % 200 == 0:
        print(f"Step {step} / {num_steps}")


# ============================================================
# 5. FINAL VELOCITY CALCULATION
# ============================================================

u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2 * dy)

v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, :-2]) / (2 * dx)


# ============================================================
# 6. CREATE VELOCITY MAGNITUDE
# ============================================================

velocity_magnitude = np.sqrt(u**2 + v**2)


# ============================================================
# 7. CREATE COORDINATES
# ============================================================

x = np.linspace(0, L, nx)
y = np.linspace(0, L, ny)

X, Y = np.meshgrid(x, y)


# ============================================================
# 8. PLOT VELOCITY MAGNITUDE
# ============================================================

plt.figure(figsize=(8, 6))

plt.contourf(X, Y, velocity_magnitude, levels=30)

plt.colorbar(label="Velocity Magnitude")

plt.xlabel("x position")
plt.ylabel("y position")

plt.title("2D Lid-Driven Cavity Flow - Velocity Magnitude")

plt.tight_layout()

plt.savefig("velocity_magnitude.png", dpi=300)

plt.show()


# ============================================================
# 9. PLOT VELOCITY VECTORS
# ============================================================

plt.figure(figsize=(8, 6))

plt.contourf(X, Y, velocity_magnitude, levels=30)

plt.colorbar(label="Velocity Magnitude")

plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])

plt.xlabel("x position")
plt.ylabel("y position")

plt.title("2D Lid-Driven Cavity Flow - Velocity Field")

plt.tight_layout()

plt.savefig("velocity_field.png", dpi=300)

plt.show()


print("Simulation complete.")
