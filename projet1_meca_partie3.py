import numpy as np
import matplotlib.pyplot as plt

g, L = 9.81, 1.0
omega_lin = np.sqrt(g/L)
T = 10
dt = 0.01

# Fonction : calcul de la solution linéaire
def theta_lin(t, theta0, omega0):
    return theta0 * np.cos(omega_lin*t) + (omega0/omega_lin) * np.sin(omega_lin*t)

# Fonction : calcul de la solution numérique avec euler
def euler2(theta0, omega0, dt, T):
    N = int(T/dt)
    theta = np.zeros(N)
    omega_t = np.zeros(N)
    theta[0] = theta0
    omega_t[0] = omega0
    
    for i in range(N - 1):
        omega_t[i+1] = omega_t[i] -(g/L) * np.sin(theta[i]) * dt
        theta[i+1] = theta[i]+omega_t[i] * dt
    t = np.linspace(0, T, N)
    return t, theta

# Grilles d'échantillonnage
theta0_values = np.linspace(0, 1.0, 30)# angles initiaux [0,1] rad
v0_values = np.linspace(0, 2.0, 30)# vitesses initiales [0,2] m/s

errors = np.zeros((len(theta0_values), len(v0_values)))
seuil = 0.05 # seuil acceptable 

# Calcul des erreurs RMS
for i, theta0 in enumerate(theta0_values):
    for j, v0 in enumerate(v0_values):
        omega0 = v0/L
        t, theta_num = euler2(theta0, omega0, dt, T)
        theta_lin_vals = theta_lin(t, theta0, omega0)
        errors[i, j] = np.sqrt(np.mean((theta_num-theta_lin_vals) ** 2))

# Détection de la zone acceptable
mask = errors < seuil
i_max, j_max = np.where(mask)
theta0_max = theta0_values[np.max(i_max)]
v0_max = v0_values[np.max(j_max)]

# Travcer figure
plt.figure(figsize=(8,6))
im = plt.imshow(errors, origin='lower',
                extent=[v0_values[0], v0_values[-1],
                        theta0_values[0],
                        theta0_values[-1]],
                aspect='auto', cmap='viridis')
plt.colorbar(im, label='Erreur RMS (rad)')
plt.xlabel('Vitesse initiale $v_0$ (m/s)')
plt.ylabel('Angle initial $\\theta_0$ (rad)')
plt.title('Zone d\'application de l\'approximation des petits angles')

# Tracer de zone admissible
plt.plot([0, v0_max, v0_max, 0, 0], [0, 0, theta0_max, theta0_max, 0], color='r', lw=2, label='Zone acceptable')

plt.legend()
plt.grid(True)
plt.show()

# Affichage des bornes numériques
print(f"Zone acceptable (erreur RMS < {seuil} rad):")
print(f"- θ_0 <= {theta0_max:.3f} rad (~{np.degrees(theta0_max):.1f}°)")
print(f"- v_0 <= {v0_max:.3f} m/s")
