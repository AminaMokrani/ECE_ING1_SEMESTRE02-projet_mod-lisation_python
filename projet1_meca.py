import numpy as np
import matplotlib.pyplot as plt

plt.close('all')  # Ferme toutes les figures précédentes

# Fonction : méthode d'Euler pour pendule amorti 
def euler2(theta0, omega0, l, m, lambd, g, t_min, t_max, dt):
    N = int((t_max - t_min) / dt)
    t = np.linspace(t_min, t_max, N)
    theta = np.zeros(N)
    omega = np.zeros(N)
    theta[0] = theta0
    omega[0] = omega0

    for i in range(N - 1):
        omega[i + 1] = omega[i] + dt * (-lambd / m * omega[i] - (g / l) * np.sin(theta[i]))
        theta[i + 1] = theta[i] + dt * omega[i]
    
    return t, theta 

# Fonction : solution analytique HPA
def theta_analytique(t, theta0, omega0, g, l):
    omega_nat = np.sqrt(g / l) 
    return theta0 * np.cos(omega_nat * t) + (omega0 / omega_nat) * np.sin(omega_nat * t)

# Paramètres de base
g = 9.81
l = 0.1
m = 1
lambda_vals = [2, 20, 50]
theta0 = np.pi / 6
omega0 = 1
t_min = 0
t_max = 10
dt = 0.01

# FIGURE 1 : Comparaison des régimes
fig1 = plt.figure(figsize=(10, 5))
for lambd in lambda_vals:
    t, theta_num = euler2(theta0, omega0, l, m, lambd, g, t_min, t_max, dt)
    theta_an = theta_analytique(t, theta0, omega0, g, l)
    plt.plot(t, theta_num, label=f"Numérique λ={lambd}")
    if lambd == lambda_vals[0]:
        plt.plot(t, theta_an, '--', label="Analytique (HPA)")
plt.title("Figure 1 : Comparaison des régimes d'amortissement")
plt.xlabel("Temps (s)")
plt.ylabel("Angle θ (rad)")
plt.grid()
plt.legend()
plt.tight_layout()

# FIGURE 2 : Influence de la masse
fig2 = plt.figure(figsize=(10, 5))
for m_test in [0.5, 1, 2]:
    t, theta_num = euler2(theta0, omega0, l, m_test, 20, g, t_min, t_max, dt)
    plt.plot(t, theta_num, label=f"m = {m_test} kg")
plt.title("Figure 2 : Influence de la masse")
plt.xlabel("Temps (s)")
plt.ylabel("Angle θ (rad)")
plt.grid()
plt.legend()
plt.tight_layout()

#FIGURE 3 : Influence de la condition initiale 
fig3 = plt.figure(figsize=(10, 5))
for theta_init in [np.pi/12, np.pi/6, np.pi/3]:
    t, theta_num = euler2(theta_init, omega0, l, m, 20, g, t_min, t_max, dt)
    plt.plot(t, theta_num, label=f"θ_0 = {round(theta_init, 2)} rad")
plt.title("Figure 3 : Influence de la position initiale θ_0")
plt.xlabel("Temps (s)")
plt.ylabel("Angle θ (rad)")
plt.grid()
plt.legend()
plt.tight_layout()

#FIGURE 4 : Influence du pas de temps dt
fig4 = plt.figure(figsize=(10, 5))
for dt_test in [0.05, 0.01, 0.001]:
    t, theta_num = euler2(theta0, omega0, l, m, 20, g, t_min, t_max, dt_test)
    plt.plot(t, theta_num, label=f"dt = {dt_test}")
plt.title("Figure 4 : Influence du pas de temps")
plt.xlabel("Temps (s)")
plt.ylabel("Angle θ (rad)")
plt.grid()
plt.legend()
plt.tight_layout()


#FIGURE 5 : Influence de la longueur du fil
fig5 = plt.figure(figsize=(10, 5))
for l_test in [0.05, 0.1, 0.2]:
    t, theta_num = euler2(theta0, omega0, l_test, m, 20, g, t_min, t_max, dt)
    plt.plot(t, theta_num, label=f"l = {l_test} m")
plt.title("Figure 5 : Influence de la longueur du fil")
plt.xlabel("Temps (s)")
plt.ylabel("Angle θ (rad)")
plt.grid()
plt.legend()
plt.tight_layout()


# Afficher toutes les figures
plt.show()
