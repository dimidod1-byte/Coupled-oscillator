#coupled osccillators
#coupled oscillators simulation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
def derivs(state, t):
    dydx = np.zeros_like(state)
    dydx[0] = state[1]
    
    # Forces based on Hooke's Law for 3 springs: F = -k * x
    f1 = -k1 * state[0] + k2 * (state[2] - state[0])
    f2 = -k2 * (state[2] - state[0]) - k3 * state[2]
    
    dydx[1] = f1 / m1
    dydx[2] = state[3]
    dydx[3] = f2 / m2
    return dydx 
def runge_kutta_4th_order(state, t, dt):
    k1 = derivs(state, t)
    k2 = derivs(state + 0.5 * dt * k1, t + 0.5 * dt)
    k3 = derivs(state + 0.5 * dt * k2, t + 0.5 * dt)
    k4 = derivs(state + dt * k3, t + dt)
    return state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
def simulate_coupled_oscillators(initial_state, t_max, dt):
    t = np.arange(0, t_max, dt)
    state = np.zeros((len(t), len(initial_state)))
    state[0] = initial_state
    for i in range(1, len(t)):
        state[i] = runge_kutta_4th_order(state[i-1], t[i-1], dt)
    return t, state
def plot_coupled_oscillators(t, state):
    # Equilibrium positions of mass 1 and mass 2 along the X axis
    eq1, eq2 = 2.0, 4.0
    x1 = eq1 + state[:, 0]
    x2 = eq2 + state[:, 2]
    y1 = np.zeros_like(x1)
    y2 = np.zeros_like(x2)
    
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_xlim(0, 6.0)
    ax.set_ylim(-1.0, 1.0)
    ax.grid()
    ax.set_title('Coupled Oscillators Simulation')
    ax.set_xlabel('X Position')
    ax.get_yaxis().set_visible(False)

    # Plot walls, masses, and springs
    ax.plot([0, 0], [-0.5, 0.5], color='black', lw=3)  # Left wall
    ax.plot([6, 6], [-0.5, 0.5], color='black', lw=3)  # Right wall
    
    springs, = ax.plot([], [], '-', lw=1, color='gray')
    masses, = ax.plot([], [], 's', markersize=15, color='blue')
    
    def init():
        springs.set_data([], [])
        masses.set_data([], [])
        return springs, masses

    def animate(i):
        # Draw lines representing the continuous connection across springs
        spring_x = [0, x1[i], x2[i], 6.0]
        spring_y = [0, 0, 0, 0]
        springs.set_data(spring_x, spring_y)
        
        mass_x = [x1[i], x2[i]]
        mass_y = [0, 0]
        masses.set_data(mass_x, mass_y)
        return springs,