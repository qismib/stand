import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # For non-interactive backend

def main():
    theta = 0.195
    delta_m2 = 20
    p = 0.1

    # Time array
    t = np.linspace(0, 4 * np.pi * p / delta_m2, 20)
    
    # Phase factor
    phase = delta_m2 * t / (4 * p)
    
    # Survival probabilities
    P_e = 1 - (np.sin(2 * theta) ** 2) * (np.sin(phase) ** 2)
    P_x = 1 - P_e  # Complement

    # Create subplots
    fig, axs = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

    # Plot P_e
    axs[0].plot(t, P_e, lw=2, color='blue')
    axs[0].set_ylabel(r"$P_e(t)$", fontsize=14)
    axs[0].set_ylim(0, 1.05)
    axs[0].set_title("Electron Neutrino Survival Probability")

    # Plot P_x
    axs[1].plot(t, P_x, lw=2, color='red')
    axs[1].set_xlabel(r"$t$", fontsize=14)
    axs[1].set_ylabel(r"$P_x(t)$", fontsize=14)
    axs[1].set_ylim(0, 1.05)
    axs[1].set_title("Complement Survival Probability")

    plt.tight_layout()
    fig.savefig(r"C:\Users\edori\Desktop\Nexus\Università\Current\MasterThesis\Images\survival_probability.pdf")

if __name__ == "__main__":
    main()