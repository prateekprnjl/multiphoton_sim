#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load_data(filename):
    """Load cross-section data from the simulation output."""
    data = np.loadtxt(filename)

    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError("Invalid Cross-section observable file.")

    return data


def plot_cross_sections(filename, states=None, all=False, output=None, show=True):
    """
    Plot cation cross sections as a function of energy.

    Parameters
    ----------
    filename : str or Path
    states : Cation states to plot. If None, all cation states are plotted.
    output : Filename for saving the figure.
    show : Whether to display the figure.
    """

    data = load_data(filename)
    # First column is energy.
    energy = data[:, 0]
    # Number of cation states.
    n_states = data.shape[1] - 1

    # If no states were specified, plot total cross section only.
    if all:
        states = list(range(1, n_states + 1))
    elif states is None:
        states = [n_states]

    # Check requested states.
    for state in states:
        if state < 1 or state > n_states:
            raise ValueError(f"Cation state {state} does not exist. "
                f"The output contains {n_states} cation states.")

    fig, ax = plt.subplots(figsize=(8, 5))

    for state in states:
        cross_section = data[:, state]

        if state == n_states:
            ax.plot(energy, cross_section, linewidth=1.5, label=f"Total")
        else:
            ax.plot(energy, cross_section, linewidth=1.5, label=f"Cation {state}")

    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(r"Cross section (Mb)")
    ax.legend()
    ax.grid(alpha=0.25)

    fig.tight_layout()

    if output is not None:
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print(f"Saved figure to {output}")

    if show:
        plt.show()

    plt.close(fig)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--states", nargs="+", type=int)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--output", "-o")

    args = parser.parse_args()

    plot_cross_sections(args.input, args.states, args.all, args.output)


if __name__ == "__main__":
    main()
