import numpy as np
import stylia
import matplotlib.pyplot as plt

import sys
import os

ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "..", "tools", "blitzgsea"))

COLORS = stylia.colors.colors.NamedColors()

# utils


def truncate_string(value, max_length=50, suffix="..."):
    string_value = str(value)
    string_truncated = string_value[
        : min(len(string_value), (max_length - len(suffix)))
    ]
    suffix = suffix if len(string_value) > max_length else ""
    return string_truncated + suffix


# plots


def contingency_table_plot(term, a, b, c, d):
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    n = a + b + c + d
    factor = 30000
    a_ = a / n
    b_ = b / n
    c_ = c / n
    d_ = d / n
    ax.scatter([0], [0], s=np.sqrt(a_) * factor, color=COLORS.red, alpha=0.5)
    ax.text(0, 0, a, ha="center", va="center", fontsize=16)
    ax.scatter([0], [1], s=np.sqrt(b_) * factor, color=COLORS.blue, alpha=0.5)
    ax.text(0, 1, b, ha="center", va="center", fontsize=16)
    ax.scatter([1], [0], s=np.sqrt(c_) * factor, color=COLORS.blue, alpha=0.5)
    ax.text(1, 0, c, ha="center", va="center", fontsize=16)
    ax.scatter([1], [1], s=np.sqrt(d_) * factor, color=COLORS.gray, alpha=0.5)
    ax.text(1, 1, d, ha="center", va="center", fontsize=16)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Yes", "No"], fontsize=16)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Yes", "No"], fontsize=16)
    ax.set_title(truncate_string(term), fontsize=16)
    ax.set_xlabel("In profile", fontsize=16)
    ax.set_ylabel("In category", fontsize=16)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(1.5, -0.5)
    return fig


def projection_plot(selected_proteins, all_proteins, X):
    selected_proteins = set(selected_proteins)
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.scatter(X[:, 0], X[:, 1], color=COLORS.gray, s=5)
    idxs = []
    for i, pid in enumerate(all_proteins):
        if pid in selected_proteins:
            idxs += [i]
    idxs = np.array(idxs)
    ax.scatter(
        X[idxs, 0], X[idxs, 1], color=COLORS.red, edgecolors="white", s=100, alpha=1
    )
    ax.set_axis_off()
    return fig
