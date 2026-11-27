import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("datos.txt")

M = 1.0
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

titulos = {
    1: "Caso 1: Fotón cae",
    2: "Caso 2: Órbita inestable",
    3: "Caso 3: Fotón escapa"
}
colores = {1: "crimson", 2: "darkorange", 3: "royalblue"}

for idx, caso in enumerate([1, 2, 3]):
    ax = axes[idx]
    sub = df[df["caso"] == caso]
    b_val = sub["b"].iloc[0]

    # Horizonte de sucesos (r = 2M) y Esfera de fotones (r = 3M)
    ax.add_patch(plt.Circle((0, 0), 2 * M, color="black", zorder=5, label="Horizonte (2M)"))
    ax.add_patch(plt.Circle((0, 0), 3 * M, color="gray", fill=False, linestyle="--", zorder=4, label="Esfera fotones (3M)"))

    # Trayectoria y punto de partida
    ax.plot(sub["x"], sub["y"], color=colores[caso], linewidth=2, label=f"b = {b_val:.4f}")
    ax.scatter(sub["x"].iloc[0], sub["y"].iloc[0], color="green", zorder=6, label="Inicio")

    ax.set_title(f"{titulos[caso]}\n$b = {b_val:.4f}$", fontweight="bold")
    ax.set_xlabel("x / M")
    ax.set_ylabel("y / M")
    ax.set_xlim(-12, 12)
    ax.set_ylim(-12, 12)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)

plt.tight_layout()
plt.savefig('trayectorias_3casos.png')