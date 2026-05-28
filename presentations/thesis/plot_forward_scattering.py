import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def draw_diagram(ax, left_up, left_down, right_up, right_down, fs):
    line_color = "black"
    box_style = dict(boxstyle="round,pad=0.4", fc="white", ec="none")
    mom_style = dict(boxstyle="round,pad=0.3", fc="white", ec="none")

    segments = [
        (-1,  1, 0,  0, r"$\vec p$", ( 0.15,  0.15)),
        (-1, -1, 0,  0, r"$\vec q$", (-0.20,  0.15)),
        ( 0,  0, 1,  1, r"$\vec p$", (-0.20,  0.15)),
        ( 0,  0, 1, -1, r"$\vec q$", ( 0.15,  0.15)),
    ]

    for x0, y0, x1, y1, mom, (ox, oy) in segments:
        ax.plot([x0, x1], [y0, y1], color=line_color, lw=2)
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = (x1 - x0) * 0.12, (y1 - y0) * 0.12
        ax.arrow(mx - dx / 2, my - dy / 2, dx, dy,
                 head_width=0.08, head_length=0.12,
                 fc=line_color, ec=line_color,
                 length_includes_head=True)
        ax.text(mx + ox, my + oy, mom,
                ha="center", va="center",
                fontsize=fs, bbox=mom_style)

    ax.text(-1.25,  1.25, left_up,   ha="center", va="center", fontsize=fs, bbox=box_style)
    ax.text(-1.25, -1.25, left_down, ha="center", va="center", fontsize=fs, bbox=box_style)
    ax.text( 1.25,  1.25, right_up,  ha="center", va="center", fontsize=fs, bbox=box_style)
    ax.text( 1.25, -1.25, right_down,ha="center", va="center", fontsize=fs, bbox=box_style)

    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.7, 1.7)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    fs = 15

    draw_diagram(
        axes[0],
        r"$\nu_e$", r"$\nu_x$",
        r"$\nu_x$", r"$\nu_e$",
        fs
    )

    draw_diagram(
        axes[1],
        r"$\nu_x$", r"$\nu_e$",
        r"$\nu_e$", r"$\nu_x$",
        fs
    )

    plt.tight_layout()
    fig.savefig(r"C:\Users\edori\Desktop\Nexus\Università\Current\MasterThesis\Images\forward_scattering_diagrams\flavor_change_diagrams.pdf")


if __name__ == "__main__":
    main()