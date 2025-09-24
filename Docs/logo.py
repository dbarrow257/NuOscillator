import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap

# --- Create figure and axis ---
fig, ax = plt.subplots(figsize=(12, 12))
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# --- Simple Earth layers ---
scale = 0.9
layers = [
    {"radius": 1.0 * scale, "color": "#4682B4", "label": "Atmosphere"},
    {"radius": 0.95 * scale, "color": "#1E3A8A", "label": "Outer Layer"},
    {"radius": 0.7 * scale, "color": "#FFB347", "label": "Mantle"},
    {"radius": 0.3 * scale, "color": "#800000", "label": "Core"}
]

# --- Draw layers ---
for layer in layers:
    radius = layer["radius"] / 2
    ax.add_patch(Circle((0.5, 0.5), radius, color=layer["color"], zorder=1))

def draw_oscillation_system(ax):
    x = np.linspace(0, 1, 1000)
    mantle_radius = layers[2]["radius"] / 2
    y_center = 0.5
    x_offset_base = 0.15
    amplitude_scale = mantle_radius - layers[3]["radius"]/2

    # --- Raw amplitude curves ---
    blue_raw  = (1 - x)**4       # blue drops faster
    red_raw   = x**1.5           # red rises faster initially
    green_raw = 0.3 * x          # green slower increase

    # --- Normalize each curve 0->1 ---
    def normalize_0_to_1(arr):
        return (arr - arr.min()) / (arr.max() - arr.min())

    blue_norm  = normalize_0_to_1(blue_raw)
    red_norm   = normalize_0_to_1(red_raw)
    green_norm = normalize_0_to_1(green_raw)

    # --- Normalize so sum = 1 ---
    total_norm = blue_norm + red_norm + green_norm
    blue_prob  = blue_norm / total_norm
    red_prob   = red_norm / total_norm
    green_prob = green_norm / total_norm

    # --- Waves configuration ---
    waves = [
        {"freq": 3, "base_amp": 0.4, "color": "blue",  "lw": 4, "offset": 0.03, "prob": blue_prob,  "phase_shift": 0.25},
        {"freq": 6, "base_amp": 0.35,"color": "red",   "lw": 3, "offset": 0.02, "prob": red_prob,   "phase_shift": 0},
        {"freq": 14,"base_amp": 0.35,"color": "green", "lw": 2, "offset": 0.01, "prob": green_prob, "phase_shift": 0}
    ]

    h = (x * 0.055) * np.sin(50*x)  # small extra oscillation

    for i, wave in enumerate(waves):
        y = np.sin(wave["freq"] * x * 2 * np.pi + wave["phase_shift"]) \
            * wave["base_amp"] * wave["prob"] * amplitude_scale

        x_shifted = x + x_offset_base * (1.2 - 0.1*i)
        mask = ((x_shifted - 0.5)**2 + (y + y_center - y_center)**2 <= mantle_radius**2)

        lw_scaled = wave["lw"] * (1.3 + 0.7 * np.mean(wave["prob"][mask]))  # scales lw by 1→1.5

        ax.plot(x_shifted[mask], y[mask] + y_center + h[mask] * wave["offset"] * i,
                color=wave["color"], lw=lw_scaled, zorder=3+i)

        indices = np.arange(0, len(x_shifted[mask]), 50)
        ax.scatter(x_shifted[mask][indices], y[mask][indices] + y_center,
                   color='magenta', s=10, alpha=0.6, zorder=6)

draw_oscillation_system(ax)

# --- Arrow ---
y_arrow = 0.52
arrow = FancyArrowPatch(
    posA=(0.5 - layers[0]["radius"]/2, y_arrow),
    posB=(0.5 - layers[2]["radius"]/2, y_arrow),
    arrowstyle='->,head_width=0.12,head_length=0.18',  # bigger arrowhead
    color='black',  # keep black
    lw=8,  # thicker line
    mutation_scale=30,  # bigger head
    zorder=4
)
ax.add_patch(arrow)

# --- Add "ν" label above arrow ---
ax.text(
    0.5 - layers[0]["radius"]/2 + 0.08, y_arrow + 0.03, r"$\nu$",
    fontsize=40, ha="center", va="bottom", weight="bold", color="black", zorder=6,
    path_effects=[path_effects.withStroke(linewidth=4, foreground="white")]  # makes it stand out
)

# --- Gradient-Colored Text ---
text_str = "NuOscillator"
colors = ["red", "green", "blue"]
cmap = LinearSegmentedColormap.from_list("osc_colors", colors)

renderer = fig.canvas.get_renderer()
total_width = 0
letter_extents = []
for letter in text_str:
    t = ax.text(0, 0, letter, fontsize=120, weight="bold", family="DejaVu Sans")
    bb = t.get_window_extent(renderer=renderer)
    width_data = (bb.width / fig.dpi) / fig.get_size_inches()[0]
    letter_extents.append(width_data)
    total_width += width_data
    t.remove()

x_pos = 0.5 - total_width / 2
for i, letter in enumerate(text_str):
    color = cmap(i / (len(text_str)-1))
    txt = ax.text(x_pos, 0.15, letter, fontsize=120, color=color,
                  ha="left", va="center", weight="bold", family="DejaVu Sans", zorder=5)
    txt.set_path_effects([path_effects.withStroke(linewidth=4, foreground="white")])
    x_pos += letter_extents[i]

plt.savefig("NuOscillatorLogo.png", bbox_inches="tight", facecolor="white")
