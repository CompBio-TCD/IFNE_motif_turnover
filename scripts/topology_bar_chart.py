# Creating barchart results from topology test 
import numpy as np
import matplotlib.pyplot as plt

# AU p-values for each test 
hypotheses = ["IFNE + IFNA", "IFNE + IFNB",
              "IFNE + IFNK", "IFNE + IFNW"]
p_values = [4.72e-050, 5.25e-006, 0.307, 0.0053]

x = np.arange(len(hypotheses))
colours = ["#4C78A8" if p >= 0.05 else "#E69F00"
           for p in p_values]

fig, (ax_top, ax_bottom) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 5),
    gridspec_kw={"height_ratios": [1, 2], "hspace": 0.06}
)

for ax in (ax_top, ax_bottom):
    ax.bar(x, p_values, color=colours, width=0.6)
    ax.axhline(0.05, color="red", linestyle="--",
               linewidth=1.4)

# Adjust ranges as needed 
ax_bottom.set_ylim(0, 0.08)
ax_top.set_ylim(0.15, 0.50)

# Hide adjoining borders
ax_top.spines["bottom"].set_visible(False)
ax_bottom.spines["top"].set_visible(False)
ax_top.tick_params(bottom=False, labelbottom=False)

# Diagonal break marks
d = 0.012
kwargs = dict(color="black", clip_on=False, linewidth=1)

for ax, y, sign in [(ax_top, 0, 1), (ax_bottom, 1, 1)]:
    ax.plot((-d, d), (y-d, y+d),
            transform=ax.transAxes, **kwargs)
    ax.plot((1-d, 1+d), (y-d, y+d),
            transform=ax.transAxes, **kwargs)

ax_bottom.set_xticks(x)
ax_bottom.set_xticklabels(hypotheses)
fig.supylabel("AU test p-value", x=0.03)
ax_bottom.text(3.45, 0.052, "p = 0.05",
               color="red", fontsize=9, ha="right")

plt.savefig("IFNE_AU.pdf", bbox_inches="tight")
plt.savefig("IFNE_AU.png", dpi=300,
            bbox_inches="tight")
plt.show()
