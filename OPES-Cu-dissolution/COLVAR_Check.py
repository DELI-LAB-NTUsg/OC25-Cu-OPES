import numpy as np

file = "COLVAR_PW1"
suptitle = "With CO, Init: solvated, Cu-Cu CN Biased Only"
filename = "cn_analysis.png"

with open(file) as f:
    header = f.readline().split()[2:]  # drop "#!" and "FIELDS"
    rows = [line.split() for line in f if line.strip() and not line.startswith("#")]

# keep only rows with the expected number of columns
ncols = len(header)
good_rows = [r for r in rows if len(r) == ncols]
if len(good_rows) < len(rows):
    print(f"dropped {len(rows) - len(good_rows)} malformed rows out of {len(rows)}")

data = np.array(good_rows, dtype=float)
col = {name: data[:, i] for i, name in enumerate(header)}

cn = col["cn"]
time = np.arange(len(cn)) * 100 * 0.5  # fs, matches STRIDE=100, timestep=0.5fs

cno = col["cno"]

height = col["height.z"]

print("min cn reached:", cn.min(), "at t =", time[cn.argmin()] / 1000, "ps")
print("max cn reached:", cn.max())
print("final cn:", cn[-1])
print("opes.bias final:", col["opes.bias"][-1])

n = len(col["opes.bias"])
print("bias, first third max:", col["opes.bias"][:n//3].max())
print("bias, middle third max:", col["opes.bias"][n//3:2*n//3].max())
print("bias, last third max:", col["opes.bias"][2*n//3:].max())

print("\nneff, start:", col["opes.neff"][0], "end:", col["opes.neff"][-1])
print("zed, start:", col["opes.zed"][0], "end:", col["opes.zed"][-1])

import matplotlib.pyplot as plt

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 11))

# --- top: histogram of cn ---
ax1.hist(cn, bins=50, edgecolor="black", alpha=0.7)
ax1.set_xlabel("CN")
ax1.set_ylabel("Count")
ax1.set_title("Histogram of Cu-Cu CN")

# --- middle1: Cu-Cu CN vs time, shaded where cn <= 3 ---
time_ps = time / 1000  # convert fs -> ps

ax2.plot(time_ps, cn, color="black", linewidth=0.8)
ax2.fill_between(time_ps, cn, cn.min(), where=(cn <= 3),
                  color="tab:red", alpha=0.3, interpolate=True,
                  label="CN ≤ 3")
ax2.axhline(3, color="tab:red", linestyle="--", linewidth=0.8)
ax2.set_xlabel("Time (ps)")
ax2.set_ylabel("Cu-Cu CN")
ax2.set_title("Cu-Cu CN")
ax2.legend()

# --- middle2: Cu height vs time
ax3.plot(time_ps, height, color="black", linewidth=0.8)
ax3.set_xlabel("Time (ps)")
ax3.set_ylabel("Cu height (Å)")
ax3.set_title("Cu height")
ax3.legend()

# --- bottom: opes.bias vs time ---
ax4.plot(time_ps, col["opes.bias"], color="tab:blue", linewidth=0.8)
ax4.set_xlabel("Time (ps)")
ax4.set_ylabel("opes.bias (eV)")
ax4.set_title("OPES Bias")

plt.suptitle(suptitle)
plt.tight_layout()
plt.savefig(filename, dpi=150)
plt.show()