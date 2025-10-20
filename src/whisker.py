import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


df = pd.read_csv(r"figs\whiskerbox\areas.csv")          # or: pd.read_feather / read_parquet / read_json
# sns.boxplot(data=df, x="Claudin", y="Area")
# plt.yscale("log")   # optional
# plt.show()

# --- style ---
sns.set_theme(style="whitegrid", context="talk")
palette = {"C2":"#4C78A8","C4":"#F58518","C5":"#54A24B","C15":"#B279A2"}

# df must have columns: ["Claudin", "Area"] (and optionally "seed", "is_outer")
# 1) (Optional) drop exterior faces if you carried a flag over from Julia
# df = df[~df["is_outer"]]

# 2) (Optional) trim extreme tails to compare cores fairly (5–95% whiskers)
upper = df["Area"].quantile(0.99)   # global cap; or do per-group transform
df["Area_clip"] = df["Area"].clip(upper=upper)
df["Area"] = df["Area"]*1e-6
fig, ax = plt.subplots(figsize=(12, 6), dpi=200)
sns.boxplot(
    data=df, x="Claudin", y="Area_clip", palette=palette,
    whis=(5, 95), showfliers=False, width=0.5, ax=ax
)

ax.set_yscale("log")
ax.set_xlabel("Claudin")
ax.set_ylabel("Area(µm²)")
ax.set_title("Mesh-area spread by Claudin (outliers hidden, 5–95% whiskers)")

# nicer log ticks
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(mticker.NullFormatter())
ax.ticklabel_format(style="plain", axis="y")

plt.tight_layout()
plt.show()

