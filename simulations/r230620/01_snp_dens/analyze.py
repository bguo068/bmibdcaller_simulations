from pathlib import Path
import pandas as pd
from subprocess import check_output
import io
import matplotlib.pyplot as plt

model = ["single", "multiple"][1]
print(model)
folder_abbr = 4 if model == "single" else 5
df_lst = []
for p in Path("..").glob(f"res/{folder_abbr}000*/vcf/*.vcf.gz"):
    r = int(p.parents[1].name.split("_")[2].replace("r", ""))
    chrno = int(p.name.replace(".vcf.gz", "").split("_")[-1])

    cmd = f"bcftools +fill-tags {p} -Ou -- -t AC,AN | bcftools query -f '%CHROM\t%POS\t%AN\t%AC\n'"
    txt = check_output(cmd, shell=True, text=True)
    f = io.StringIO(txt)
    df = pd.read_csv(f, sep="\t", names=["Chrom", "Pos", "AN", "AC"])
    df["Rec"] = r
    df_lst.append(df)


df = pd.concat(df_lst, axis=0)
df = (
    df[["Rec", "Chrom", "Pos", "AN", "AC"]]
    .sort_values(["Rec", "Chrom", "Pos"])
    .reset_index(drop=True)
)

df["MAC"] = df["AC"]
sel = df["AC"] * 2 > df["AN"]
df.loc[sel, "MAC"] = df.loc[sel, "AN"] - df.loc[sel, "AC"]

df["Maf"] = df["MAC"] / df["AN"]


df2 = (
    df[df["Maf"] > 0.01]
    .groupby(["Rec", "Chrom"])["Pos"]
    .count()
    .rename("Count")
    .reset_index()
    .groupby(["Rec"])["Count"]
    .agg(["mean", "std"])
    .reset_index()
    .assign(Rec=lambda df: df.Rec * 1e-9)
)

"""
%config inlinebackend.figure_format="retina"
"""
df2.to_csv("data_to_plot.csv")

fig, ax = plt.subplots(constrained_layout=True)
# ax.plot(df2["Rec"], df2["mean"], marker=".", ms=3, color="k")
ax.errorbar(df2["Rec"], y=df2["mean"], yerr=df2["std"], marker=".", ms=20, color="k")
ax.set_xlabel("Recombination Rate ($10^{-9}$)")
ax.set_ylabel("No. common variants (100cM)")
ax.axvline(x=667e-9, label="Pf", linestyle=":", color="r")
ax.axvline(x=10e-9, label="Hs", linestyle=":", color="b")
ax.legend()
ax.set_xscale("log")
ax.set_yscale("log")
fig.savefig(f"snp_dns_via_{model}_model.png", dpi=600)

