
from pathlib import Path
import pandas as pd

R=Path.home()/"Influenza_RSV_Project"
O=R/"results/regulatory_driver_analysis/qc"
O.mkdir(parents=True,exist_ok=True)

net=pd.read_csv(R/"results/regulatory_driver_analysis/design/COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz",sep="\t")
tot=net.groupby("source")["target"].nunique().rename("unique_targets_total")

rows=[]
for p in ["GPL6884","GPL10558"]:
    f=R/f"data/processed/DS001_GSE38900/analysis_ready/{p}_gene_expression.tsv.gz"
    genes=set(pd.read_csv(f,sep="\t",usecols=["gene_symbol"])["gene_symbol"].dropna().astype(str))
    m=net[net["target"].isin(genes)].groupby("source")["target"].nunique().rename("measurable_targets")
    d=pd.concat([tot,m],axis=1).fillna({"measurable_targets":0}).reset_index()
    d["measurable_targets"]=d["measurable_targets"].astype(int)
    d["coverage_fraction"]=d["measurable_targets"]/d["unique_targets_total"]
    d["eligible_min5"]=d["measurable_targets"]>=5
    d["eligible_min10"]=d["measurable_targets"]>=10
    d.insert(0,"platform",p)
    d.to_csv(O/f"{p}_COLLECTRI_TARGET_COVERAGE_v1.0.tsv",sep="\t",index=False)
    rows.append(d)

all_df=pd.concat(rows,ignore_index=True)
all_df.to_csv(O/"COLLECTRI_TARGET_COVERAGE_ALL_PLATFORMS_v1.0.tsv",sep="\t",index=False)

summary=(all_df.groupby("platform")
         .agg(expression_sources=("source","size"),
              sources_with_at_least_1_target=("measurable_targets",lambda x:(x>=1).sum()),
              sources_eligible_min5=("eligible_min5","sum"),
              sources_eligible_min10=("eligible_min10","sum"),
              median_measurable_targets=("measurable_targets","median"),
              median_coverage_fraction=("coverage_fraction","median"))
         .reset_index())
summary.to_csv(O/"COLLECTRI_TARGET_COVERAGE_SUMMARY_v1.0.tsv",sep="\t",index=False)

pv=all_df.pivot(index="source",columns="platform",values="measurable_targets").reset_index()
pv["eligible_min5_both"]=(pv["GPL6884"]>=5)&(pv["GPL10558"]>=5)
pv["eligible_min10_both"]=(pv["GPL6884"]>=10)&(pv["GPL10558"]>=10)
pv.to_csv(O/"COLLECTRI_CROSSPLATFORM_EVALUABILITY_v1.0.tsv",sep="\t",index=False)

print(summary.to_string(index=False))
print(">=5 both:",int(pv["eligible_min5_both"].sum()))
print(">=10 both:",int(pv["eligible_min10_both"].sum()))
print("COVERAGE QC COMPLETE")
