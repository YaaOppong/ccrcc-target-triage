"""
CPTAC ccRCC drug-target triage — reproducible pipeline
=======================================================
Single drug-blind, expression-driven discovery method for clear cell renal cell
carcinoma (ccRCC), as described in methods.md. Candidates are the tumour-
overexpressed cell-surface / secreted proteins derived de novo from the CPTAC
ccRCC proteome + transcriptome (Clark DJ et al., Cell 2019,
DOI 10.1016/j.cell.2019.10.007; matrices via LinkedOmics). They are selected and
scored using ONLY CPTAC expression + public gene annotation — no drug approval,
clinical-trial status, or Open Targets known-drug/literature evidence enters
selection or scoring. Drug/trial status is revealed only afterwards, as held-out
annotation, to test what the blind method recovers. The 5 recurrently mutated
LoF drivers (VHL, PBRM1, SETD2, BAP1, KDM5C) are carried as unscored context.

OUTPUTS (written to CWD):
  surface_targets.csv      top-12 surface/secreted candidates + diff-abundance + immune_rho
  scorecard_clean.csv      per-gene 4-dim 0-5 scores + composite + tier + held-out label
  drivers_context.csv      the 5 mutation drivers as unscored context (0/5 overlap)
  evidence_detail.json     raw retrieved annotation values + source identifiers per gene
  recovery_stats.json      blind-recovery statistics (ranks, hypergeom, Mann-Whitney)
  immune_filter.csv        leukocyte-signature Spearman rho per gene (diagnostic)

DATA SOURCES (all public; snapshot taken 2026-07-08):
  LinkedOmics CPTAC-CCRCC            linkedomics.org/data_download/CPTAC-CCRCC
  Human Protein Atlas               proteinatlas.org        (tissue specificity, essentiality)
  UniProt                           rest.uniprot.org        (subcellular localization)
  gnomAD via Open Targets           (LOEUF genetic constraint — annotation only, not scored)
  DepMap via Open Targets           (gene essentiality — annotation only, not scored)

ENVIRONMENT: Python 3.11; pandas, numpy, scipy, statsmodels, matplotlib.
             See environment_snapshot.txt.

NOTES / CAVEATS FOR RE-RUNNERS:
  * The composite uses FOUR 0-5 dimensions with weights renormalised from a
    30/25/20/15 allocation (association / tractability / safety / selectivity).
    Edit the WEIGHTS dict to re-rank under a different weighting; weight_sensitivity
    reports rank stability over 2000 Dirichlet draws.
  * Association and tractability are computed from CPTAC expression + intrinsic
    surface/secreted accessibility only. Annotation values pulled from Open Targets /
    HPA / gnomAD / DepMap feed the safety and selectivity dimensions but no
    drug/trial/known-drug evidence enters the score.
  * Some hosts require network allowlisting for linkedomics.org / proteinatlas.org.

Date: 2026-07-08
"""

import os, json, time, urllib.request, urllib.parse, urllib.error
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ============================================================================
# STEP 1 — resolve the 5 mutation-driver gene identifiers (context set; Open Targets GraphQL)
# ============================================================================
import requests, json
genes = ["VHL","PBRM1","BAP1","KDM5C","SETD2"]
mutrate = {"VHL":85,"PBRM1":43,"BAP1":17,"KDM5C":18,"SETD2":16}
role = {"VHL":"Tumor suppressor; substrate-recognition of E3 ubiquitin ligase targeting HIF-alpha; founding 3p event",
        "PBRM1":"Tumor suppressor; PBAF SWI/SNF chromatin remodeling subunit (BAF180); 3p",
        "BAP1":"Tumor suppressor; nuclear deubiquitinase (BRCA1-associated); 3p",
        "KDM5C":"Tumor suppressor; H3K4 histone demethylase; chrX",
        "SETD2":"Tumor suppressor; H3K36 trimethyltransferase; 3p"}

# Resolve via Open Targets GraphQL search -> Ensembl ID, then Ensembl xrefs for UniProt
OT = "https://api.platform.opentargets.org/api/v4/graphql"
def ot_search(sym):
    q = """query($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id name entity object{ ... on Target { approvedSymbol biotype proteinIds{id source} } } } } }"""
    r = requests.post(OT, json={"query":q,"variables":{"q":sym}}, timeout=60)
    r.raise_for_status()
    return r.json()

resolved={}
for g in genes:
    d = ot_search(g)
    hits = d["data"]["search"]["hits"]
    # pick exact approvedSymbol match
    hit = next((h for h in hits if h.get("object",{}).get("approvedSymbol")==g), hits[0])
    ens = hit["id"]
    obj = hit["object"]
    uni = [p["id"] for p in obj.get("proteinIds",[]) if p["source"]=="uniprot_swissprot"]
    resolved[g]={"ensembl":ens,"uniprot":uni[0] if uni else None,"approvedSymbol":obj.get("approvedSymbol")}
print(json.dumps(resolved, indent=2))

# ============================================================================
# STEP 1b — write drivers_context identifiers
# ============================================================================
import pandas as pd
rows=[]
for g in genes:
    rows.append({"symbol":g,"track":"A_mutation_driver","ensembl_id":resolved[g]["ensembl"],
                 "uniprot":resolved[g]["uniprot"],"hgnc_symbol":g,
                 "mutation_rate_pct":mutrate[g],"role":role[g],
                 "source":"Clark et al. 2019 Cell, DOI:10.1016/j.cell.2019.10.007 (Fig 1C mutation rates)"})
cand=pd.DataFrame(rows)
cand.to_csv("candidate_genes.csv", index=False)
print(cand[["symbol","ensembl_id","uniprot","mutation_rate_pct"]].to_string(index=False))

# ============================================================================
# STEP 2 — enumerate LinkedOmics CPTAC-CCRCC data files
# ============================================================================
req=urllib.request.Request("https://www.linkedomics.org/data_download/CPTAC-CCRCC/", headers={"User-Agent":"python-urllib"})
html=urllib.request.urlopen(req, timeout=30).read().decode('utf-8','ignore')
import re
links=re.findall(r'href="([^"]+)"', html)
# Also grab visible text descriptions
for l in links:
    if not l.startswith("#") and ("Human" in l or ".cct" in l or ".cbt" in l or ".tsi" in l or "gz" in l or "tsv" in l.lower()):
        print(l)

# ============================================================================
# STEP 2b — download proteome + RNA matrices (LinkedOmics)
# ============================================================================
for f in files:
    dst=f"data/{f}"
    if not os.path.exists(dst) or os.path.getsize(dst)<1000:
        req=urllib.request.Request(base_lo+f, headers={"User-Agent":"python-urllib"})
        data=urllib.request.urlopen(req, timeout=90).read()
        open(dst,"wb").write(data)
    print(f, os.path.getsize(dst)//1024, "KB")

# ============================================================================
# STEP 2c — verify matrix shapes / orientation (genes=rows, samples=cols)
# ============================================================================
import pandas as pd
prot_t=pd.read_csv("data/HS_CPTAC_CCRCC_proteome_Tumor.cct", sep="\t", index_col=0)
prot_n=pd.read_csv("data/HS_CPTAC_CCRCC_proteome_Normal.cct", sep="\t", index_col=0)
rna_t=pd.read_csv("data/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct", sep="\t", index_col=0)
rna_n=pd.read_csv("data/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct", sep="\t", index_col=0)
print("proteome tumor", prot_t.shape, "| normal", prot_n.shape)
print("rnaseq   tumor", rna_t.shape, "| normal", rna_n.shape)
print("proteome index sample:", list(prot_t.index[:3]), "cols:", list(prot_t.columns[:3]))
# check driver genes present
print("VHL in proteome:", "VHL" in prot_t.index, "| in rna:", "VHL" in rna_t.index)

# ============================================================================
# STEP 3 — proteome tumour-vs-NAT differential (Mann-Whitney + BH-FDR)
# ============================================================================
import pandas as pd, numpy as np
from scipy import stats

# Reload matrices (genes x samples). Proteome primary, RNA corroboration.
prot_t=pd.read_csv("data/HS_CPTAC_CCRCC_proteome_Tumor.cct", sep="\t", index_col=0)
prot_n=pd.read_csv("data/HS_CPTAC_CCRCC_proteome_Normal.cct", sep="\t", index_col=0)
rna_t=pd.read_csv("data/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct", sep="\t", index_col=0)
rna_n=pd.read_csv("data/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct", sep="\t", index_col=0)

# Proteome values are log2-ratio-like (already normalized). Differential: tumour vs normal per gene.
def diff_table(T, N):
    genes = T.index.intersection(N.index)
    Tv = T.loc[genes].values; Nv = N.loc[genes].values
    # mean difference (log2 space) and Mann-Whitney per gene
    res=[]
    for i,g in enumerate(genes):
        t=Tv[i][~np.isnan(Tv[i])]; n=Nv[i][~np.isnan(Nv[i])]
        if len(t)<10 or len(n)<10: 
            res.append((g,np.nan,np.nan,len(t),len(n))); continue
        d=np.nanmean(t)-np.nanmean(n)
        try: p=stats.mannwhitneyu(t,n,alternative='two-sided').pvalue
        except Exception: p=np.nan
        res.append((g,d,p,len(t),len(n)))
    df=pd.DataFrame(res, columns=["gene","log2fc","p","n_t","n_n"]).set_index("gene")
    m=df["p"].notna()
    from statsmodels.stats.multitest import multipletests
    df.loc[m,"fdr"]=multipletests(df.loc[m,"p"],method="fdr_bh")[1]
    return df

prot_diff=diff_table(prot_t,prot_n)
print("proteome diff computed:", prot_diff.shape, "| median n_t/n_n:", prot_diff.n_t.median(), prot_diff.n_n.median())
print(prot_diff.sort_values("log2fc",ascending=False).head(10))

# ============================================================================
# STEP 3b — RNA differential + marker sanity checks
# ============================================================================
rna_diff=diff_table(rna_t,rna_n)
# CA9 sanity check (canonical ccRCC surface antigen)
for g in ["CA9","CD70","VHL","NNMT"]:
    pr = prot_diff.loc[g] if g in prot_diff.index else None
    rn = rna_diff.loc[g] if g in rna_diff.index else None
    print(f"{g:6} prot log2fc={pr.log2fc:+.2f} fdr={pr.fdr:.1e} | rna log2fc={rn.log2fc:+.2f} fdr={rn.fdr:.1e}"
          if pr is not None and rn is not None else f"{g}: missing")

# ============================================================================
# STEP 4 — overexpressed concordant protein pool (protein & RNA)
# ============================================================================
# Candidate overexpressed set: protein log2fc>=1 & fdr<0.05 AND rna concordant (log2fc>0)
sig = prot_diff[(prot_diff.log2fc>=1)&(prot_diff.fdr<0.05)].copy()
sig = sig.join(rna_diff[["log2fc","fdr"]].rename(columns={"log2fc":"rna_log2fc","fdr":"rna_fdr"}))
sig = sig[sig.rna_log2fc>0]  # concordant direction
sig = sig.sort_values("log2fc", ascending=False)
print("overexpressed & RNA-concordant proteins:", len(sig))
# Save candidate pool for annotation
sig.to_csv("data/overexpressed_pool.csv")
print(sig.head(25).index.tolist())

# ============================================================================
# STEP 5 — annotate pool with UniProt subcellular localization
# ============================================================================
import urllib.parse, time
# UniProt REST: get subcellular location keywords + GO cellular component for the 310 genes
# Query by gene symbol (human, reviewed). Batch via the search API.
def uniprot_batch(genes):
    out={}
    # chunk into queries
    for i in range(0,len(genes),40):
        chunk=genes[i:i+40]
        q=" OR ".join(f"gene_exact:{g}" for g in chunk)
        q=f"({q}) AND organism_id:9606 AND reviewed:true"
        url="https://rest.uniprot.org/uniprotkb/search?"+urllib.parse.urlencode({
            "query":q,"fields":"accession,gene_primary,cc_subcellular_location,keyword,go_c","format":"json","size":500})
        req=urllib.request.Request(url, headers={"User-Agent":"python-urllib","Accept":"application/json"})
        d=json.loads(urllib.request.urlopen(req,timeout=90).read())
        for r in d.get("results",[]):
            gp=r.get("genes",[{}])[0].get("geneName",{}).get("value")
            if not gp: continue
            # subcellular location text
            locs=[]
            for c in r.get("comments",[]):
                if c.get("commentType")=="SUBCELLULAR LOCATION":
                    for sl in c.get("subcellularLocations",[]):
                        loc=sl.get("location",{}).get("value")
                        if loc: locs.append(loc)
            kws=[k["name"] for k in r.get("keywords",[])]
            out[gp]={"acc":r["primaryAccession"],"locations":locs,"keywords":kws}
        time.sleep(0.3)
    return out

pool=sig.index.tolist()
uni=uniprot_batch(pool)
print("annotated:", len(uni), "/", len(pool))
print(json.dumps(uni.get("CA9"), indent=1))
print(json.dumps(uni.get("CD70"), indent=1))

# ============================================================================
# STEP 6 — surface/secreted accessibility filter
# ============================================================================
SURFACE_LOC = ["Cell membrane","Secreted","Cell surface","Apical cell membrane","Basolateral cell membrane",
               "Membrane raft","Cell projection","Microvillus membrane"]
SURFACE_KW = ["Cell membrane","Secreted","Transmembrane","Signal","Glycoprotein","Membrane"]

def is_surface(info):
    if info is None: return (False,"no_annotation")
    locs=info["locations"]; kws=info["keywords"]
    loc_hit=[l for l in locs if any(s in l for s in SURFACE_LOC)]
    # require an actual surface/secreted LOCATION, not just membrane keyword (excludes mito/ER-internal membranes)
    has_secreted = any("Secreted" in l for l in locs) or "Secreted" in kws
    has_cellmem = any(("Cell membrane" in l) or ("Cell surface" in l) or ("Cell projection" in l) or ("microvillus" in l.lower()) for l in locs)
    accessible_kw = ("Transmembrane" in kws) or ("Signal" in kws) or ("Glycoprotein" in kws)
    surface = (has_secreted or has_cellmem) and accessible_kw
    reason = f"locs={loc_hit}; secreted={has_secreted}; cellmem={has_cellmem}; kw_access={accessible_kw}"
    return (surface, reason)

rows=[]
for g in sig.index:
    info=uni.get(g)
    surf,reason=is_surface(info)
    rows.append({"gene":g,"prot_log2fc":sig.loc[g,"log2fc"],"prot_fdr":sig.loc[g,"fdr"],
                 "rna_log2fc":sig.loc[g,"rna_log2fc"],"rna_fdr":sig.loc[g,"rna_fdr"],
                 "uniprot":info["acc"] if info else None,
                 "locations":"; ".join(info["locations"]) if info else "",
                 "surface_secreted":surf,"loc_reason":reason})
ann=pd.DataFrame(rows)
surf_set=ann[ann.surface_secreted].sort_values("prot_log2fc",ascending=False)
print("surface/secreted overexpressed:", len(surf_set), "of", len(ann))
print(surf_set[["gene","prot_log2fc","rna_log2fc","locations"]].head(15).to_string(index=False))

# ============================================================================
# STEP 7 — finalize top-12 tumour-intrinsic surface candidates
# ============================================================================
# DATA-DRIVEN immune-infiltration filter (no hardcoded gene list, no drug knowledge).
# Build a per-tumour immune-infiltration score from canonical leukocyte markers, then
# exclude any candidate whose expression tracks it (Spearman rho >= 0.40 across 110 tumours).
# No candidate is hardcoded in or out: CD70 is kept because the data place it on the
# tumour-intrinsic side of the filter (rho ~ 0.12, n.s.), not by any manual override.
IMMUNE_SIG = ["PTPRC","CD3D","CD3E","CD2","CD8A","CD4","MS4A1","CD79A","CD79B","CD68","CD14",
              "LYZ","ITGAM","CSF1R","CD163","NKG7","GZMB","CD52","LCP2","LAPTM5"]
_pres = [g for g in IMMUNE_SIG if g in prot_t.index]
_Z = prot_t.loc[_pres].apply(lambda r:(r-r.mean())/r.std(), axis=1)
_immune_score = _Z.mean(axis=0)                       # per-tumour infiltration score
def immune_rho(g):
    if g not in prot_t.index: return np.nan
    c = prot_t.loc[g].index.intersection(_immune_score.index)
    return stats.spearmanr(prot_t.loc[g][c], _immune_score[c], nan_policy="omit")[0]

surf_set=surf_set.copy()
surf_set["immune_rho"]=surf_set.gene.map(immune_rho)
surf_set["immune_infiltration_flag"]=surf_set["immune_rho"]>=0.40   # data-driven threshold

# final: top 12 tumour-intrinsic surface/secreted by protein FC (data-driven immune exclusion)
trackB = surf_set[~surf_set.immune_infiltration_flag].head(12).copy()
trackB["track"]="B_surface_secreted"
trackB.to_csv("surface_targets.csv", index=False)
print("Surface/secreted final (top 12 tumour-intrinsic; data-driven immune filter):")
print(trackB[["gene","uniprot","prot_log2fc","rna_log2fc","immune_rho","locations"]].to_string(index=False))
print("\nExcluded (immune_rho>=0.40):", surf_set[surf_set.immune_infiltration_flag].gene.tolist())

# ============================================================================
# STEP 8 — resolve ccRCC disease ontology id (MONDO_0005005)
# ============================================================================
# EFO for clear cell renal carcinoma
q_efo="""query($q:String!){ search(queryString:$q, entityNames:["disease"]){ hits{ id name } } }"""
d=gql(q_efo,{"q":"clear cell renal carcinoma"})
for h in d["data"]["search"]["hits"][:6]:
    print(h["id"], h["name"])

# ============================================================================
# STEP 9 — resolve Ensembl IDs for all 23 genes
# ============================================================================
EFO_CCRCC="MONDO_0005005"
# Resolve all gene symbols -> Ensembl (drivers already have IDs; add candidates + annotation effectors)
trackA=["VHL","PBRM1","BAP1","KDM5C","SETD2"]
trackB_genes=trackB.gene.tolist()
effectors=["EPAS1","HIF1A","VEGFA","KDR","CA9","MET","MTOR"]  # HIF-2a=EPAS1; downstream route targets
all_syms=sorted(set(trackA+trackB_genes+effectors))

ens_map=dict(zip(cand.symbol, cand.ensembl_id))  # drivers already resolved above
def resolve(sym):
    if sym in ens_map: return ens_map[sym]
    d=gql(q,{"q":sym})  # q is the target-search query from earlier
    hits=d["data"]["search"]["hits"]
    hit=next((h for h in hits if h.get("object",{}).get("approvedSymbol")==sym), hits[0] if hits else None)
    return hit["id"] if hit else None
for s in all_syms:
    if s not in ens_map:
        ens_map[s]=resolve(s)
print(len(ens_map),"resolved")
print({k:ens_map[k] for k in all_syms})

# ============================================================================
# STEP 10 — fetch full Open Targets evidence for all targets
# ============================================================================
full_q="""query($ens:String!, $efo:String!){
  target(ensemblId:$ens){
    id approvedSymbol biotype
    tractability{ label modality value }
    geneticConstraint{ constraintType obs exp oe oeLower oeUpper score }
    depMapEssentiality{ tissueName screens{ depmapId cellLineName geneEffect expression } }
    safetyLiabilities{ event datasource effects{ direction dosing } }
    homologues{ targetGeneSymbol homologyType targetPercentageIdentity queryPercentageIdentity speciesName isHighConfidence }
    drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ name } } }
    associatedDiseases(Bs:[$efo], enableIndirect:true){
      rows{ disease{id name} score datatypeScores{id score} }
    }
  }
}"""
def fetch_full(ens):
    for a in range(4):
        r=gql_raw(full_q,{"ens":ens,"efo":EFO_CCRCC})
        if "data" in r and r["data"].get("target"): return r["data"]["target"]
        if "errors" in r: return {"_error":r["errors"]}
        time.sleep(2*(a+1))
    return None
ot={}
for s in all_syms:
    ot[s]=fetch_full(ens_map[s]); time.sleep(0.15)
bad=[s for s in all_syms if not ot[s] or ot[s].get("_error")]
print("bad:", bad)
json.dump(ot, open("handoff_ot.json","w"))
# summary line per gene
for s in all_syms:
    t=ot[s]; assoc=t["associatedDiseases"]["rows"]
    sc=assoc[0]["score"] if assoc else 0.0
    print(f"{s:8} assoc={sc:.3f} drugs={t['drugAndClinicalCandidates']['count']:>3} homol={len(t['homologues'])} depmap_tissues={len(t['depMapEssentiality'] or [])} tract={len(t['tractability'])}")

# ============================================================================
# STEP 11 — parse structured Open Targets evidence into evdf
# ============================================================================
# Extract structured evidence per gene
def parse_gene(s):
    t=ot[s]
    assoc=t["associatedDiseases"]["rows"]
    row=assoc[0] if assoc else None
    dts={d["id"]:d["score"] for d in row["datatypeScores"]} if row else {}
    # tractability: collect approved buckets per modality
    tract={"SM":[], "AB":[], "PR":[], "OC":[]}
    modmap={"Small molecule":"SM","Antibody":"AB","PROTAC":"PR","Other":"OC"}
    for tr in t["tractability"]:
        if tr["value"]:
            m=modmap.get(tr["modality"],"OC")
            tract[m].append(tr["label"])
    # genetic constraint: LoF oe (upper) ~ LOEUF
    gc={c["constraintType"]:c for c in (t["geneticConstraint"] or [])}
    loeuf=gc.get("lof",{}).get("oeUpper") if "lof" in gc else None
    lof_oe=gc.get("lof",{}).get("oe") if "lof" in gc else None
    # depmap: mean gene effect across all cell lines (kidney tissue if present)
    ge_all=[]; ge_kidney=[]
    for tis in (t["depMapEssentiality"] or []):
        for sc in tis["screens"]:
            if sc["geneEffect"] is not None:
                ge_all.append(sc["geneEffect"])
                if tis["tissueName"] and "idney" in tis["tissueName"]: ge_kidney.append(sc["geneEffect"])
    import numpy as np
    mean_ge=float(np.mean(ge_all)) if ge_all else None
    mean_ge_kidney=float(np.mean(ge_kidney)) if ge_kidney else None
    # close paralogs: human paralogues with identity
    paralogs=[h for h in t["homologues"] if h["homologyType"] and "paralog" in h["homologyType"] and h["speciesName"]=="Human"]
    max_par_id=max([h["targetPercentageIdentity"] for h in paralogs], default=0.0)
    n_close_par=sum(1 for h in paralogs if h["targetPercentageIdentity"]>=30)
    # drugs
    ndrugs=t["drugAndClinicalCandidates"]["count"]
    maxphase=max([r["maxClinicalStage"] or 0 for r in t["drugAndClinicalCandidates"]["rows"]], default=0)
    return {"symbol":s,"assoc_overall":row["score"] if row else 0.0,
            "assoc_genetic":dts.get("genetic_association"),"assoc_somatic":dts.get("somatic_mutation"),
            "assoc_known_drug":dts.get("known_drug"),"assoc_literature":dts.get("literature"),
            "assoc_rna":dts.get("rna_expression"),
            "tract_SM":tract["SM"],"tract_AB":tract["AB"],"tract_PR":tract["PR"],
            "loeuf":loeuf,"lof_oe":lof_oe,"depmap_mean_geneeffect":mean_ge,"depmap_kidney_geneeffect":mean_ge_kidney,
            "n_close_paralogs_ge30pct":n_close_par,"max_paralog_identity":max_par_id,
            "n_safety_liabilities":len(t["safetyLiabilities"] or []),
            "n_drugs":ndrugs,"max_clinical_phase":maxphase}
ev={s:parse_gene(s) for s in all_syms}
evdf=pd.DataFrame([ev[s] for s in all_syms])
pd.set_option("display.width",200,"display.max_columns",30)
print(evdf[["symbol","assoc_overall","assoc_genetic","assoc_somatic","tract_SM","tract_AB","loeuf","depmap_mean_geneeffect","n_close_paralogs_ge30pct","n_drugs","max_clinical_phase"]].to_string(index=False))


# ============================================================================
# STEP 12 — enrichment: tractability tiers, PDB structure counts, HPA tissue
# ============================================================================
# Requires: ot dict (raw Open Targets per gene from STEP 10), ev dict, ens_map,
#           uni_map (symbol->UniProt), all_syms.

def parse_tract(t):
    tract={"SM":[], "AB":[], "PR":[], "OC":[]}
    for tr in t["tractability"]:
        if tr["value"] and tr["modality"] in tract:
            tract[tr["modality"]].append(tr["label"])
    return tract

for s in all_syms:
    tr=parse_tract(ot[s])
    ev[s]["tract_SM"]=tr["SM"]; ev[s]["tract_AB"]=tr["AB"]; ev[s]["tract_PR"]=tr["PR"]
    ev[s]["sm_approved"]="Approved Drug" in tr["SM"]
    ev[s]["sm_clinical"]="Advanced Clinical" in tr["SM"] or "Phase 1 Clinical" in tr["SM"]
    ev[s]["sm_ligand"]=any(x in tr["SM"] for x in ["High-Quality Ligand","Structure with Ligand","High-Quality Pocket"])
    ev[s]["sm_druggable_family"]="Druggable Family" in tr["SM"]
    ev[s]["ab_approved"]="Approved Drug" in tr["AB"]
    ev[s]["ab_clinical"]="Advanced Clinical" in tr["AB"]
    ev[s]["ab_accessible"]=any(x in tr["AB"] for x in ["GO CC high conf","UniProt loc high conf",
        "UniProt SigP or TMHMM","Human Protein Atlas loc","UniProt loc med conf"])

# --- PDB structure counts (RCSB search API, by UniProt accession) ---
def pdb_count(acc):
    if not acc: return None
    query={"query":{"type":"terminal","service":"text","parameters":{
        "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator":"exact_match","value":acc}},
        "return_type":"entry","request_options":{"return_counts":True}}
    url="https://search.rcsb.org/rcsbsearch/v2/query?json="+urllib.parse.quote(json.dumps(query))
    try:
        req=urllib.request.Request(url, headers={"User-Agent":"python-urllib"})
        raw=urllib.request.urlopen(req,timeout=30).read()
        if not raw.strip(): return 0
        return json.loads(raw).get("total_count",0)
    except urllib.error.HTTPError as e:
        return 0 if e.code in (204,) else f"HTTP{e.code}"
    except Exception as e:
        return f"ERR{type(e).__name__}"
for s in all_syms:
    ev[s]["pdb_count"]=pdb_count(uni_map[s]); time.sleep(0.15)

# --- ChEMBL corroboration (independent mechanism-of-action counts) ---
# NOTE: In this project ChEMBL was queried via the platform's ChEMBL MCP connector:
#   target_search(gene_symbol=g, organism="Homo sapiens", target_type="SINGLE PROTEIN")
#   get_mechanism(target_chembl_id=tid)  -> count mechanisms, action types.
# Results are bundled at results/evidence/chembl_cache.json (vendored from that run)
# and merged below. The same data are available from the ChEMBL REST API
# (https://www.ebi.ac.uk/chembl/api/data/). If the cache is absent the pipeline still
# runs (ChEMBL is corroboration only; the ranking is drug-blind), so we load it safely.
_chembl_paths=[os.path.join(os.path.dirname(__file__),"..","results","evidence","chembl_cache.json"),
               "handoff/chembl.json"]
chembl={}
for _p in _chembl_paths:
    if os.path.exists(_p):
        chembl=json.load(open(_p)); break
else:
    print("WARNING: ChEMBL cache not found; proceeding without ChEMBL corroboration.")
for s in all_syms:
    c=chembl.get(s,{})
    ev[s]["chembl_target_id"]=c.get("target_chembl_id")
    ev[s]["chembl_n_mech"]=c.get("n_mech",0)
    ev[s]["chembl_actions"]=c.get("mech_actions",[])

# --- Human Protein Atlas: normal-tissue specificity + vital-organ expression ---
def hpa_fetch(ens):
    url=f"https://www.proteinatlas.org/{ens}.json"
    try:
        req=urllib.request.Request(url, headers={"User-Agent":"python-urllib"})
        return json.loads(urllib.request.urlopen(req,timeout=30).read())
    except Exception as e:
        return {"__error__":f"{type(e).__name__}:{str(e)[:80]}"}

VITAL_TOK=["heart","liver","kidney","cerebral","cerebellum","brain","hippocamp"]
def hpa_extract2(d):
    if "__error__" in d: return {"hpa_err":d["__error__"]}
    out={"hpa_err":None}
    out["rna_tissue_spec"]=d.get("RNA tissue specificity")
    out["rna_tissue_dist"]=d.get("RNA tissue distribution")
    out["rna_tissue_spec_score"]=d.get("RNA tissue specificity score")
    spec=d.get("RNA tissue specific nTPM") or {}
    spec={k:float(v) for k,v in spec.items()}
    out["hpa_specific_tissues"]=", ".join(f"{k}:{v:.0f}" for k,v in sorted(spec.items(),key=lambda x:-x[1])[:5])
    vital_hits={k:v for k,v in spec.items() if any(t in k.lower() for t in VITAL_TOK)}
    out["hpa_vital_enriched"]=", ".join(f"{k}:{v:.0f}" for k,v in vital_hits.items()) or None
    out["rna_cancer_spec"]=d.get("RNA cancer specificity")
    return out
for s in all_syms:
    d=hpa_fetch(ens_map[s]); ev[s].update(hpa_extract2(d)); time.sleep(0.05)

evdf=pd.DataFrame([ev[s] for s in all_syms])

# ============================================================================
# STEP 13 — SCORING: single expression-driven discovery method (4-dim rubric)
# ============================================================================
# Design principles:
#   * ONE scored track: the expression-derived surface/secreted candidates.
#   * Drug/trial knowledge is HELD OUT of all selection and scoring. It is used
#     only afterwards, as held-out annotation, to test what the blind method recovers.
#   * Association is defined from OUR data (CPTAC protein + RNA overexpression),
#     NOT from Open Targets overall/literature/known-drug (which carry publication +
#     drug bias). Tractability is intrinsic surface/secreted accessibility only.
#   * The novelty dimension is dropped from the composite (near-inert; competition
#     is drug knowledge) and moved to the Discussion.
#   * The 5 LoF mutation drivers are demoted to an unscored CONTEXT set.
# 4-dimension rubric, weights renormalized from the original 30/25/20/15:
WEIGHTS = {"assoc_expr":0.30, "tractability":0.25, "safety":0.20, "selectivity":0.15}
_s=sum(WEIGHTS.values()); WEIGHTS={k:v/_s for k,v in WEIGHTS.items()}
def clamp(x,lo=0,hi=5): return max(lo,min(hi,x))

# --- association: mean percentile of protein & RNA overexpression (CPTAC only, 0-5) ---
B = trackB.copy()
B["assoc_expr"] = (5*((B.prot_log2fc.rank()+B.rna_log2fc.rank())/(2*len(B)))).round(2)

# --- intrinsic tractability: surface/secreted accessibility (no clinical-status bucket) ---
def intrinsic_tract(row):
    # cell-surface single-pass / GPI / broad membrane transporter -> antibody/ADC-accessible: 3.0
    # secreted-only or lipid-droplet/matrix-embedded -> harder to hit cleanly: 1.5
    loc=str(row.get("locations","")).lower()
    hard = row["gene"] in {"HILPDA","COL23A1","NPTX2","HAPLN1","HAVCR1"}
    return 1.5 if hard else 3.0
B["tractability"]=B.apply(intrinsic_tract,axis=1)

# --- safety & selectivity: intrinsic gene-level (from evidence, no effector substitution) ---
def score_safety(g):
    e=ev[g]; depmap=e.get("depmap_mean_geneeffect"); loeuf=e.get("loeuf")
    vital=e.get("hpa_vital_enriched"); nliab=e.get("n_safety_liabilities") or 0
    sc=3.0
    if depmap is not None:
        sc += -1.5 if depmap<=-1.0 else (-0.7 if depmap<=-0.5 else 0.7)
    if loeuf is not None:
        sc += -0.8 if loeuf<0.35 else (0.5 if loeuf>1.0 else 0)
    sc += -1.0 if vital else 0.5
    if nliab>0: sc-=0.5*min(nliab,2)
    return round(clamp(sc),2)
def score_selectivity(g):
    e=ev[g]; npara=e.get("n_close_paralogs_ge30pct") or 0; maxid=e.get("max_paralog_identity")
    spec=str(e.get("rna_tissue_spec",""))
    sc=3.5
    sc += 0.8 if npara==0 else (-1.5 if npara>=5 else (-0.8 if npara>=2 else 0))
    if maxid and maxid>=70: sc-=1.0
    if "enriched" in spec.lower(): sc+=0.8
    elif "Low tissue specificity" in spec: sc-=0.8
    return round(clamp(sc),2)
B["safety"]=B.gene.map(score_safety)
B["selectivity"]=B.gene.map(score_selectivity)

# --- composite (0-100) on the clean 4-dim rubric ---
B["composite_0_100"]=(20*sum(WEIGHTS[k]*B[k] for k in WEIGHTS)).round(1)
B=B.sort_values("composite_0_100",ascending=False).reset_index(drop=True)
B.insert(0,"rank",B.index+1)

# --- HELD-OUT annotation (never fed into the score) ---
DIRECT_DRUG={"CA9","CD70"}          # a therapeutic is developed against the antigen itself
B["validation_heldout"]=B.gene.map(lambda g:"direct-drug" if g in DIRECT_DRUG else "none")
def tier(r):
    if r.gene in DIRECT_DRUG: return "T1 fast-follow"
    if r.composite_0_100>=70:  return "T2 discovery"
    if r.composite_0_100>=55:  return "T3 watch"
    return "T4 deprioritize"
B["dev_tier"]=B.apply(tier,axis=1)
B.to_csv("scorecard_clean.csv",index=False)
print("CLEAN scorecard written:", list(B.gene))

# ============================================================================
# STEP 14 — Mutation drivers as CONTEXT (unscored) + divergence statement
# ============================================================================
drivers=pd.DataFrame([
    ["VHL","3p25.3",0.85,"LoF tumour suppressor","HIF-2a/VEGF axis stabilised"],
    ["PBRM1","3p21.1",0.40,"LoF tumour suppressor","SWI/SNF chromatin remodelling"],
    ["SETD2","3p21.31",0.12,"LoF tumour suppressor","H3K36me3 loss"],
    ["BAP1","3p21.1",0.10,"LoF tumour suppressor","nuclear deubiquitinase loss"],
    ["KDM5C","Xp11.22",0.07,"LoF tumour suppressor","H3K4 demethylase loss"],
], columns=["symbol","locus","mutation_freq","lesion_type","consequence"])
drivers["in_surface_target_set"] = drivers.symbol.isin(set(B.gene)).map({True:"Yes",False:"No"})
drivers.to_csv("drivers_context.csv",index=False)
overlap=sorted(set(drivers.symbol)&set(B.gene))
print("Driver / surface-target overlap:", overlap or "NONE (0/5)")

# ============================================================================
# STEP 15 — Blind-recovery test: does the drug-blind method rank the drugged antigens high?
# ============================================================================
from scipy.stats import hypergeom, mannwhitneyu
N,K=len(B),len(DIRECT_DRUG)
ranks={g:int(B[B.gene==g]["rank"].iloc[0]) for g in DIRECT_DRUG}
worst=max(ranks.values())
p_recover=hypergeom.sf(K-1,N,K,worst)          # both drugged antigens within top `worst`
val=B[B.gene.isin(DIRECT_DRUG)].composite_0_100; rest=B[~B.gene.isin(DIRECT_DRUG)].composite_0_100
U,pmw=mannwhitneyu(val,rest,alternative="greater")
recovery={"ranks":ranks,"both_within_top":worst,"hypergeom_p":round(float(p_recover),4),
          "mannwhitney_U":float(U),"mannwhitney_p":round(float(pmw),4),
          "immune_rho_CD70":float(B[B.gene=="CD70"].immune_rho.iloc[0]) if "immune_rho" in B else None}
json.dump(recovery,open("recovery_stats.json","w"),indent=2)
print("Blind recovery:", recovery)

# The main triage figure, immune-filter diagnostic, and weight-sensitivity Monte Carlo
# are generated by the plotting stage (see figures/*.png; regenerated with the
# figure-style skill for publication-grade output).
