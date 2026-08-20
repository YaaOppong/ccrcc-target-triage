"""
CPTAC ccRCC drug-target triage — reproducible pipeline
=======================================================
Single drug-blind, expression-driven discovery method for clear cell renal cell
carcinoma (ccRCC), as described in docs/methods.md. Candidates are the tumour-
overexpressed cell-surface / secreted proteins derived de novo from the CPTAC
ccRCC proteome + transcriptome (Clark DJ et al., Cell 2019,
DOI 10.1016/j.cell.2019.10.007; matrices via LinkedOmics). They are selected and
scored using ONLY CPTAC expression + public gene annotation — no drug approval,
clinical-trial status, or Open Targets known-drug/literature evidence enters
selection or scoring. Drug/trial status is revealed only afterwards, as held-out
annotation, to test what the blind method recovers. The 5 recurrently mutated
LoF drivers (VHL, PBRM1, SETD2, BAP1, KDM5C) are carried as unscored context.

RUNNING
-------
    pip install -r requirements.txt
    python code/pipeline.py

Runs end to end from a clean checkout. Input matrices (~80 MB) download to data/
on first run; every API response is memoised under _cache/ so re-runs take
seconds and do not re-hit UniProt / Open Targets / HPA / RCSB / ClinicalTrials.gov.
Delete _cache/ to force a live refresh.

OUTPUTS (written under results/)
--------------------------------
  scorecards/scorecard_clean.csv      the 12 selected candidates, 4 dims + composite + tier
  scorecards/scorecard_101.csv        ALL 101 surface candidates scored on the same rubric
  scorecards/drivers_context.csv      the 5 mutation drivers as unscored context (0/5 overlap)
  scorecards/recovery_stats.json      held-out recovery statistics (see EVALUATION below)
  scorecards/surface_candidates_all.csv  the 101, with why each was dropped (audits 101 -> 12)
  evidence/evidence_detail.json       raw retrieved annotation values + sources per gene
  evidence/clinical_precedent.csv     held-out drug/clinical status per gene
  enrichment/immune_filter.csv        leukocyte-signature Spearman rho per gene (diagnostic)
  enrichment/weight_sensitivity.csv   rank stability over 2000 Dirichlet weight draws
  curated/clinical_agents.csv         held-out positive set for all 101 (see CURATED LAYER)

EVALUATION — what the recovery statistics do and do not test
------------------------------------------------------------
An earlier version of this pipeline tested recovery with a hypergeometric over the
12 FINAL candidates. That null is conditional on selection: both drugged antigens
have already survived the 11,710 -> 364 -> 310 -> 101 -> 12 cascade before the
statistic runs, so it could only ever measure the last ranking step, and a positive
dropped anywhere in the cascade was invisible to it. It is reported here as
`ranking_given_selection` and is NOT the headline number.

The pipeline now scores all 101 surface candidates and evaluates two distinct
questions against a held-out positive set defined by rule over all 101:
  * cascade_test  — did selection concentrate positives into the retained 12?
  * ranking_test  — does the composite rank positives above non-positives across
                    the whole 101, unconditional on selection?
Both can record a miss. See docs/results.md for the interpretation.

CURATED LAYER
-------------
curated/clinical_agents.csv is the held-out positive set. It is generated
mechanically here (Open Targets mechanism-of-action + ClinicalTrials.gov), but
the `agent_modality` column is known-incomplete: ChEMBL's text search does not
reliably identify antibody-drug conjugates (it returns "Unknown" for MDX-1411 and
AGS-16C3F, both of which are biologics). Treat that column as a curated input open
to correction; the pipeline re-runs against an edited file without code changes.
The two positive definitions used for the statistics do NOT depend on modality.

DATA SOURCES (all public)
  LinkedOmics CPTAC-CCRCC   linkedomics.org/data_download/CPTAC-CCRCC
  UniProt                   rest.uniprot.org       (subcellular localization)
  Open Targets              api.platform.opentargets.org (tractability, constraint, DepMap)
  Human Protein Atlas       proteinatlas.org       (tissue specificity)
  RCSB PDB                  search.rcsb.org        (structure counts)
  ClinicalTrials.gov v2     clinicaltrials.gov/api (held-out trial annotation)

ENVIRONMENT: Python 3.11; pandas, numpy, scipy, statsmodels. See environment_snapshot.txt.
"""

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CACHE = os.path.join(ROOT, "_cache")
CACHE_CT = os.path.join(ROOT, "_cache_ct")
RES = os.path.join(ROOT, "results")
CURATED = os.path.join(ROOT, "curated")
for d in (DATA, CACHE, CACHE_CT, CURATED,
          f"{RES}/scorecards", f"{RES}/evidence", f"{RES}/enrichment", f"{RES}/trials"):
    os.makedirs(d, exist_ok=True)

UA = {"User-Agent": "python-urllib"}


def cached(key, fn, cache_dir=CACHE):
    """Memoise one network call to JSON. Keeps the pipeline re-runnable offline."""
    p = os.path.join(cache_dir, f"{key}.json")
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    v = fn()
    with open(p, "w") as fh:
        json.dump(v, fh)
    return v


# ============================================================================
# STEP 1 — Open Targets GraphQL helpers
# ============================================================================
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"


def gql_raw(query, variables):
    """POST a GraphQL query; return the raw envelope (callers may inspect `errors`)."""
    body = json.dumps({"query": query, "variables": variables}).encode()
    req = urllib.request.Request(
        OT_API, data=body, headers={"Content-Type": "application/json", **UA})
    return json.loads(urllib.request.urlopen(req, timeout=60).read())


def gql(query, variables):
    """POST a GraphQL query; raise on API-level errors."""
    r = gql_raw(query, variables)
    if "errors" in r:
        raise RuntimeError(r["errors"])
    return r


Q_TARGET = """query($q:String!){ search(queryString:$q, entityNames:["target"]){
    hits{ id name object{ ... on Target { approvedSymbol proteinIds{id source} } } } } }"""

EFO_CCRCC = "MONDO_0005005"   # clear cell renal carcinoma

FULL_Q = """query($ens:String!, $efo:String!){
  target(ensemblId:$ens){
    id approvedSymbol biotype
    tractability{ label modality value }
    geneticConstraint{ constraintType obs exp oe oeLower oeUpper score }
    depMapEssentiality{ tissueName screens{ depmapId cellLineName geneEffect expression } }
    safetyLiabilities{ event datasource effects{ direction dosing } }
    homologues{ targetGeneSymbol homologyType targetPercentageIdentity queryPercentageIdentity speciesName isHighConfidence }
    drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ name id } } }
    associatedDiseases(Bs:[$efo], enableIndirect:true){
      rows{ disease{id name} score datatypeScores{id score} }
    }
  }
}"""


def resolve_target(sym):
    """Gene symbol -> (Ensembl id, SwissProt accession)."""
    d = cached(f"search_{sym}", lambda: gql(Q_TARGET, {"q": sym}))
    hits = d["data"]["search"]["hits"]
    hit = next((h for h in hits if (h.get("object") or {}).get("approvedSymbol") == sym),
               hits[0] if hits else None)
    if not hit:
        return None, None
    prot = [p["id"] for p in (hit.get("object") or {}).get("proteinIds", [])
            if p["source"] == "uniprot_swissprot"]
    return hit["id"], (prot[0] if prot else None)


# ============================================================================
# STEP 2 — enumerate + download CPTAC-CCRCC matrices (LinkedOmics)
# ============================================================================
base_lo = "https://www.linkedomics.org/data_download/CPTAC-CCRCC/"
_listing = cached("linkedomics_listing", lambda: {"html": urllib.request.urlopen(
    urllib.request.Request(base_lo, headers=UA), timeout=60).read().decode("utf-8", "ignore")})
import re  # noqa: E402  (used only for the listing scrape)

links = re.findall(r'href="([^"]+)"', _listing["html"])
WANTED = ["proteome_Tumor", "proteome_Normal",
          "RNAseq_fpkm_log2_Tumor", "RNAseq_fpkm_log2_Normal"]
files = sorted({os.path.basename(l) for l in links
                if l.endswith(".cct") and any(w in l for w in WANTED)})
assert len(files) == 4, f"expected 4 matrices in the LinkedOmics listing, found {files}"

for f in files:
    dst = os.path.join(DATA, f)
    if not os.path.exists(dst) or os.path.getsize(dst) < 1000:
        req = urllib.request.Request(base_lo + f, headers=UA)
        with open(dst, "wb") as fh:
            fh.write(urllib.request.urlopen(req, timeout=300).read())
    print(f"  {f} {os.path.getsize(dst)//1024} KB")

prot_t = pd.read_csv(f"{DATA}/HS_CPTAC_CCRCC_proteome_Tumor.cct", sep="\t", index_col=0)
prot_n = pd.read_csv(f"{DATA}/HS_CPTAC_CCRCC_proteome_Normal.cct", sep="\t", index_col=0)
rna_t = pd.read_csv(f"{DATA}/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Tumor.cct", sep="\t", index_col=0)
rna_n = pd.read_csv(f"{DATA}/HS_CPTAC_CCRCC_RNAseq_fpkm_log2_Normal.cct", sep="\t", index_col=0)
print(f"proteome tumour {prot_t.shape} | normal {prot_n.shape}")
print(f"rnaseq   tumour {rna_t.shape} | normal {rna_n.shape}")


# ============================================================================
# STEP 3 — tumour-vs-NAT differential abundance (Mann-Whitney + BH-FDR)
# ============================================================================
def diff_table(T, N):
    genes = T.index.intersection(N.index)
    Tv, Nv = T.loc[genes].values, N.loc[genes].values
    res = []
    for i, g in enumerate(genes):
        t = Tv[i][~np.isnan(Tv[i])]
        n = Nv[i][~np.isnan(Nv[i])]
        if len(t) < 10 or len(n) < 10:
            res.append((g, np.nan, np.nan, len(t), len(n)))
            continue
        d = np.nanmean(t) - np.nanmean(n)
        try:
            p = stats.mannwhitneyu(t, n, alternative="two-sided").pvalue
        except ValueError:
            p = np.nan
        res.append((g, d, p, len(t), len(n)))
    df = pd.DataFrame(res, columns=["gene", "log2fc", "p", "n_t", "n_n"]).set_index("gene")
    m = df["p"].notna()
    df.loc[m, "fdr"] = multipletests(df.loc[m, "p"], method="fdr_bh")[1]
    return df


prot_diff = diff_table(prot_t, prot_n)
rna_diff = diff_table(rna_t, rna_n)
print(f"proteins measured: {len(prot_diff)}")

# ============================================================================
# STEP 4 — overexpressed, RNA-concordant pool
# ============================================================================
sig = prot_diff[(prot_diff.log2fc >= 1) & (prot_diff.fdr < 0.05)].copy()
print(f"  -> protein-overexpressed (log2FC>=1, FDR<0.05): {len(sig)}")
sig = sig.join(rna_diff[["log2fc", "fdr"]].rename(
    columns={"log2fc": "rna_log2fc", "fdr": "rna_fdr"}))
sig = sig[sig.rna_log2fc > 0].sort_values("log2fc", ascending=False)
print(f"  -> + RNA-concordant: {len(sig)}")


# ============================================================================
# STEP 5 — UniProt subcellular localization for the pool
# ============================================================================
def uniprot_batch(genes):
    out = {}
    for i in range(0, len(genes), 40):
        chunk = genes[i:i + 40]
        q = " OR ".join(f"gene_exact:{g}" for g in chunk)
        q = f"({q}) AND organism_id:9606 AND reviewed:true"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
            "query": q,
            "fields": "accession,gene_primary,cc_subcellular_location,keyword,go_c",
            "format": "json", "size": 500})
        req = urllib.request.Request(url, headers={"Accept": "application/json", **UA})
        d = json.loads(urllib.request.urlopen(req, timeout=90).read())
        for r in d.get("results", []):
            gp = r.get("genes", [{}])[0].get("geneName", {}).get("value")
            if not gp:
                continue
            locs = []
            for c in r.get("comments", []):
                if c.get("commentType") == "SUBCELLULAR LOCATION":
                    for sl in c.get("subcellularLocations", []):
                        loc = sl.get("location", {}).get("value")
                        if loc:
                            locs.append(loc)
            out[gp] = {"acc": r["primaryAccession"], "locations": locs,
                       "keywords": [k["name"] for k in r.get("keywords", [])]}
        time.sleep(0.3)
    return out


pool = sig.index.tolist()
uni = cached("uniprot_pool", lambda: uniprot_batch(pool))
print(f"UniProt annotated: {len(uni)} / {len(pool)}")

# ============================================================================
# STEP 6 — surface / secreted accessibility filter
# ============================================================================
SURFACE_LOC = ["Cell membrane", "Secreted", "Cell surface", "Apical cell membrane",
               "Basolateral cell membrane", "Membrane raft", "Cell projection",
               "Microvillus membrane"]


def is_surface(info):
    """Require a genuine surface/secreted LOCATION, not merely a membrane keyword
    (which would admit mitochondrial and ER-internal membranes)."""
    if info is None:
        return False, "no_annotation"
    locs, kws = info["locations"], info["keywords"]
    loc_hit = [l for l in locs if any(s in l for s in SURFACE_LOC)]
    has_secreted = any("Secreted" in l for l in locs) or "Secreted" in kws
    has_cellmem = any(("Cell membrane" in l) or ("Cell surface" in l)
                      or ("Cell projection" in l) or ("microvillus" in l.lower())
                      for l in locs)
    accessible_kw = any(k in kws for k in ("Transmembrane", "Signal", "Glycoprotein"))
    surface = (has_secreted or has_cellmem) and accessible_kw
    return surface, (f"locs={loc_hit}; secreted={has_secreted}; "
                     f"cellmem={has_cellmem}; kw_access={accessible_kw}")


rows = []
for g in sig.index:
    info = uni.get(g)
    surf, reason = is_surface(info)
    rows.append({"gene": g, "prot_log2fc": sig.loc[g, "log2fc"], "prot_fdr": sig.loc[g, "fdr"],
                 "rna_log2fc": sig.loc[g, "rna_log2fc"], "rna_fdr": sig.loc[g, "rna_fdr"],
                 "uniprot": info["acc"] if info else None,
                 "locations": "; ".join(info["locations"]) if info else "",
                 "surface_secreted": surf, "loc_reason": reason})
ann = pd.DataFrame(rows)
surf_set = ann[ann.surface_secreted].sort_values("prot_log2fc", ascending=False).copy()
print(f"  -> + surface/secreted: {len(surf_set)}")

# ============================================================================
# STEP 7 — data-driven immune-infiltration filter, then the top-12 cut
# ============================================================================
# Bulk-tumour overexpression screens routinely surface markers of infiltrating
# leukocytes rather than tumour cells. Build a per-tumour infiltration score from
# canonical leukocyte markers and exclude candidates whose abundance tracks it.
# No candidate is hardcoded in or out.
IMMUNE_SIG = ["PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD4", "MS4A1", "CD79A", "CD79B",
              "CD68", "CD14", "LYZ", "ITGAM", "CSF1R", "CD163", "NKG7", "GZMB", "CD52",
              "LCP2", "LAPTM5"]
_pres = [g for g in IMMUNE_SIG if g in prot_t.index]
_Z = prot_t.loc[_pres].apply(lambda r: (r - r.mean()) / r.std(), axis=1)
_immune_score = _Z.mean(axis=0)


def immune_rho(g):
    if g not in prot_t.index:
        return np.nan
    c = prot_t.loc[g].index.intersection(_immune_score.index)
    return stats.spearmanr(prot_t.loc[g][c], _immune_score[c], nan_policy="omit")[0]


surf_set["immune_rho"] = surf_set.gene.map(immune_rho)
surf_set["immune_infiltration_flag"] = surf_set["immune_rho"] >= 0.40

trackB = surf_set[~surf_set.immune_infiltration_flag].head(12).copy()
trackB["track"] = "B_surface_secreted"
n_immune_dropped = int(surf_set.immune_infiltration_flag.sum())
print(f"  -> immune filter drops {n_immune_dropped}; "
      f"{len(surf_set) - n_immune_dropped} survive; top-12 by protein FC retained")
print("  final 12:", trackB.gene.tolist())

# Persist the full 101 so the 101 -> 12 step is auditable.
audit = surf_set.assign(
    dropped_by=np.where(surf_set.immune_infiltration_flag, "immune_filter",
                        np.where(surf_set.gene.isin(trackB.gene), "", "fc_rank_cutoff")),
    retained=surf_set.gene.isin(trackB.gene))
audit[["gene", "uniprot", "prot_log2fc", "rna_log2fc", "immune_rho",
       "immune_infiltration_flag", "dropped_by", "retained", "locations"]].to_csv(
    f"{RES}/scorecards/surface_candidates_all.csv", index=False)

# Immune-filter diagnostic: every surface candidate, plus tumour-intrinsic controls
# that are not themselves candidates. Candidates above the threshold are labelled
# `immune_filtered` — the filter's claim is that their abundance tracks infiltration,
# which is not the same as asserting they are leukocyte-lineage genes.
imm_rows = [{"group": ("immune_filtered" if r.immune_infiltration_flag else "surface_candidate"),
             "gene": r.gene, "immune_rho": round(float(r.immune_rho), 3)}
            for r in surf_set.itertuples() if pd.notna(r.immune_rho)]
_seen = {r["gene"] for r in imm_rows}
for g in ["NNMT", "VIM", "VEGFA", "ENO2", "ANGPTL4"]:      # tumour-intrinsic controls
    rho = immune_rho(g)
    if g not in _seen and pd.notna(rho):
        imm_rows.append({"group": "tumour_intrinsic_control", "gene": g,
                         "immune_rho": round(float(rho), 3)})
pd.DataFrame(imm_rows).sort_values("immune_rho", ascending=False).to_csv(
    f"{RES}/enrichment/immune_filter.csv", index=False)

# ============================================================================
# STEP 8 — evidence retrieval for ALL 101 surface candidates (+ drivers, effectors)
# ============================================================================
# The whole 101 is scored, not just the 12, so the selection cascade itself can be
# evaluated in STEP 12 rather than assumed correct.
DRIVERS = ["VHL", "PBRM1", "BAP1", "KDM5C", "SETD2"]
MUTRATE = {"VHL": 85, "PBRM1": 43, "BAP1": 17, "KDM5C": 18, "SETD2": 16}
EFFECTORS = ["EPAS1", "HIF1A", "VEGFA", "KDR", "MET", "MTOR"]
surface101 = surf_set.gene.tolist()
all_syms = sorted(set(surface101 + DRIVERS + EFFECTORS))
print(f"evidence universe: {len(all_syms)} genes")

ens_map, uni_map = {}, {}
for s in all_syms:
    ens_map[s], uni_map[s] = resolve_target(s)
unresolved = [s for s in all_syms if not ens_map[s]]
if unresolved:
    print("  WARNING unresolved:", unresolved)

ot = {}
for s in all_syms:
    if not ens_map[s]:
        ot[s] = None
        continue
    ot[s] = cached(f"ot_{s}", lambda s=s: (gql_raw(FULL_Q, {"ens": ens_map[s], "efo": EFO_CCRCC})
                                           .get("data") or {}).get("target"))


def hpa_fetch(ens):
    try:
        req = urllib.request.Request(f"https://www.proteinatlas.org/{ens}.json", headers=UA)
        return json.loads(urllib.request.urlopen(req, timeout=30).read())
    except Exception as e:                                    # noqa: BLE001
        return {"__error__": f"{type(e).__name__}:{str(e)[:80]}"}


hpa = {s: (cached(f"hpa_{s}", lambda s=s: hpa_fetch(ens_map[s])) if ens_map[s]
           else {"__error__": "no_ensembl"}) for s in all_syms}


def pdb_count(acc):
    if not acc:
        return None
    query = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": ("rcsb_polymer_entity_container_identifiers"
                      ".reference_sequence_identifiers.database_accession"),
        "operator": "exact_match", "value": acc}},
        "return_type": "entry", "request_options": {"return_counts": True}}
    url = ("https://search.rcsb.org/rcsbsearch/v2/query?json="
           + urllib.parse.quote(json.dumps(query)))
    try:
        raw = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30).read()
        return json.loads(raw).get("total_count", 0) if raw.strip() else 0
    except urllib.error.HTTPError as e:
        return 0 if e.code == 204 else None
    except Exception:                                         # noqa: BLE001
        return None


pdb = {s: cached(f"pdb_{s}", lambda s=s: {"n": pdb_count(uni_map.get(s))})["n"]
       for s in all_syms}

# ============================================================================
# STEP 9 — parse structured evidence
# ============================================================================
# NOTE: the Open Targets API reports tractability `modality` as short codes
# ("SM"/"AB"/"PR"). Long-name forms are accepted too so the parse survives either
# shape — a long-name-only map silently yields empty tractability lists.
MODMAP = {"Small molecule": "SM", "Antibody": "AB", "PROTAC": "PR", "Other": "OC",
          "SM": "SM", "AB": "AB", "PR": "PR", "OC": "OC"}
VITAL_TOK = ["heart", "liver", "kidney", "cerebral", "cerebellum", "brain", "hippocamp"]
PHASE = {"PHASE_1": 1, "PHASE_1_2": 1, "PHASE_2": 2, "PHASE_2_3": 2, "PHASE_3": 3,
         "PHASE_4": 4, "EARLY_PHASE_1": 1, "PRECLINICAL": 0, "NA": 0, None: 0}


def parse_gene(s):
    t = ot[s]
    assoc = t["associatedDiseases"]["rows"]
    row = assoc[0] if assoc else None
    dts = {d["id"]: d["score"] for d in row["datatypeScores"]} if row else {}
    tract = {"SM": [], "AB": [], "PR": [], "OC": []}
    for tr in t["tractability"]:
        if tr["value"]:
            tract[MODMAP.get(tr["modality"], "OC")].append(tr["label"])
    gc = {c["constraintType"]: c for c in (t["geneticConstraint"] or [])}
    ge_all, ge_kid = [], []
    for tis in (t["depMapEssentiality"] or []):
        for sc in tis["screens"]:
            if sc["geneEffect"] is not None:
                ge_all.append(sc["geneEffect"])
                if tis["tissueName"] and "idney" in tis["tissueName"]:
                    ge_kid.append(sc["geneEffect"])
    paralogs = [h for h in t["homologues"] if h["homologyType"]
                and "paralog" in h["homologyType"] and h["speciesName"] == "Human"]
    dc = t["drugAndClinicalCandidates"]
    d = {"symbol": s, "ensembl_id": ens_map[s], "uniprot": uni_map[s],
         "assoc_overall": row["score"] if row else 0.0,
         "assoc_genetic": dts.get("genetic_association"),
         "assoc_somatic": dts.get("somatic_mutation"),
         "tract_SM": tract["SM"], "tract_AB": tract["AB"], "tract_PR": tract["PR"],
         "loeuf": gc.get("lof", {}).get("oeUpper") if "lof" in gc else None,
         "depmap_mean_geneeffect": float(np.mean(ge_all)) if ge_all else None,
         "depmap_kidney_geneeffect": float(np.mean(ge_kid)) if ge_kid else None,
         "n_close_paralogs_ge30pct": sum(1 for h in paralogs
                                         if h["targetPercentageIdentity"] >= 30),
         "max_paralog_identity": max([h["targetPercentageIdentity"] for h in paralogs],
                                     default=0.0),
         "n_safety_liabilities": len(t["safetyLiabilities"] or []),
         "pdb_count": pdb.get(s),
         # ---- held out of scoring; used only for the recovery test ----
         "n_drugs": dc["count"],
         "max_clinical_phase": max([PHASE.get(r["maxClinicalStage"], 0)
                                    for r in dc["rows"]], default=0),
         "agents": sorted({r["drug"]["name"] for r in dc["rows"]})}
    h = hpa[s] or {}
    if "__error__" not in h:
        spec = {k: float(v) for k, v in (h.get("RNA tissue specific nTPM") or {}).items()}
        vital = {k: v for k, v in spec.items() if any(tk in k.lower() for tk in VITAL_TOK)}
        d.update({"rna_tissue_spec": h.get("RNA tissue specificity"),
                  "rna_tissue_dist": h.get("RNA tissue distribution"),
                  "hpa_vital_enriched": ", ".join(f"{k}:{v:.0f}" for k, v in vital.items()) or None,
                  "hpa_specific_tissues": ", ".join(
                      f"{k}:{v:.0f}" for k, v in sorted(spec.items(), key=lambda x: -x[1])[:5]),
                  "rna_sc_spec": h.get("RNA single cell type specificity")})
    else:
        d.update({"rna_tissue_spec": None, "rna_tissue_dist": None,
                  "hpa_vital_enriched": None, "hpa_specific_tissues": None,
                  "rna_sc_spec": None, "hpa_err": h.get("__error__")})
    return d


ev = {s: parse_gene(s) for s in all_syms if ot[s]}
with open(f"{RES}/evidence/evidence_detail.json", "w") as fh:
    json.dump(ev, fh, indent=1, default=str)

# ============================================================================
# STEP 10 — SCORING: drug-blind 4-dimension rubric, applied to all 101
# ============================================================================
# Design principles:
#   * ONE scored track: the expression-derived surface/secreted candidates.
#   * Drug/trial knowledge is HELD OUT of all selection and scoring. It is revealed
#     only in STEP 11-12 as held-out annotation.
#   * Association comes from OUR data (CPTAC protein + RNA), not from Open Targets
#     overall/literature/known-drug scores, which carry publication and drug bias.
WEIGHTS = {"assoc_expr": 0.30, "tractability": 0.25, "safety": 0.20, "selectivity": 0.15}
_s = sum(WEIGHTS.values())
WEIGHTS = {k: v / _s for k, v in WEIGHTS.items()}


def clamp(x, lo=0, hi=5):
    return max(lo, min(hi, x))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# --- association: percentile of protein & RNA overexpression --------------------
# Ranked against the FULL RNA-concordant pool, not against the survivors. Ranking the
# 12 against themselves would re-rank a list already sorted on that same axis, and a
# within-set percentile discards magnitude: the retained candidates span a narrow
# log2FC band that a self-referential rank stretches across the whole 0-5 range.
_p_pct = sig.log2fc.rank(pct=True)
_r_pct = sig.rna_log2fc.rank(pct=True)

# --- intrinsic tractability: antibody-accessibility evidence grade ---------------
# Open Targets AB-tractability buckets with the clinical-status buckets removed, so
# no drug knowledge enters a drug-blind score. The top grade requires two orthogonal
# high-confidence localisation calls (UniProt curated AND GO cellular component).
CLINICAL_BUCKETS = {"Advanced Clinical", "Phase 1 Clinical", "Approved Drug"}


def intrinsic_tract(g):
    ab = set(ev[g].get("tract_AB") or []) - CLINICAL_BUCKETS
    return 3.0 if ("UniProt loc high conf" in ab and "GO CC high conf" in ab) else 1.5


def score_safety(g):
    e = ev[g]
    depmap, loeuf = _f(e.get("depmap_mean_geneeffect")), _f(e.get("loeuf"))
    vital, nliab = e.get("hpa_vital_enriched"), int(e.get("n_safety_liabilities") or 0)
    sc = 3.0
    if depmap is not None:
        sc += -1.5 if depmap <= -1.0 else (-0.7 if depmap <= -0.5 else 0.7)
    if loeuf is not None:
        sc += -0.8 if loeuf < 0.35 else (0.5 if loeuf > 1.0 else 0)
    sc += -1.0 if vital not in (None, "None", "") else 0.5
    if nliab > 0:
        sc -= 0.5 * min(nliab, 2)
    return round(clamp(sc), 2)


def score_selectivity(g):
    e = ev[g]
    npara = int(e.get("n_close_paralogs_ge30pct") or 0)
    maxid = _f(e.get("max_paralog_identity"))
    spec = str(e.get("rna_tissue_spec", ""))
    sc = 3.5
    sc += 0.8 if npara == 0 else (-1.5 if npara >= 5 else (-0.8 if npara >= 2 else 0))
    if maxid and maxid >= 70:
        sc -= 1.0
    if "enriched" in spec.lower():
        sc += 0.8
    elif "Low tissue specificity" in spec:
        sc -= 0.8
    return round(clamp(sc), 2)


def score_frame(genes_df):
    B = genes_df.copy()
    B["assoc_expr"] = (5 * (B.gene.map(_p_pct).values + B.gene.map(_r_pct).values) / 2).round(2)
    B["tractability"] = B.gene.map(intrinsic_tract)
    B["safety"] = B.gene.map(score_safety)
    B["selectivity"] = B.gene.map(score_selectivity)
    B["composite_0_100"] = (20 * sum(WEIGHTS[k] * B[k] for k in WEIGHTS)).round(1)
    B = B.sort_values("composite_0_100", ascending=False).reset_index(drop=True)
    B.insert(0, "rank", B.index + 1)
    return B


scored101 = score_frame(surf_set[surf_set.gene.isin(ev)])
retained12 = set(trackB.gene)
scored101["retained_top12"] = scored101.gene.isin(retained12)
print(f"scored all {len(scored101)} surface candidates")

# ============================================================================
# STEP 11 — HELD-OUT clinical annotation for all 101 (never enters the score)
# ============================================================================
# Positive definitions, fixed before inspecting which genes they flag:
#   any_indication_positive : >=1 agent whose Open Targets mechanism target is this
#                             gene, having reached >= Phase 1 in ANY disease.
#   ccrcc_direct_positive   : >=1 such agent appears as an intervention in a
#                             ClinicalTrials.gov study whose condition is RCC.
def ct_gov(drug):
    key = urllib.parse.quote_plus(drug.lower())[:80]

    def _fetch():
        url = "https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode({
            "query.cond": "renal cell carcinoma", "query.intr": drug,
            "countTotal": "true", "pageSize": 20,
            "fields": "NCTId,BriefTitle,Phase,OverallStatus,StudyType"})
        try:
            return json.loads(urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=45).read())
        except Exception as e:                                # noqa: BLE001
            return {"__error__": f"{type(e).__name__}:{str(e)[:60]}"}

    r = cached(key, _fetch, cache_dir=CACHE_CT)
    time.sleep(0.0 if os.path.exists(os.path.join(CACHE_CT, f"{key}.json")) else 0.25)
    return r


clin_rows = []
for g in scored101.gene:
    e = ev[g]
    agents = e["agents"]
    any_pos = bool(agents) and e["max_clinical_phase"] >= 1
    ncts, hit_agents, interventional = [], [], False
    if any_pos:
        for a in agents:
            d = ct_gov(a)
            if "__error__" in d or not d.get("totalCount") or not d.get("studies"):
                continue
            hit_agents.append(a)
            for s_ in d["studies"]:
                ps = s_["protocolSection"]
                ncts.append(ps["identificationModule"]["nctId"])
                if ps.get("designModule", {}).get("studyType") == "INTERVENTIONAL":
                    interventional = True
    clin_rows.append({
        "gene": g, "n_agents": len(agents), "max_clinical_phase": e["max_clinical_phase"],
        "any_indication_positive": any_pos, "ccrcc_direct_positive": bool(hit_agents),
        "ccrcc_interventional": interventional,
        "agent_modality": "",     # curated layer: ChEMBL misses ADCs; see module docstring
        "agents": "; ".join(agents)[:300], "ccrcc_agents": "; ".join(hit_agents)[:300],
        "example_ncts": ";".join(sorted(set(ncts))[:5])})
clin = pd.DataFrame(clin_rows)

# Preserve any hand-curated agent_modality edits across re-runs.
_cur_path = f"{CURATED}/clinical_agents.csv"
if os.path.exists(_cur_path):
    prev = pd.read_csv(_cur_path).set_index("gene")["agent_modality"].fillna("")
    clin["agent_modality"] = clin.gene.map(prev).fillna("")
clin.to_csv(_cur_path, index=False)

scored101 = scored101.merge(clin, on="gene", how="left")
scored101.to_csv(f"{RES}/scorecards/scorecard_101.csv", index=False)
print(f"held-out positives: {int(clin.ccrcc_direct_positive.sum())} ccRCC-direct, "
      f"{int(clin.any_indication_positive.sum())} any-indication (of {len(clin)})")

# clinical precedent table (held-out; drivers + effectors included for context)
pd.DataFrame([{"symbol": s, "n_drugs": ev[s]["n_drugs"],
               "max_phase": ev[s]["max_clinical_phase"],
               "agents": "; ".join(ev[s]["agents"])[:200],
               "track": ("B_surface_secreted" if s in set(surf_set.gene)
                         else ("A_mutation_driver" if s in DRIVERS else "effector"))}
              for s in sorted(ev)]).to_csv(
    f"{RES}/evidence/clinical_precedent.csv", index=False)

# ============================================================================
# STEP 12 — held-out recovery evaluation
# ============================================================================
N = len(scored101)
recovery = {
    "method": ("single expression-driven discovery; drug/trial data held out of all "
               "selection and scoring"),
    "n_surface_candidates": N,
    "n_retained_top12": int(scored101.retained_top12.sum()),
}

for label, col in [("ccrcc_direct", "ccrcc_direct_positive"),
                   ("any_indication", "any_indication_positive")]:
    pos = scored101[scored101[col]]
    P, n_ret = len(pos), int(scored101.retained_top12.sum())
    k = int(pos.retained_top12.sum())
    ranks = sorted(int(r) for r in pos["rank"])
    U, pmw = mannwhitneyu(pos.composite_0_100,
                          scored101[~scored101[col]].composite_0_100, alternative="greater")
    miss = pos[~pos.retained_top12]
    lost_immune = sorted(miss[miss.gene.isin(
        surf_set[surf_set.immune_infiltration_flag].gene)].gene)
    recovery[label] = {
        "n_positives": P,
        "positives": sorted(pos.gene),
        # (a) does SELECTION concentrate positives into the retained 12?
        "cascade_test": {
            "positives_retained": k, "expected_by_chance": round(P * n_ret / N, 2),
            "hypergeom_p": round(float(hypergeom.sf(k - 1, N, P, n_ret)), 4)},
        # (b) does the COMPOSITE rank positives highly across all candidates?
        "ranking_test": {
            "positive_ranks": ranks,
            "top_k_p": {str(t): round(float(hypergeom.sf(
                sum(1 for r in ranks if r <= t) - 1, N, P, t)), 4) for t in (12, 20, 30)},
            "mannwhitney_U": float(U), "mannwhitney_p": round(float(pmw), 4)},
        "misses": {"n": len(miss),
                   "dropped_by_immune_filter": lost_immune,
                   "dropped_by_fc_cutoff": sorted(set(miss.gene) - set(lost_immune))},
    }

# The legacy statistic, retained for comparability and explicitly labelled as
# conditional on selection: it permutes only the 12 survivors, so it cannot see a
# positive dropped anywhere in the 11,710 -> 12 cascade.
_t12 = scored101[scored101.retained_top12].sort_values("composite_0_100", ascending=False)
_t12 = _t12.assign(rank12=range(1, len(_t12) + 1))
_leg = _t12[_t12.ccrcc_direct_positive]
if len(_leg):
    worst = int(_leg.rank12.max())
    recovery["ranking_given_selection_LEGACY"] = {
        "caveat": ("conditional on selection: permutes only the 12 survivors, so a "
                   "positive dropped earlier in the cascade is invisible. Not the "
                   "headline result; see cascade_test and ranking_test."),
        "ranks_within_12": {r.gene: int(r.rank12) for r in _leg.itertuples()},
        "both_within_top": worst,
        "hypergeom_p": round(float(hypergeom.sf(len(_leg) - 1, len(_t12),
                                                len(_leg), worst)), 4)}

recovery["immune_filter"] = {
    "definition": (f"Spearman rho vs {len(_pres)}-marker leukocyte signature across "
                   f"{prot_t.shape[1]} tumours; exclude rho>=0.40"),
    "n_dropped": n_immune_dropped,
    "CD70_rho": round(float(surf_set.loc[surf_set.gene == "CD70", "immune_rho"].iloc[0]), 3)
    if (surf_set.gene == "CD70").any() else None}
recovery["driver_overlap"] = (f"{len(set(DRIVERS) & set(surf_set.gene))}/{len(DRIVERS)} "
                              "mutation drivers appear in the surface/secreted set")

with open(f"{RES}/scorecards/recovery_stats.json", "w") as fh:
    json.dump(recovery, fh, indent=2)
for lab in ("ccrcc_direct", "any_indication"):
    r = recovery[lab]
    print(f"  {lab}: {r['n_positives']} positives | cascade p="
          f"{r['cascade_test']['hypergeom_p']} | ranking MW p="
          f"{r['ranking_test']['mannwhitney_p']}")

# ============================================================================
# STEP 13 — the selected 12: scorecard + tiers
# ============================================================================
B = scored101[scored101.retained_top12].copy()
B = B.sort_values("composite_0_100", ascending=False).reset_index(drop=True)
B["rank"] = B.index + 1
B["validation_heldout"] = np.where(B.ccrcc_direct_positive, "ccrcc-direct-agent", "none")


def tier(r):
    if r.ccrcc_direct_positive:
        return "T1 fast-follow"
    if r.composite_0_100 >= 70:
        return "T2 discovery"
    if r.composite_0_100 >= 55:
        return "T3 watch"
    return "T4 deprioritize"


B["dev_tier"] = B.apply(tier, axis=1)
B[["rank", "gene", "prot_log2fc", "rna_log2fc", "assoc_expr", "immune_rho", "tractability",
   "safety", "selectivity", "composite_0_100", "validation_heldout", "dev_tier"]].to_csv(
    f"{RES}/scorecards/scorecard_clean.csv", index=False)
print("scorecard (12):", B.gene.tolist())

# ============================================================================
# STEP 14 — mutation drivers as unscored CONTEXT
# ============================================================================
pd.DataFrame([{"symbol": g, "track": "A_mutation_driver", "ensembl_id": ens_map.get(g),
               "uniprot": uni_map.get(g), "mutation_rate_pct": MUTRATE[g],
               "lesion_type": "LoF tumour suppressor",
               "in_surface_target_set": "Yes" if g in set(surf_set.gene) else "No",
               "source": "Clark et al. 2019 Cell, DOI:10.1016/j.cell.2019.10.007"}
              for g in DRIVERS]).to_csv(f"{RES}/scorecards/drivers_context.csv", index=False)

# ============================================================================
# STEP 15 — grade / stage association for the scored candidates
# ============================================================================
# Does the candidate's tumour abundance track histologic grade or pathological stage?
# Purely descriptive context: it is computed AFTER scoring and feeds nothing upstream.
_cli = pd.read_csv(f"{DATA}/HS_CPTAC_CCRCC_CLI.tsi", sep="\t", index_col=0)
_cli = _cli[_cli.index.notna()]
_grade = pd.to_numeric(_cli["Histologic_Grade"].str.extract(r"G(\d)")[0], errors="coerce")
_stage = _cli["Tumor_Stage_Pathological"].map(
    {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})


def _clin_assoc(g, series):
    if g not in prot_t.index:
        return np.nan, np.nan
    c = prot_t.columns.intersection(series.dropna().index)
    if len(c) < 10:
        return np.nan, np.nan
    rho, p = stats.spearmanr(prot_t.loc[g, c], series[c], nan_policy="omit")
    return rho, p


_ga = []
for g in scored101.gene:
    gr, gp = _clin_assoc(g, _grade)
    sr, sp = _clin_assoc(g, _stage)
    _ga.append({"symbol": g, "grade_rho": _f(gr), "grade_p": _f(gp),
                "stage_rho": _f(sr), "stage_p": _f(sp)})
ga = pd.DataFrame(_ga)
for _pre in ("grade", "stage"):
    m = ga[f"{_pre}_p"].notna()
    ga.loc[m, f"{_pre}_fdr"] = multipletests(ga.loc[m, f"{_pre}_p"], method="fdr_bh")[1]
    ga[f"{_pre}_sig"] = ga[f"{_pre}_fdr"] < 0.05
ga = ga.merge(scored101[["gene", "composite_0_100", "retained_top12"]],
              left_on="symbol", right_on="gene", how="left").drop(columns="gene")
ga.round(4).to_csv(f"{RES}/enrichment/grade_stage_association.csv", index=False)
print(f"grade/stage association: {int(ga.grade_sig.sum())} grade-significant, "
      f"{int(ga.stage_sig.sum())} stage-significant (BH-FDR<0.05, n={len(ga)})")

# ============================================================================
# STEP 16 — weight-sensitivity Monte Carlo over the selected 12
# ============================================================================
rng = np.random.default_rng(0)
dims = list(WEIGHTS)
draws = rng.dirichlet(np.array([WEIGHTS[d] for d in dims]) * 20, size=2000)
mat = B[dims].to_numpy()
mc = (draws @ mat.T) * 20
order = (-mc).argsort(axis=1).argsort(axis=1) + 1
pd.DataFrame({"gene": B.gene, "base_composite": B.composite_0_100,
              "mc_mean_rank": order.mean(axis=0).round(2),
              "mc_rank_std": order.std(axis=0).round(2),
              "p_rank1": (order == 1).mean(axis=0).round(3)}).to_csv(
    f"{RES}/enrichment/weight_sensitivity.csv", index=False)
print("weight sensitivity: P(rank1) top =",
      B.gene.iloc[int((order == 1).mean(axis=0).argmax())])

print("\nDONE — all outputs written under results/ and curated/")
