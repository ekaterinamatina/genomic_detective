# hpo_data.py
"""
Слой знаний о фенотипах.
"""
from __future__ import annotations

import csv
import math
import os
import re
import time
import urllib.request
from functools import lru_cache

DATA_DIR = os.environ.get("GD_DATA_DIR", "data")
G2P_URL = "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
HPOA_URL = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
G2P_FILE = os.path.join(DATA_DIR, "genes_to_phenotype.txt")
HPOA_FILE = os.path.join(DATA_DIR, "phenotype.hpoa")
OMIM_FILE = os.path.join(DATA_DIR, "genemap2.txt")
MAX_AGE_DAYS = 90

HPO_ID_RE = re.compile(r"^HP:\d{7}$")

FREQ_TERMS = {
    "HP:0040280": 1.0, "HP:0040281": 0.90, "HP:0040282": 0.55,
    "HP:0040283": 0.17, "HP:0040284": 0.025, "HP:0040285": 0.0,
}

MOI_TERMS = {
    "HP:0000006": "AD", "HP:0000007": "AR", "HP:0032113": "SD",
    "HP:0001417": "XL", "HP:0001419": "XLR", "HP:0001423": "XLD",
    "HP:0001450": "YL", "HP:0001427": "MT",
}

_ont_ready = False

# --------------------------------------------------------------------------
# Онтология (PyHPO)
# --------------------------------------------------------------------------
def init_ontology():
    global _ont_ready
    if not _ont_ready:
        from pyhpo import Ontology
        Ontology()
        _ont_ready = True

@lru_cache(maxsize=None)
def term(hid: str):
    if not isinstance(hid, str) or not HPO_ID_RE.match(hid.strip()):
        return None
    init_ontology()
    from pyhpo import Ontology
    try:
        return Ontology.get_hpo_object(hid.strip())
    except Exception:
        return None

def canonical_id(hid: str):
    t = term(hid)
    return t.id if t is not None else None

def term_name(hid: str) -> str:
    t = term(hid)
    return t.name if t is not None else hid

@lru_cache(maxsize=None)
def ancestors(hid: str) -> frozenset:
    t = term(hid)
    if t is None:
        return frozenset()
    return frozenset({p.id for p in t.all_parents} | {t.id})

def phenotype_options() -> dict:
    init_ontology()
    from pyhpo import Ontology
    root = term("HP:0000118")
    out = {}
    for t in Ontology:
        if t.id == root.id: continue
        try:
            if t.child_of(root):
                out[f"{t.name} ({t.id})"] = t.id
        except Exception: continue
    return out

# --------------------------------------------------------------------------
# Загрузка файлов
# --------------------------------------------------------------------------
def _download(url: str, path: str):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".part"
    req = urllib.request.Request(url, headers={"User-Agent": "genomic-detective/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
        f.write(r.read())
    os.replace(tmp, path)

def ensure_files(warnings: list | None = None):
    for url, path in ((G2P_URL, G2P_FILE), (HPOA_URL, HPOA_FILE)):
        exists = os.path.exists(path)
        stale = exists and (time.time() - os.path.getmtime(path)) > MAX_AGE_DAYS * 86400
        if exists and not stale: continue
        try:
            _download(url, path)
        except Exception as e:
            if not exists:
                raise RuntimeError(f"Не удалось скачать {url}: {e}. Положите файл вручную в {path}.") from e
            if warnings is not None:
                warnings.append(f"Не удалось обновить {os.path.basename(path)} ({e}); используется старая копия.")

def parse_freq(s: str, default: float) -> float:
    s = (s or "").strip()
    if not s: return default
    if s in FREQ_TERMS: return FREQ_TERMS[s]
    m = re.fullmatch(r"(\d+)\s*/\s*(\d+)", s)
    if m and int(m.group(2)) > 0: return int(m.group(1)) / int(m.group(2))
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*%", s)
    if m: return min(float(m.group(1)) / 100.0, 1.0)
    return default

# --------------------------------------------------------------------------
# OMIM
# --------------------------------------------------------------------------
def load_omim_inheritance(genemap_path=OMIM_FILE):
    """Парсит genemap2.txt для получения типов наследования из OMIM."""
    inheritance_map = {"1": "AD", "2": "AR", "3": "AD/AR", "4": "XL", "5": "YL", "6": "MT"}
    gene_inheritance = {}
    if not os.path.exists(genemap_path):
        return gene_inheritance
        
    with open(genemap_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.strip().split("|")
            if len(parts) < 12: continue
            gene_symbols = parts[5]
            phenotypes = parts[11]
            if not phenotypes: continue
            official_gene = gene_symbols.split(",")[0].strip()
            first_phenotype = phenotypes.split(";")[0].strip()
            matches = re.findall(r"\((\d)\)", first_phenotype)
            if matches:
                code = matches[0]
                inh = inheritance_map.get(code, "unknown")
                gene_inheritance[official_gene] = {"inheritance": inh, "phenotype": first_phenotype}
    return gene_inheritance

# --------------------------------------------------------------------------
# База знаний
# --------------------------------------------------------------------------
class KB:
    def __init__(self):
        self.diseases = {}   # disease_id -> {"name", "ann": {hpo: freq}, "moi": set}
        self.gene2dis = {}   # gene symbol -> [disease_id]
        self.version = "unknown"
        self.omim_inh = {}   # 🚨 Новое поле для OMIM

_KB = None

def kb(default_freq: float = 0.5, warnings: list | None = None) -> KB:
    global _KB
    if _KB is not None:
        return _KB
    ensure_files(warnings)
    k = KB()

    # 1) OMIM
    k.omim_inh = load_omim_inheritance()

    # 2) phenotype.hpoa
    with open(HPOA_FILE, encoding="utf-8") as f:
        header = None
        for line in f:
            if line.startswith("#"):
                m = re.match(r"#\s*version:\s*(\S+)", line)
                if m: k.version = m.group(1)
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = {name.strip().lower(): i for i, name in enumerate(parts)}
                need = {"database_id", "disease_name", "hpo_id", "aspect"}
                if not need.issubset(header):
                    raise RuntimeError(f"Неожиданный формат phenotype.hpoa: {parts[:6]}")
                continue

            def col(name):
                i = header.get(name)
                return parts[i].strip() if i is not None and i < len(parts) else ""

            did, hid, aspect = col("database_id"), col("hpo_id"), col("aspect")
            d = k.diseases.setdefault(did, {"name": col("disease_name"), "ann": {}, "moi": set()})
            if aspect == "I":
                code = MOI_TERMS.get(hid)
                if code: d["moi"].add(code)
            elif aspect == "P":
                if col("qualifier").upper() == "NOT": continue
                fr = parse_freq(col("frequency"), default_freq)
                if fr > d["ann"].get(hid, -1): d["ann"][hid] = fr

    # 3) genes_to_phenotype.txt
    with open(G2P_FILE, encoding="utf-8") as f:
        first = f.readline().rstrip("\n")
        cols = [c.strip().lower().lstrip("#") for c in first.split("\t")]
        has_new = "gene_symbol" in cols and "disease_id" in cols
        gi = cols.index("gene_symbol") if has_new else None
        di = cols.index("disease_id") if has_new else None
        seen = set()
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if has_new:
                if len(parts) <= max(gi, di): continue
                gene, did = parts[gi].strip(), parts[di].strip()
                if gene and did in k.diseases and (gene, did) not in seen:
                    seen.add((gene, did))
                    k.gene2dis.setdefault(gene, []).append(did)
        if not has_new:
            raise RuntimeError("genes_to_phenotype.txt в старом формате. Скачайте актуальную версию: " + G2P_URL)
    _KB = k
    return k

def _p_disease(ann: dict, hid: str, discount: float):
    anc_q = ancestors(hid)
    best = None
    for a, fr in ann.items():
        if hid in ancestors(a):
            best = fr if best is None else max(best, fr)
    if best is not None: return best, "matched", None
    via = None
    for a, fr in ann.items():
        if a in anc_q:
            if best is None or fr * discount > best:
                best, via = fr * discount, a
    if best is not None: return best, "ancestor", via
    return None, "unexplained", None

@lru_cache(maxsize=None)
def background_prob(hid: str) -> float:
    k = kb()
    dids = {d for ds in k.gene2dis.values() for d in ds}
    n = 0
    for did in dids:
        if any(hid in ancestors(a) for a in k.diseases[did]["ann"]):
            n += 1
    return (n + 0.5) / (len(dids) + 1.0)

def phenotype_lr(disease: dict, present, excluded, cfg: dict):
    total, details = 0.0, []
    for t in present:
        p, kind, via = _p_disease(disease["ann"], t, cfg["ancestor_discount"])
        pb = background_prob(t)
        if kind == "unexplained":
            lr = cfg["unexplained_lr"]
        else:
            lr = min(max(p / pb, 1.0 / cfg["term_lr_cap"]), cfg["term_lr_cap"])
        total += math.log10(lr)
        details.append({"term": t, "name": term_name(t), "kind": kind, "p_disease": p,
                        "p_background": round(pb, 5), "log10_lr": round(math.log10(lr), 3),
                        "via": via})
    for e in excluded:
        best = None
        for a, fr in disease["ann"].items():
            if e in ancestors(a): best = fr if best is None else max(best, fr)
        if best is None:
            lr = 1.0
            kind = "excluded_not_expected"
        else:
            lr = max(1.0 - best, cfg["excluded_lr_floor"])
            kind = "excluded_conflict"
        total += math.log10(lr)
        details.append({"term": e, "name": term_name(e), "kind": kind, "p_disease": best,
                        "p_background": None, "log10_lr": round(math.log10(lr), 3)})
    return total, details