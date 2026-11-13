# core_pipeline.py
from __future__ import annotations
import json, math, os, re, time, urllib.error, urllib.request
from collections import defaultdict
import pandas as pd
import hpo_data as H

CACHE_DIR = os.environ.get("GD_CACHE_DIR", "cache")

CFG = dict(
    default_freq=0.5, ancestor_discount=0.25, unexplained_lr=0.2, term_lr_cap=200.0, excluded_lr_floor=0.02,
    single_het_ar=0.1, comphet_unphased=0.5, xlr_female_het=0.2, xl_unspecified_het=0.5, moi_mismatch=0.01, known_moi_mismatch=0.05,
    af_thr_dominant=0.001, af_thr_recessive=0.02, af_thr_unknown=0.01, min_an=2000,
    min_dp=8, min_gq=20, top_n=10, tier_strong=6.0, tier_probable=4.0, tier_possible=2.0,
    competing_margin=1.0, generic_bg=0.2, no_match_cap=1.0, max_variants_per_gene=8, max_variants_total=30000,
    vep_batch=100, gnomad_delay=2.0, gnomad_max_queries=60,
)

LR_STRONG, LR_MODERATE, LR_SUPPORTING = 18.7, 4.33, 2.08
def _lr_needed(prior: float, post: float) -> float: return (post / (1 - post)) / (prior / (1 - prior))
CLINVAR_LR = {"pathogenic": _lr_needed(0.10, 0.99), "likely_pathogenic": _lr_needed(0.10, 0.90)}
REVEL_PATH = [(0.932, LR_STRONG), (0.773, LR_MODERATE), (0.644, LR_SUPPORTING)]
REVEL_BENIGN = [(0.016, 1 / LR_STRONG), (0.183, 1 / LR_MODERATE), (0.290, 1 / LR_SUPPORTING)]
LOF_TERMS = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost"}
SUPPORT_TERMS = {"inframe_insertion", "inframe_deletion", "stop_lost", "protein_altering_variant"}
QUALIFYING = LOF_TERMS | SUPPORT_TERMS | {"missense_variant"}
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

ACMG_SF_V2 = {"ACTA2", "ACTC1", "APC", "APOB", "ATP7B", "BMPR1A", "BRCA1", "BRCA2", "CACNA1S", "COL3A1", "DSC2", "DSG2", "DSP", "FBN1", "GLA", "KCNH2", "KCNQ1", "LDLR", "LMNA", "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH", "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3", "NF2", "OTC", "PCSK9", "PKP2", "PMS2", "PRKAG2", "PTEN", "RB1", "RET", "RYR1", "RYR2", "SCN5A", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD3", "SMAD4", "STK11", "TGFBR1", "TGFBR2", "TMEM43", "TNNI3", "TNNT2", "TP53", "TPM1", "TSC1", "TSC2", "VHL", "WT1"}

POPULATIONS = {"nfe", "afr", "eas", "sas", "amr", "asj", "fin", "mid"}
VEP_SERVERS = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
GNOMAD_DATASET = {"GRCh38": "gnomad_r4", "GRCh37": "gnomad_r2_1"}
VEP_EXTRAS = {"hgvs": 1, "canonical": 1, "mane": 1, "REVEL": 1}
CLINVAR_RU = {"pathogenic": "патогенный", "likely_pathogenic": "вероятно патогенный", "uncertain": "VUS", "conflicting": "противоречивые", "benign": "доброкачественный", "likely_benign": "вероятно доброкачественный", "other": "прочее", "none": "нет записи"}
ZYG_RU = {"het": "гетерозигота", "hom": "гомозигота", "hemi": "гемизигота"}

class AnnotationError(RuntimeError): pass
_CACHES = {}

def _cache(name):
    if name not in _CACHES:
        path = os.path.join(CACHE_DIR, name)
        data = {}
        if os.path.exists(path):
            try:
                with open(path, encoding="utf-8") as f: data = json.load(f)
            except Exception: data = {}
        _CACHES[name] = (path, data)
    return _CACHES[name]

def _save_cache(name):
    path, data = _CACHES[name]
    os.makedirs(CACHE_DIR, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f: json.dump(data, f)
    os.replace(tmp, path)

def norm_chrom(c) -> str:
    c = str(c).strip()
    if c.lower().startswith("chr"): c = c[3:]
    return "MT" if c in ("M", "MT", "m", "mt") else c

def parse_gt(gt) -> str:
    s = str(gt).strip().lower()
    simple = {"het": "het", "hom": "hom", "hemi": "hemi", "heterozygous": "het", "homozygous": "hom"}
    if s in simple: return simple[s]
    if s in ("", ".", "./.", ".|.", "nan", "none"): return "missing"
    alleles = re.split(r"[/|]", s)
    if any(a in (".", "") for a in alleles): return "missing"
    if all(a == "0" for a in alleles): return "ref"
    if len(alleles) == 1: return "hemi"
    return "hom" if len(set(alleles)) == 1 else "het"

def parse_phenopacket(path_or_dict):
    if isinstance(path_or_dict, dict): data = path_or_dict
    else:
        with open(path_or_dict, encoding="utf-8") as f: data = json.load(f)
    if "phenopacket" in data and isinstance(data["phenopacket"], dict): data = data["phenopacket"]
    present, excluded, unknown = [], [], []
    for ft in data.get("phenotypicFeatures", []) or []:
        raw = (ft.get("type") or {}).get("id", "")
        cid = H.canonical_id(raw)
        if cid is None:
            unknown.append(raw)
            continue
        (excluded if ft.get("excluded") else present).append(cid)
    present, excluded = list(dict.fromkeys(present)), list(dict.fromkeys(excluded))
    conflicts = sorted(set(present) & set(excluded))
    present = [t for t in present if t not in conflicts]
    excluded = [t for t in excluded if t not in conflicts]
    sex_raw = str((data.get("subject") or {}).get("sex", "")).upper()
    sex = {"MALE": "male", "FEMALE": "female"}.get(sex_raw, "unknown")
    return dict(present=present, excluded=excluded, unknown=unknown, conflicts=conflicts, sex=sex)

def normalize_table(df: pd.DataFrame, warnings: list) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    for col in ("chrom", "pos", "ref", "alt"):
        if col not in df.columns: raise ValueError(f"Ошибка валидации: нет колонки '{col}'")
    if df[["chrom", "pos", "ref", "alt"]].isnull().any().any(): raise ValueError("Ошибка: пустые координаты.")
    if "genotype" not in df.columns:
        warnings.append("Нет колонки genotype: все варианты считаются гетерозиготными.")
        df["genotype"] = "0/1"
    df["chrom"] = df["chrom"].map(norm_chrom)
    df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(int)
    df["ref"] = df["ref"].astype(str).str.upper().str.strip()
    df["alt"] = df["alt"].astype(str).str.upper().str.strip()
    bad = df["alt"].str.contains(r"[,<>\[\]\*]", regex=True) | df["alt"].isin([".", ""])
    if bad.any():
        warnings.append(f"Пропущено {int(bad.sum())} строк с мультиаллельным ALT.")
        df = df[~bad]
    return df

class _ParamError(Exception): pass

def _vep_post(server, payload, extras):
    body = dict(payload)
    body.update(extras)
    data = json.dumps(body).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(server + "/vep/human/region", data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=90) as r: return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}"
            if e.code == 429: time.sleep(float(e.headers.get("Retry-After", 2 ** attempt)))
            elif e.code == 400: raise _ParamError(e.read().decode("utf-8", "replace")[:300])
            elif e.code >= 500: time.sleep(2 ** attempt)
            else: raise AnnotationError(f"VEP: {last}")
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = str(e)
            time.sleep(2 ** attempt)
    raise AnnotationError(f"VEP не ответил после 5 попыток ({last})")

def classify_clinvar(sigs) -> str:
    toks = set()
    for s in sigs:
        for t in re.split(r"[,/;]", str(s).lower().replace(" ", "_")):
            if t: toks.add(t.strip())
    if not toks: return "none"
    if any(t.startswith("conflicting") for t in toks): return "conflicting"
    p, lp = "pathogenic" in toks, "likely_pathogenic" in toks
    b, lb = "benign" in toks, "likely_benign" in toks
    if (p or lp) and (b or lb): return "conflicting"
    if p: return "pathogenic"
    if lp: return "likely_pathogenic"
    if b: return "benign"
    if lb: return "likely_benign"
    if "uncertain_significance" in toks: return "uncertain"
    return "other"

def extract_clinvar(res: dict, vep_alt: str):
    sigs, ambiguous = [], False
    for cv in res.get("colocated_variants") or []:
        csa = cv.get("clin_sig_allele")
        if csa:
            for part in str(csa).split(";"):
                allele, _, sig = part.partition(":")
                if allele.strip() == vep_alt: sigs.append(sig)
        elif cv.get("clin_sig"):
            alts = [a for a in str(cv.get("allele_string", "")).split("/")[1:]]
            if alts == [vep_alt]: sigs.extend(cv["clin_sig"] if isinstance(cv["clin_sig"], list) else [cv["clin_sig"]])
            else: ambiguous = True
    return classify_clinvar(sigs), ambiguous

def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def parse_vep_result(res: dict) -> list:
    best = {}
    for tc in res.get("transcript_consequences") or []:
        gene = (tc.get("gene_symbol") or "").strip()
        terms = tc.get("consequence_terms") or []
        if not gene or not any(t in QUALIFYING for t in terms): continue
        rank = (1 if tc.get("mane_select") or tc.get("mane") else 0, 1 if tc.get("canonical") else 0, IMPACT_RANK.get(tc.get("impact"), 0))
        if gene not in best or rank > best[gene][0]: best[gene] = (rank, tc)
    out = []
    for gene, (_, tc) in best.items():
        vep_alt = tc.get("variant_allele", "")
        clin, amb = extract_clinvar(res, vep_alt)
        out.append(dict(gene=gene, consequence=[t for t in tc.get("consequence_terms", []) if t in QUALIFYING] or tc.get("consequence_terms", []),
            all_terms=tc.get("consequence_terms", []), impact=tc.get("impact", ""), transcript=tc.get("transcript_id", ""),
            hgvsc=tc.get("hgvsc", ""), hgvsp=tc.get("hgvsp", ""), clinvar=clin, clinvar_ambiguous=amb,
            revel=_f(tc.get("revel")), sift=tc.get("sift_prediction", ""), polyphen=tc.get("polyphen_prediction", "")))
    return out

def vkey(chrom, pos, ref, alt) -> str: return f"{chrom}-{pos}-{ref}-{alt}"

def annotate_with_vep(variants, build, warnings, progress=None):
    cname = f"vep_{build}.json"
    path, cache = _cache(cname)
    server = VEP_SERVERS[build]
    todo = [v for v in variants if vkey(*v) not in cache]
    extras = dict(VEP_EXTRAS)
    failed = 0
    for i in range(0, len(todo), CFG["vep_batch"]):
        batch = todo[i:i + CFG["vep_batch"]]
        payload = {"variants": [f"{c} {p} . {r} {a} . . ." for c, p, r, a in batch]}
        try:
            try: results = _vep_post(server, payload, extras)
            except _ParamError as e:
                if not extras: raise AnnotationError(f"VEP отклонил запрос: {e}")
                extras = {}
                warnings.append("Сервер VEP не принял расширенные параметры (REVEL): аннотация без REVEL.")
                results = _vep_post(server, payload, extras)
        except AnnotationError as e:
            failed += len(batch)
            warnings.append(f"Не аннотировано {len(batch)} вариантов: {e}")
            continue
        seen = set()
        for res in results:
            p = str(res.get("input", "")).split()
            if len(p) < 5: continue
            k = vkey(p[0], p[1], p[3], p[4])
            cache[k] = parse_vep_result(res)
            seen.add(k)
        for v in batch:
            if vkey(*v) not in seen: failed += 1
        if progress: progress(min(i + len(batch), len(todo)), len(todo))
        time.sleep(0.2)
    if todo: _save_cache(cname)
    if failed: warnings.append(f"VEP не вернул аннотацию для {failed} вариантов.")
    return {vkey(*v): cache.get(vkey(*v)) for v in variants}, failed

def gnomad_id(v) -> str: return f"{'M' if v['chrom'] == 'MT' else v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}"

def summarize_gnomad(vdata) -> dict:
    min_an = CFG["min_an"]
    tot_ac = tot_an = 0
    pops = defaultdict(lambda: [0, 0])
    for src in ("exome", "genome"):
        s = (vdata or {}).get(src)
        if not s: continue
        tot_ac += s.get("ac") or 0
        tot_an += s.get("an") or 0
        for p in s.get("populations") or []:
            pid = p.get("id", "")
            if "_" in pid or pid in ("XX", "XY", "remaining"): continue
            pops[pid][0] += p.get("ac") or 0
            pops[pid][1] += p.get("an") or 0
    if tot_ac == 0: return dict(status="absent")
    group_af = {pid: ac / an for pid, (ac, an) in pops.items() if an >= min_an}
    af_global = tot_ac / tot_an if tot_an else None
    if group_af:
        gpop = max(group_af, key=group_af.get)
        af_filter, gname = group_af[gpop], gpop
    else: af_filter, gname = af_global, "global"
    return dict(status="found", ac=tot_ac, an=tot_an, af_global=af_global, af_filter=af_filter, af_filter_group=gname, group_af=group_af)

def fetch_gnomad(v, build) -> dict:
    cname = f"gnomad_{GNOMAD_DATASET[build]}.json"
    path, cache = _cache(cname)
    vid = gnomad_id(v)
    if vid in cache: return cache[vid]
    query = """query Variant($id: String!, $ds: DatasetId!) { variant(variantId: $id, dataset: $ds) { exome { ac an populations { id ac an } } genome { ac an populations { id ac an } } } }"""
    payload = json.dumps({"query": query, "variables": {"id": vid, "ds": GNOMAD_DATASET[build]}}).encode()
    info = dict(status="error", message="не запрашивалось")
    for attempt in range(4):
        try:
            req = urllib.request.Request("https://gnomad.broadinstitute.org/api", data=payload, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r: res = json.loads(r.read().decode())
            errs = res.get("errors") or []
            vdata = (res.get("data") or {}).get("variant")
            if vdata: info = summarize_gnomad(vdata)
            elif errs and not any("not found" in str(e.get("message", "")).lower() for e in errs): info = dict(status="error", message=str(errs[0].get("message", ""))[:200])
            else: info = dict(status="absent")
            break
        except urllib.error.HTTPError as e:
            info = dict(status="error", message=f"HTTP {e.code}")
            if e.code == 429: time.sleep(5 * 2 ** attempt); continue
            break
        except Exception as e: info = dict(status="error", message=str(e)[:200]); time.sleep(2 ** attempt)
    if info["status"] in ("found", "absent"):
        cache[vid] = info
        _save_cache(cname)
    return info

def variant_base_lr(v):
    parts = []
    cv = v["clinvar"]
    if cv in CLINVAR_LR:
        lr = CLINVAR_LR[cv]
        parts.append((f"ClinVar: {CLINVAR_RU[cv]} (Tavtigian 2018)", lr))
        return lr, parts, "clinvar"
    terms = set(v["all_terms"])
    lr = 1.0
    if terms & LOF_TERMS:
        parts.append(("Потеря функции (VEP)", LR_MODERATE))
        lr *= LR_MODERATE
    elif terms & SUPPORT_TERMS:
        parts.append(("Вставка/делеция без сдвига рамки", LR_SUPPORTING))
        lr *= LR_SUPPORTING
    if "missense_variant" in terms and v["revel"] is not None:
        r = v["revel"]
        for thr, val in REVEL_PATH:
            if r >= thr: parts.append((f"REVEL={r:.3f} (Pejaver 2022)", val)); lr *= val; break
        else:
            for thr, val in REVEL_BENIGN:
                if r <= thr: parts.append((f"REVEL={r:.3f} в пользу доброкачественности", val)); lr *= val; break
    return lr, parts, "computational"

def prepare_case(df: pd.DataFrame, build: str = "GRCh38", progress=None) -> dict:
    if build not in VEP_SERVERS: raise ValueError(f"Сборка {build} не поддерживается.")
    warnings, stats = [], {}
    df = normalize_table(df, warnings)
    stats["n_input"] = len(df)
    if "filter" in df.columns:
        m = ~df["filter"].astype(str).isin(["PASS", ".", "nan", ""])
        stats["failed_filter"] = int(m.sum()); df = df[~m]
    for col, thr, name in (("dp", CFG["min_dp"], "DP"), ("gq", CFG["min_gq"], "GQ")):
        if col in df.columns:
            vals = pd.to_numeric(df[col], errors="coerce")
            m = vals < thr
            stats[f"low_{col}"] = int(m.sum()); df = df[~m.fillna(False)]
    df = df.assign(zyg=df["genotype"].map(parse_gt))
    stats["not_carried"] = int(df["zyg"].isin(["ref", "missing"]).sum())
    df = df[~df["zyg"].isin(["ref", "missing"])]
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    stats["n_after_qc"] = len(df)
    if len(df) > CFG["max_variants_total"]: raise ValueError("Слишком много вариантов для REST VEP.")
    if df.empty: raise ValueError("После фильтров не осталось вариантов.")
    items = [(r.chrom, r.pos, r.ref, r.alt) for r in df.itertuples()]
    ann, n_failed = annotate_with_vep(items, build, warnings, progress)
    stats["vep_failed"] = n_failed
    genes = defaultdict(list)
    n_benign = n_noncoding = 0
    for r in df.itertuples():
        k = vkey(r.chrom, r.pos, r.ref, r.alt)
        a = ann.get(k)
        if a is None: continue
        if not a: n_noncoding += 1; continue
        for g in a:
            if g["clinvar"] in ("benign", "likely_benign"): n_benign += 1; continue
            v = dict(g, key=k, chrom=r.chrom, pos=int(r.pos), ref=r.ref, alt=r.alt, zyg=r.zyg)
            v["lr0"], v["lr0_parts"], v["lr0_source"] = variant_base_lr(v)
            genes[g["gene"]].append(v)
    stats["no_qualifying_consequence"] = n_noncoding
    stats["clinvar_benign_removed"] = n_benign
    stats["n_genes"] = len(genes)
    stats["n_variants_kept"] = sum(len(x) for x in genes.values())
    if stats["no_qualifying_consequence"] > 0: warnings.append("Оцениваются только кодирующие варианты со значимыми последствиями.")
    return dict(build=build, genes=dict(genes), stats=stats, warnings=warnings, afc={})

def af_threshold(moi):
    if moi in ("AR", "XLR"): return CFG["af_thr_recessive"]
    if moi is None: return CFG["af_thr_unknown"]
    return CFG["af_thr_dominant"]

def eff_zyg(v, ctx):
    if v["chrom"] == "X" and ctx["sex"] == "male" and v["zyg"] in ("hom", "hemi"): return "hemi"
    return v["zyg"]

KNOWN_MOI_COMPAT = {"AD": {"AD", "SD"}, "AR": {"AR"}, "XL": {"XL", "XLR", "XLD"}, "de_novo": {"AD", "SD", "XL", "XLD", "XLR"}}

def eval_moi(moi, vs, ctx, afc):
    cfg = CFG
    notes, thr = [], af_threshold(moi)
    usable = []
    for v in vs:
        info = afc.get(v["key"])
        if info and info.get("status") == "found" and info.get("af_filter") is not None and info["af_filter"] > thr:
            notes.append(f"{v['key']} исключён: AF = {info['af_filter']:.4g} > {thr} для наследования {moi or 'не указано'}"); continue
        usable.append(v)
    if not usable: return dict(filtered=True, used=[], log10_moi=0.0, log10_var=0.0, notes=notes)
    chrom = usable[0]["chrom"]
    auto = chrom not in ("X", "Y", "MT")
    log_moi, used, weight = 0.0, [usable[0]], 1
    ez = [(v, eff_zyg(v, ctx)) for v in usable]
    def mismatch(why):
        nonlocal log_moi
        log_moi += math.log10(cfg["moi_mismatch"]); notes.append(why)
    if moi in ("AD", "SD"):
        if not auto: mismatch(f"Аутосомно-доминантная болезнь, а ген на хромосоме {chrom}")
    elif moi == "AR":
        if not auto: mismatch(f"Аутосомно-рецессивная болезнь, а ген на хромосоме {chrom}")
        else:
            homs = [v for v, z in ez if z in ("hom", "hemi")]
            if homs: used, weight = [homs[0]], 2; notes.append("Гомозигота - соответствует рецессивному наследованию")
            elif len(ez) >= 2:
                used = [ez[0][0], ez[1][0]]; log_moi += math.log10(cfg["comphet_unphased"]); notes.append("Два гетерозиготных варианта: возможна компаунд-гетерозигота (фаза неизвестна)")
            else: log_moi += math.log10(cfg["single_het_ar"]); notes.append("Один гетерозиготный вариант при рецессивной болезни: пациент может быть только носителем")
    elif moi in ("XLR", "XLD", "XL"):
        if chrom != "X": mismatch(f"X-сцепленная болезнь, а ген на хромосоме {chrom}")
        else:
            v, z = ez[0]
            if z in ("hemi", "hom"): notes.append("Гемизигота/гомозигота на X - соответствует X-сцепленному наследованию")
            elif ctx["sex"] == "male": mismatch("Гетерозиготный вызов у мужчины на X")
            elif moi == "XLR": log_moi += math.log10(cfg["xlr_female_het"]); notes.append("Гетерозиготная женщина при X-рецессивной болезни: возможна лишь манифестирующая носительница")
            elif moi == "XL": log_moi += math.log10(cfg["xl_unspecified_het"]); notes.append("X-сцепленное наследование без уточнения")
    elif moi == "YL":
        if chrom != "Y" or ctx["sex"] == "female": mismatch("Y-сцепленное наследование несовместимо с хромосомой/полом")
    elif moi == "MT":
        if chrom != "MT": mismatch("Митохондриальное наследование, а вариант не в mtDNA")
        else: notes.append("Гетероплазмия в этом режиме не оценивается")
    else: notes.append("Тип наследования болезни в HPOA не указан - совместимость не проверялась")
    km = ctx["known_moi"]
    if km != "unknown" and moi is not None and moi not in KNOWN_MOI_COMPAT[km]:
        log_moi += math.log10(cfg["known_moi_mismatch"]); notes.append(f"Заданное пользователем наследование ({km}) противоречит типу наследования болезни ({moi})")
    log_var = 0.0
    for v in used:
        lr = v["lr0"]
        info = afc.get(v["key"])
        if v["lr0_source"] == "computational" and info and info.get("status") == "absent": lr *= LR_SUPPORTING
        log_var += weight * math.log10(lr)
    return dict(filtered=False, used=used, log10_moi=log_moi, log10_var=log_var, notes=notes)

def evaluate_gene(gene, vs, ctx, afc):
    k = H.kb(CFG["default_freq"])
    vs = sorted(vs, key=lambda v: -v["lr0"])[:CFG["max_variants_per_gene"]]
    dids = k.gene2dis.get(gene, [])
    options = [(d, k.diseases[d]) for d in dids] or [(None, None)]
    best, all_filtered = None, True
    for did, d in options:
        pl, pdet = (0.0, []) if d is None else H.phenotype_lr(d, ctx["present"], ctx["excluded"], CFG)
        mois = (d["moi"] if d and d["moi"] else {None})
        for moi in sorted(mois, key=str):
            r = eval_moi(moi, vs, ctx, afc)
            if r["filtered"]: continue
            all_filtered = False
            total = pl + r["log10_moi"] + r["log10_var"]
            explains = (d is None or not ctx["present"] or any(x["kind"] == "matched" or (x["kind"] == "ancestor" and H.background_prob(x["via"]) <= CFG["generic_bg"]) for x in pdet))
            capped = (not explains) and total > CFG["no_match_cap"]
            if capped: total = CFG["no_match_cap"]
            if best is None or total > best["total"]: best = dict(total=total, did=did, d=d, moi=moi, pl=pl, pdet=pdet, r=r, capped=capped)
    if best is None: return dict(gene=gene, filtered=True, considered=vs, af_checked=True)
    return dict(gene=gene, filtered=False, considered=vs, best=best, af_checked=all(v["key"] in afc and afc[v["key"]]["status"] != "error" for v in vs) or False)

def _vdesc(v): return f"chr{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}"

def build_candidate(ev, ctx, build, population):
    b, r = ev["best"], ev["best"]["r"]
    used, afc_all, pdet = r["used"], ctx["afc"], b["pdet"]
    db, pred, infer = [], [], []
    for v in used:
        db.append(f"ClinVar (по данным Ensembl VEP) для {_vdesc(v)}: {CLINVAR_RU[v['clinvar']]}" + (" [есть другие аллели]" if v["clinvar_ambiguous"] else ""))
        info = afc_all.get(v["key"])
        if not info: db.append(f"gnomAD: частота для {_vdesc(v)} не проверялась")
        elif info["status"] == "found":
            pa = info["group_af"].get(population)
            db.append(f"gnomAD {GNOMAD_DATASET[build]}: {_vdesc(v)} глобальная AF={info['af_global']:.3g}, максимум по группам {info['af_filter_group']}={info['af_filter']:.3g}" + (f", выбранная популяция {population}={pa:.3g}" if pa is not None else ""))
        elif info["status"] == "absent": db.append(f"gnomAD {GNOMAD_DATASET[build]}: {_vdesc(v)} не найден")
        else: db.append(f"gnomAD: запрос для {_vdesc(v)} не удался ({info.get('message', '')}) - частота НЕИЗВЕСТНА")
    if b["d"] is not None: 
        db.append(f"HPOA: {b['did']} {b['d']['name']}; наследование: {', '.join(sorted(b['d']['moi'])) or 'не указано'}; аннотировано признаков: {len(b['d']['ann'])}")
    else: 
        db.append("HPOA: для гена нет аннотированных болезней")
        
    # 🚨 ДОБАВЛЯЕМ OMIM ИНФОРМАЦИЮ
    k = H.kb()
    if ev["gene"] in k.omim_inh:
        omim_data = k.omim_inh[ev["gene"]]
        db.append(f"OMIM (genemap2): наследование {omim_data['inheritance']}; фенотип: {omim_data['phenotype']}")
    # ---------------------------------
    for v in used:
        s = f"VEP: {', '.join(v['consequence'])} ({v['impact']}), транскрипт {v['transcript']}"
        if v["hgvsc"]: s += f", {v['hgvsc']}"
        if v["hgvsp"]: s += f", {v['hgvsp']}"
        pred.append(s)
        if v["revel"] is not None: pred.append(f"REVEL = {v['revel']:.3f}")
        if v["sift"] or v["polyphen"]: pred.append(f"SIFT: {v['sift'] or 'н/д'}; PolyPhen: {v['polyphen'] or 'н/д'}")
    matched = [x for x in pdet if x["kind"] == "matched"]
    anc = [x for x in pdet if x["kind"] == "ancestor"]
    unexpl = [x for x in pdet if x["kind"] == "unexplained"]
    confl = [x for x in pdet if x["kind"] == "excluded_conflict"]
    infer.append(f"Фенотип: log10 LR = {b['pl']:+.2f}. Совпало {len(matched)} из {len(ctx['present'])} признаков" + (f" ({', '.join(x['name'] for x in matched)})" if matched else "") + (f"; только по общему предку: {', '.join(x['name'] for x in anc)}" if anc else "") + (f"; не объяснены: {', '.join(x['name'] for x in unexpl)}" if unexpl else ""))
    if confl: infer.append("Исключённые у пациента признаки, типичные для болезни (снижают балл): " + ", ".join(f"{x['name']} (частота {x['p_disease']:.0%})" for x in confl))
    infer.append(f"Наследование: log10 LR = {r['log10_moi']:+.2f}")
    infer.extend(r["notes"])
    for v in used:
        for txt, lr in v["lr0_parts"]: infer.append(f"Вариант {_vdesc(v)}: {txt}; LR={lr:.3g}")
        info = afc_all.get(v["key"])
        if v["lr0_source"] == "computational" and info and info.get("status") == "absent": infer.append(f"Вариант {_vdesc(v)} отсутствует в gnomAD: +PM2 (supporting), LR={LR_SUPPORTING}")
    if b.get("capped"): infer.append(f"Болезнь гена не объясняет ни одного наблюдаемого признака: балл ограничен {CFG['no_match_cap']} (вероятно, вторичная находка).")
    infer.append("Балл = log10 LR (шкала силы свидетельств, не вероятность); пороги уверенности эвристические.")
    score = b["total"]
    af_info = afc_all.get(used[0]["key"])
    return dict(gene=ev["gene"], variant=_vdesc(used[0]), variants="; ".join(_vdesc(v) for v in used), zygosity=", ".join(ZYG_RU.get(eff_zyg(v, ctx), v["zyg"]) for v in used), consequence=",".join(used[0]["consequence"]), clinvar=used[0]["clinvar"], disease=(b["d"]["name"] if b["d"] else "нет известной болезни"), disease_id=b["did"] or "", inheritance=("/".join(sorted(b["d"]["moi"])) if b["d"] and b["d"]["moi"] else "Unknown"), score=round(score, 2), pheno_log10=round(b["pl"], 2), moi_log10=round(r["log10_moi"], 2), variant_log10=round(r["log10_var"], 2), af=(af_info.get("af_filter") if af_info and af_info["status"] == "found" else None), af_status=(af_info["status"] if af_info else "not_checked"), n_matched=len(matched), n_unexplained=len(unexpl), evidence_db=db, evidence_pred=pred, evidence_infer=infer)

def tier_of(score):
    if score >= CFG["tier_strong"]: return "сильное свидетельство"
    if score >= CFG["tier_probable"]: return "вероятный кандидат"
    if score >= CFG["tier_possible"]: return "возможный кандидат"
    return "недостаточно данных"

def next_steps(c, ctx):
    s = ["Подтвердить вариант независимым методом (Сэнгер)."]
    if c["af_status"] != "found" and c["af_status"] != "absent": s.append("Проверить популяционную частоту варианта (gnomAD).")
    if len(ctx["present"]) < 3: s.append("Дополнить HPO-описание: признаков меньше трёх.")
    if c["n_unexplained"]: s.append("Часть признаков не объясняется этой болезнью: рассмотреть второй диагноз.")
    if c["inheritance"].startswith("AR") and "гетерозигота" in c["zygosity"] and ";" not in c["variants"]: s.append("Искать второй аллель: CNV, глубокие интронные варианты.")
    if ";" in c["variants"]: s.append("Секвенировать родителей (фазирование компаунд-гетерозиготы).")
    elif ctx["known_moi"] == "unknown": s.append("Данные родителей/семьи (сегрегация) повысят надёжность.")
    if c["clinvar"] in ("uncertain", "none", "conflicting", "other"): s.append("Вариант без однозначной классификации: функциональные данные, сегрегация.")
    if "stop_gained" in c["consequence"] or "frameshift_variant" in c["consequence"]: s.append("Проверить NMD-escape и механизм болезни.")
    return s

def rank_case(case, pheno, sex="unknown", known_moi="unknown", population="nfe", progress=None, acmg_genes=None):
    if known_moi not in ("unknown", "AD", "AR", "XL", "de_novo"): raise ValueError(f"Неизвестный тип наследования: {known_moi}")
    if population not in POPULATIONS: raise ValueError(f"Неизвестная популяция gnomAD: {population}")
    if sex == "unknown" and pheno.get("sex") in ("male", "female"): sex = pheno["sex"]
    warnings = list(case["warnings"])
    if pheno["unknown"]: warnings.append("Не распознаны HPO-ID (пропущены): " + ", ".join(pheno["unknown"]))
    if pheno["conflicts"]: warnings.append("Термины одновременно в 'есть' и 'исключено': " + ", ".join(pheno["conflicts"]))
    if not pheno["present"]: raise ValueError("Нужен хотя бы один распознанный HPO-признак.")
    ctx = dict(present=pheno["present"], excluded=pheno["excluded"], sex=sex, known_moi=known_moi, afc=case["afc"])
    build = case["build"]
    H.kb(CFG["default_freq"], warnings)
    cands = {g: evaluate_gene(g, vs, ctx, case["afc"]) for g, vs in case["genes"].items()}
    queries, budget = 0, CFG["gnomad_max_queries"]
    while True:
        live = sorted((c for c in cands.values() if not c["filtered"]), key=lambda c: -c["best"]["total"])
        target = [c for c in live[:CFG["top_n"]] if not c["af_checked"]]
        if not target or queries >= budget: break
        c = target[0]
        for v in c["considered"]:
            if v["key"] not in case["afc"] or case["afc"][v["key"]]["status"] == "error":
                case["afc"][v["key"]] = fetch_gnomad(v, build)
                queries += 1
                time.sleep(CFG["gnomad_delay"])
        newc = evaluate_gene(c["gene"], case["genes"][c["gene"]], ctx, case["afc"])
        newc["af_checked"] = True
        cands[c["gene"]] = newc
        if progress: progress(queries, budget)
    if queries >= budget: warnings.append(f"Достигнут лимит запросов к gnomAD ({budget}).")
    live = [c for c in cands.values() if not c["filtered"]]
    n_af_filtered = sum(1 for c in cands.values() if c["filtered"])
    rows = [build_candidate(c, ctx, build, population) for c in live]
    rows.sort(key=lambda r: -r["score"])
    for i, r in enumerate(rows, 1): r["rank"] = i; r["tier"] = tier_of(r["score"]); r["next_steps"] = next_steps(r, ctx)
    cols = ["rank", "gene", "variant", "variants", "zygosity", "consequence", "clinvar", "disease", "disease_id", "inheritance", "score", "pheno_log10", "moi_log10", "variant_log10", "af", "af_status", "tier", "n_matched", "n_unexplained", "evidence_db", "evidence_pred", "evidence_infer", "next_steps"]
    primary = pd.DataFrame(rows, columns=cols)
    verdict, top_score = "ok", (rows[0]["score"] if rows else None)
    if not rows or top_score < CFG["tier_possible"]: verdict = "no_convincing_candidate"
    elif len(rows) > 1 and rows[0]["score"] - rows[1]["score"] < CFG["competing_margin"]: verdict = "competing_hypotheses"
    carrier, secondary = carrier_and_secondary(case, ctx, acmg_genes)
    k = H.kb()
    meta = dict(verdict=verdict, build=build, population=population, sex=sex, known_moi=known_moi, n_present=len(ctx["present"]), n_excluded=len(ctx["excluded"]), hpoa_version=k.version, vep_server=VEP_SERVERS[build], gnomad=GNOMAD_DATASET[build], stats=dict(case["stats"], af_filtered_genes=n_af_filtered, gnomad_queries=queries), warnings=warnings, cfg=CFG)
    return dict(primary=primary, carrier=carrier, secondary=secondary, meta=meta)

def carrier_and_secondary(case, ctx, acmg_genes=None):
    k = H.kb()
    acmg = set(acmg_genes) if acmg_genes else ACMG_SF_V2
    carriers, secondary = [], []
    for gene, vs in case["genes"].items():
        dids = k.gene2dis.get(gene, [])
        mois = [k.diseases[d]["moi"] for d in dids]
        ar_only = bool(mois) and all(m and m <= {"AR"} for m in mois)
        for v in vs:
            if v["clinvar"] not in ("pathogenic", "likely_pathogenic"): continue
            if not (ar_only or gene in acmg): continue
            info = case["afc"].get(v["key"])
            if info is None:
                info = case["afc"][v["key"]] = fetch_gnomad(v, case["build"])
                time.sleep(CFG["gnomad_delay"])
            if info.get("status") == "found" and info["af_filter"] > (CFG["af_thr_recessive"] if ar_only else CFG["af_thr_dominant"]): continue
            dname = "; ".join(sorted({k.diseases[d]["name"] for d in dids}))[:300]
            row = dict(gene=gene, variant=_vdesc(v), zygosity=ZYG_RU.get(eff_zyg(v, ctx), v["zyg"]), clinvar=v["clinvar"], disease=dname, af_status=info.get("status"))
            if ar_only and eff_zyg(v, ctx) == "het" and v["chrom"] not in ("X", "Y", "MT"): carriers.append(row)
            if gene in acmg:
                note = "Гены ACMG SF: требуется консультация врача-генетика."
                if gene == "MUTYH": note += " MUTYH клинически значим при биаллельных вариантах."
                secondary.append(dict(row, recommendation=note))
    return (pd.DataFrame(carriers, columns=["gene", "variant", "zygosity", "clinvar", "disease", "af_status"]), pd.DataFrame(secondary, columns=["gene", "variant", "zygosity", "clinvar", "disease", "af_status", "recommendation"]))

def run_pipeline(tsv_path, json_path, population="nfe", build="GRCh38", sex="unknown", known_moi="unknown", progress=None, acmg_genes=None):
    df = pd.read_csv(tsv_path, sep="\t", dtype={"chrom": str})
    case = prepare_case(df, build=build, progress=progress)
    pheno = parse_phenopacket(json_path)
    return rank_case(case, pheno, sex=sex, known_moi=known_moi, population=population, acmg_genes=acmg_genes)