"""Тест логики скоринга на синтетических данных (без сети и pyhpo)."""
import math, os, tempfile
os.environ["GD_CACHE_DIR"] = tempfile.mkdtemp()
import hpo_data as H
import core_pipeline as C

# Мини-онтология: ROOT > NERV > SEIZ > FOCAL ; NERV > ATAXIA ; ROOT > SKEL > ARACH
DAG = {"HP:0000001": {"HP:0000001"},
       "HP:0000010": {"HP:0000010", "HP:0000001"},            # NERV
       "HP:0000011": {"HP:0000011", "HP:0000010", "HP:0000001"},  # SEIZ
       "HP:0000012": {"HP:0000012", "HP:0000011", "HP:0000010", "HP:0000001"},  # FOCAL
       "HP:0000013": {"HP:0000013", "HP:0000010", "HP:0000001"},  # ATAXIA
       "HP:0000020": {"HP:0000020", "HP:0000001"},            # SKEL
       "HP:0000021": {"HP:0000021", "HP:0000020", "HP:0000001"}}  # ARACH
H.ancestors = lambda h: frozenset(DAG.get(h, {h}))
H.term_name = lambda h: h
kb = H.KB()
kb.diseases = {
 "OMIM:1": {"name": "Epilepsy AD", "ann": {"HP:0000012": 0.9, "HP:0000013": 0.5}, "moi": {"AD"}},
 "OMIM:2": {"name": "Ataxia AR",  "ann": {"HP:0000013": 0.9}, "moi": {"AR"}},
 "OMIM:3": {"name": "Marfanoid",   "ann": {"HP:0000021": 0.9}, "moi": {"AD"}},
 "OMIM:4": {"name": "Other",       "ann": {"HP:0000021": 0.5}, "moi": {"AD"}},
 "OMIM:5": {"name": "Other2",      "ann": {"HP:0000021": 0.5}, "moi": {"AR"}},
}
kb.gene2dis = {"GENE_E": ["OMIM:1"], "GENE_A": ["OMIM:2"], "GENE_M": ["OMIM:3"], "GENE_X": ["OMIM:4"], "GENE_Y": ["OMIM:5"]}
H._KB = kb
H.background_prob.cache_clear()
H.kb = lambda *a, **k: kb

def var(gene, pos, zyg, clin="uncertain", revel=None, terms=("missense_variant",), chrom="1"):
    v = dict(gene=gene, consequence=list(terms), all_terms=list(terms), impact="MODERATE", transcript="T", hgvsc="", hgvsp="",
             clinvar=clin, clinvar_ambiguous=False, revel=revel, sift="", polyphen="", key=f"{chrom}-{pos}-A-T",
             chrom=chrom, pos=pos, ref="A", alt="T", zyg=zyg)
    v["lr0"], v["lr0_parts"], v["lr0_source"] = C.variant_base_lr(v)
    return v

# 1) classify_clinvar: подстрока 'pathogenic' не должна срабатывать на conflicting
assert C.classify_clinvar(["conflicting_interpretations_of_pathogenicity"]) == "conflicting"
assert C.classify_clinvar(["pathogenic/likely_pathogenic"]) == "pathogenic"
assert C.classify_clinvar(["likely_pathogenic"]) == "likely_pathogenic"
assert C.classify_clinvar(["uncertain_significance"]) == "uncertain"
# 2) genotype
assert C.parse_gt("0/0") == "ref" and C.parse_gt("./.") == "missing" and C.parse_gt("1|0") == "het"
assert C.parse_gt("1/1") == "hom" and C.parse_gt("1") == "hemi"
# 3) LR ClinVar
print("ClinVar LR:", {k: round(v, 1) for k, v in C.CLINVAR_LR.items()})

pheno = dict(present=["HP:0000012", "HP:0000013"], excluded=[], unknown=[], conflicts=[], sex="unknown")
ctx = dict(present=pheno["present"], excluded=[], sex="unknown", known_moi="unknown", afc={})
afc = {}

# 4) AD ген с подходящим фенотипом vs. ClinVar P в «чужом» гене
eE = C.evaluate_gene("GENE_E", [var("GENE_E", 10, "het", revel=0.95)], ctx, afc)
eM = C.evaluate_gene("GENE_M", [var("GENE_M", 20, "het", clin="pathogenic")], ctx, afc)
print("GENE_E", round(eE["best"]["total"], 2), "GENE_M(P, нет совпадения)", round(eM["best"]["total"], 2))
assert eE["best"]["total"] > eM["best"]["total"], "фенотипически подходящий ген должен обгонять несовпадающий P-вариант"

# 5) AR: один het << hom
eA1 = C.evaluate_gene("GENE_A", [var("GENE_A", 30, "het", revel=0.95)], ctx, afc)
eA2 = C.evaluate_gene("GENE_A", [var("GENE_A", 30, "hom", revel=0.95)], ctx, afc)
eA3 = C.evaluate_gene("GENE_A", [var("GENE_A", 30, "het", revel=0.95), var("GENE_A", 31, "het", revel=0.95)], ctx, afc)
t = lambda e: round(e["best"]["total"], 2)
print("AR het/hom/comphet:", t(eA1), t(eA2), t(eA3))
assert t(eA2) > t(eA3) > t(eA1)

# 6) Фильтр по частоте: AF 5% исключает AD-вариант, но не рецессивный
afc_hi = {"1-10-A-T": dict(status="found", af_filter=0.05, af_filter_group="nfe", af_global=0.04, group_af={"nfe": 0.05}, ac=1, an=1)}
eF = C.evaluate_gene("GENE_E", [var("GENE_E", 10, "het", revel=0.95)], ctx, afc_hi)
assert eF["filtered"], "AD-вариант с AF 5% должен быть отфильтрован"
afc_hi2 = {"1-30-A-T": dict(status="found", af_filter=0.005, af_filter_group="nfe", af_global=0.004, group_af={"nfe": 0.005}, ac=1, an=1)}
eR = C.evaluate_gene("GENE_A", [var("GENE_A", 30, "hom", revel=0.95)], ctx, afc_hi2)
assert not eR["filtered"], "AF 0.5% не выше порога 2% для AR"

# 7) Исключённые признаки реально снижают балл болезни, у которой признак типичен, и не влияют на посторонние
ctx_ex = dict(ctx, excluded=["HP:0000013"])
base = C.evaluate_gene("GENE_E", [var("GENE_E", 10, "het", revel=0.95)], ctx, afc)["best"]["pl"]
pen = C.evaluate_gene("GENE_E", [var("GENE_E", 10, "het", revel=0.95)], ctx_ex, afc)["best"]["pl"]
unrel = C.evaluate_gene("GENE_M", [var("GENE_M", 20, "het")], ctx, afc)["best"]["pl"]
unrel2 = C.evaluate_gene("GENE_M", [var("GENE_M", 20, "het")], ctx_ex, afc)["best"]["pl"]
print("excluded penalty:", round(base, 2), "->", round(pen, 2), "| unrelated:", round(unrel, 2), round(unrel2, 2))
assert pen < base and abs(unrel - unrel2) < 1e-9

# 8) Известное наследование: AD заданное пользователем снижает AR-болезнь
ctx_ad = dict(ctx, known_moi="AD")
a_ad = C.evaluate_gene("GENE_A", [var("GENE_A", 30, "hom", revel=0.95)], ctx_ad, afc)["best"]["total"]
assert a_ad < eA2["best"]["total"] - 1.0

# 9) X-сцепленность: мужской пол, пол неизвестен
H_xl = H.KB.__new__(H.KB)
kb.diseases["OMIM:6"] = {"name": "XLR dis", "ann": {"HP:0000013": 0.9}, "moi": {"XLR"}}
kb.gene2dis["GENE_XR"] = ["OMIM:6"]
vx = var("GENE_XR", 40, "hom", chrom="X", revel=0.95)
male = C.evaluate_gene("GENE_XR", [vx], dict(ctx, sex="male"), afc)["best"]["total"]
vxf = var("GENE_XR", 41, "het", chrom="X", revel=0.95)
fem = C.evaluate_gene("GENE_XR", [vxf], dict(ctx, sex="female"), afc)["best"]["total"]
print("XLR male hemi", round(male, 2), "female het", round(fem, 2))
assert male > fem

# 10) tier / параметры
assert C.tier_of(7) == "сильное свидетельство" and C.tier_of(0.5) == "недостаточно данных"
print("ALL TESTS PASSED")
