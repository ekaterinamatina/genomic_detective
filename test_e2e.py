"""E2E на заглушках VEP/gnomAD: проверка потока prepare_case -> rank_case, VCF-ридера и определения сборки."""
import os, tempfile, json, gzip
os.environ["GD_CACHE_DIR"] = tempfile.mkdtemp()
import pandas as pd
import hpo_data as H
import core_pipeline as C
from vcf_utils import vcf_to_table
from build_detector import detect_build_from_vcf

# --- заглушки знаний (как в test_logic)
DAG = {"HP:0000001": {"HP:0000001"}, "HP:0000010": {"HP:0000010", "HP:0000001"},
       "HP:0000011": {"HP:0000011", "HP:0000010", "HP:0000001"},
       "HP:0000012": {"HP:0000012", "HP:0000011", "HP:0000010", "HP:0000001"}}
H.ancestors = lambda h: frozenset(DAG.get(h, {h})); H.term_name = lambda h: h
H.canonical_id = lambda h: h if h in DAG else None
kb = H.KB(); kb.version = "test"
kb.diseases = {"OMIM:1": {"name": "Epilepsy AD", "ann": {"HP:0000012": 0.9}, "moi": {"AD"}},
               "OMIM:2": {"name": "CF-like AR", "ann": {"HP:0000010": 0.5}, "moi": {"AR"}},
               "OMIM:3": {"name": "Other", "ann": {"HP:0000001": 0.5}, "moi": {"AD"}}}
kb.gene2dis = {"GENE_E": ["OMIM:1"], "BRCA2": ["OMIM:3"], "CFTR": ["OMIM:2"]}
H._KB = kb; H.background_prob.cache_clear(); H.kb = lambda *a, **k: kb

def fake_vep(variants, build, warnings, progress=None):
    out = {}
    for c, p, r, a in variants:
        k = C.vkey(c, p, r, a)
        gene = {10: "GENE_E", 20: "BRCA2", 30: "CFTR", 40: "GENE_E"}.get(p)
        if p == 99: out[k] = []; continue
        if gene is None: continue
        clin = {20: "pathogenic", 30: "pathogenic", 40: "likely_benign"}.get(p, "uncertain")
        g = dict(gene=gene, consequence=["missense_variant"], all_terms=["missense_variant"], impact="MODERATE",
                 transcript="T", hgvsc="", hgvsp="", clinvar=clin, clinvar_ambiguous=False,
                 revel=0.95 if p == 10 else None, sift="", polyphen="")
        out[k] = [g]
    return out, 0
C.annotate_with_vep = fake_vep
calls = []
def fake_gn(v, build):
    calls.append(v["key"])
    return dict(status="absent") if v["pos"] == 10 else dict(status="found", ac=5, an=100000, af_global=5e-5,
            af_filter=5e-5, af_filter_group="nfe", group_af={"nfe": 5e-5})
C.fetch_gnomad = fake_gn
C.CFG["gnomad_delay"] = 0

df = pd.DataFrame({"chrom": ["chr1", "1", "1", "1", "1", "1"], "pos": [10, 20, 30, 40, 99, 50],
                   "ref": ["A"] * 6, "alt": ["T"] * 6, "genotype": ["0/1", "0/1", "0/1", "0/1", "0/1", "0/0"]})
case = C.prepare_case(df, "GRCh38")
print("stats:", case["stats"])
assert case["stats"]["not_carried"] == 1                 # 0/0 не анализируется
assert case["stats"]["clinvar_benign_removed"] == 1      # likely_benign отброшен
pheno = dict(present=["HP:0000012"], excluded=[], unknown=[], conflicts=[], sex="unknown")
res = C.rank_case(case, pheno, population="nfe")
p = res["primary"]
print(p[["rank", "gene", "score", "tier", "af_status", "inheritance", "zygosity"]].to_string())
assert p.iloc[0]["gene"] == "GENE_E"
brca = p[p.gene == "BRCA2"].iloc[0]
assert brca["score"] <= C.CFG["no_match_cap"] and any("не объясняет" in x for x in brca["evidence_infer"])
assert p.iloc[0]["af_status"] == "absent"
# CFTR: P-вариант, het при AR -> носитель, но не топ-кандидат
assert "CFTR" in res["carrier"]["gene"].tolist()
assert "BRCA2" in res["secondary"]["gene"].tolist()
print("verdict:", res["meta"]["verdict"], "| gnomAD queries:", len(calls))

# изменение фенотипа/наследования пересчитывается без новой аннотации (VEP-заглушка не вызывается повторно)
C.annotate_with_vep = lambda *a, **k: (_ for _ in ()).throw(AssertionError("повторная аннотация"))
res2 = C.rank_case(case, dict(pheno, present=["HP:0000010"]), known_moi="AR", population="nfe")
print(res2["primary"][["rank", "gene", "score"]].to_string())
assert res2["primary"].iloc[0]["gene"] in ("CFTR", "GENE_E", "BRCA2")

# --- VCF-ридер
vcf = os.path.join(tempfile.mkdtemp(), "t.vcf")
open(vcf, "w").write("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:60
chr1\t200\t.\tC\tT,G\t50\t.\t.\tGT:DP:GQ\t1/2:30:60
chr1\t300\t.\tC\tT\t50\tLowQual\t.\tGT:DP:GQ\t0/1:30:60
chr1\t400\t.\tC\t<DEL>\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:60
chr1\t500\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/0:30:60
chr1\t600\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:3:60
""")
t, st = vcf_to_table(vcf, min_dp=8, min_gq=20)
print(t.to_string(), st)
assert len(t) == 3 and all(C.parse_gt(g) == "het" for g in t["genotype"]) and detect_build_from_vcf(vcf) == "GRCh38"
print("E2E PASSED")
