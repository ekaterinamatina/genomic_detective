import pandas as pd
import core_pipeline as core
import hpo_data as H

# 1. Загружаем варианты
df = pd.read_csv("case_05_variants.tsv", sep="\t")

# 2. Готовим кейс (запускаем VEP)
print("Запуск VEP (подождите ~20 секунд)...")
case = core.prepare_case(df, build="GRCh38")

# 3. Ищем LMNA
print("\n=== Что VEP выдал для chr1:156138613 ===")
found = False
for gene, variants in case["genes"].items():
    for v in variants:
        if v["pos"] == 156138613:
            print(f"Ген: {gene}")
            print(f"Мутация: {v['ref']}>{v['alt']}")
            print(f"Последствия: {v['consequence']}")
            print(f"ClinVar: {v['clinvar']}")
            print(f"Базовый LR варианта: {v['lr0']}")
            found = True

if not found:
    print("❌ Вариант chr1:156138613 вообще не попал в список генов (отсеян VEP или фильтрами)!")

# 4. Смотрим, как система его оценила (если нашел)
if found:
    pheno = core.parse_phenopacket("case_05_phenotype.json")
    res = core.rank_case(case, pheno, population="nfe")
    print("\n=== Итоговый рейтинг (топ-10) ===")
    print(res["primary"][["rank", "gene", "variant", "score", "pheno_log10", "variant_log10"]].head(10).to_string(index=False))