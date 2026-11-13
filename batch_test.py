# batch_test.py
"""
Пакетная проверка на размеченных случаях с метриками.
Ожидаемые файлы в папках development/ и hidden/:
    <case>_variants.tsv, <case>_phenotype.json  (как раньше)
    truth.tsv : два столбца без заголовка - case_id <TAB> ген-причина   (формат предложен здесь; при другом
                формате правится функция load_truth)
Пороги и веса (core_pipeline.CFG) подбираются ТОЛЬКО на development; hidden - для финальной проверки.
"""
import os
import sys

import core_pipeline as core


def load_truth(folder):
    path = os.path.join(folder, "truth.tsv")
    truth = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[0] and not p[0].startswith("#"):
                    truth[p[0]] = {g.strip() for g in p[1].split(",") if g.strip()}
    return truth


def evaluate_folder(folder, population="nfe", build="GRCh38"):
    truth = load_truth(folder)
    ranks = []
    for fn in sorted(f for f in os.listdir(folder) if f.endswith("_variants.tsv")):
        case = fn[:-len("_variants.tsv")]
        pheno = os.path.join(folder, f"{case}_phenotype.json")
        if not os.path.exists(pheno):
            print(f"[{case}] нет файла фенотипа"); continue
        try:
            res = core.run_pipeline(os.path.join(folder, fn), pheno, population=population, build=build)
        except Exception as e:
            print(f"[{case}] ОШИБКА: {e}"); continue
        df = res["primary"]
        print(f"--- {case}: вердикт={res['meta']['verdict']}")
        for _, r in df.head(3).iterrows():
            af = f"AF={r['af']:.2g}" if r["af"] == r["af"] and r["af"] is not None else f"AF:{r['af_status']}"
            print(f"  {r['rank']}. {r['gene']:<8} score={r['score']:<6} {r['variant']} {af} [{r['tier']}]")
        if case in truth:
            pos = [int(r["rank"]) for _, r in df.iterrows() if r["gene"] in truth[case]]
            ranks.append(min(pos) if pos else None)
            print(f"  причинный ген {sorted(truth[case])}: ранг {ranks[-1] or 'не найден'}")
    return ranks


def summarize(name, ranks):
    n = len(ranks)
    if not n:
        print(f"{name}: нет размеченных случаев (truth.tsv)"); return
    hit = lambda k: sum(1 for r in ranks if r and r <= k)
    mrr = sum(1 / r for r in ranks if r) / n
    print(f"{name}: N={n}  top1={hit(1)}/{n}  top3={hit(3)}/{n}  top5={hit(5)}/{n}  MRR={mrr:.3f}")


if __name__ == "__main__":
    build = sys.argv[1] if len(sys.argv) > 1 else "GRCh38"
    for folder in ("development", "hidden"):
        if os.path.isdir(folder):
            summarize(folder, evaluate_folder(folder, build=build))
        else:
            print(f"[WARN] папка {folder} не найдена")
