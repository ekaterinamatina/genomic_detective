# vcf_utils.py
"""Чтение VCF в таблицу вариантов (без внешних зависимостей)."""
import bisect
import gzip

import pandas as pd

SKIP_ALT_CHARS = ("<", "[", "]", "*")


def _open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith((".gz", ".bgz")) \
        else open(path, encoding="utf-8", errors="replace")


def norm_chrom(c):
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return "MT" if c in ("M", "MT", "m", "mt") else c


def load_bed(path):
    """BED -> {chrom: (starts, ends)} для быстрой проверки попадания позиции в регионы."""
    tmp = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.split("\t")
            tmp.setdefault(norm_chrom(p[0]), []).append((int(p[1]), int(p[2])))
    out = {}
    for c, iv in tmp.items():
        iv.sort()
        out[c] = ([a for a, _ in iv], [b for _, b in iv])
    return out


def in_bed(bed, chrom, pos1):
    """pos1 - 1-based позиция."""
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    i = bisect.bisect_right(starts, pos1 - 1) - 1
    # интервалы могут перекрываться - проверяем несколько предыдущих
    for j in range(i, max(i - 50, -1), -1):
        if ends[j] >= pos1:
            return True
    return False


def vcf_to_table(vcf_path, sample=None, regions_bed=None, require_pass=True,
                 min_dp=None, min_gq=None, max_variants=None):
    """
    Возвращает (DataFrame[chrom,pos,ref,alt,genotype,dp,gq,filter], stats).
    Мультиаллельные сайты разбиваются; GT перекодируется относительно каждого альтернативного аллеля.
    Берётся первый образец (или sample).
    """
    bed = load_bed(regions_bed) if regions_bed else None
    rows = []
    st = {"records": 0, "non_pass": 0, "symbolic": 0, "ref_only": 0, "low_dp": 0, "low_gq": 0,
          "outside_regions": 0, "kept": 0}
    sample_idx = None
    with _open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 10:
                    raise ValueError("В VCF нет столбцов с генотипами (нужен хотя бы один образец).")
                names = cols[9:]
                sample_idx = names.index(sample) if sample else 0
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            st["records"] += 1
            chrom, pos, ref, alts, flt = norm_chrom(p[0]), int(p[1]), p[3].upper(), p[4].upper(), p[6]
            if require_pass and flt not in ("PASS", "."):
                st["non_pass"] += 1
                continue
            fmt = p[8].split(":")
            vals = p[9 + sample_idx].split(":")
            fd = dict(zip(fmt, vals))
            gt = fd.get("GT", "./.")

            def num(k):
                try:
                    return float(fd[k])
                except Exception:
                    return None
            dp, gq = num("DP"), num("GQ")
            if min_dp is not None and dp is not None and dp < min_dp:
                st["low_dp"] += 1
                continue
            if min_gq is not None and gq is not None and gq < min_gq:
                st["low_gq"] += 1
                continue
            if bed is not None and not in_bed(bed, chrom, pos):
                st["outside_regions"] += 1
                continue
            gt_alleles = gt.replace("|", "/").split("/")
            for k, alt in enumerate(alts.split(","), start=1):
                if alt == "." or any(c in alt for c in SKIP_ALT_CHARS):
                    st["symbolic"] += 1
                    continue
                recoded = "/".join("1" if a == str(k) else ("0" if a not in (".",) else ".") for a in gt_alleles)
                if "1" not in recoded:
                    st["ref_only"] += 1
                    continue
                rows.append((chrom, pos, ref, alt, recoded, dp, gq, flt))
                st["kept"] += 1
                if max_variants and len(rows) > max_variants:
                    raise ValueError(
                        f"После фильтров осталось более {max_variants} вариантов. Аннотация через публичный "
                        "REST-сервис VEP для такого объёма непрактична: ограничьте анализ экзонами "
                        "(параметр regions_bed, например BED экзонов MANE) или используйте локальный VEP."
                    )
    df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "genotype", "dp", "gq", "filter"])
    return df, st
