# haplogroups.py
"""
Гаплогруппы Y (Yleaf) и mtDNA (HaploGrep3). Результат не влияет на ранжирование.
Ключи командной строки зависят от версии программ - проверьте их через `<tool> --help`.
Вместо выдуманных значений всегда возвращается причина, почему гаплогруппа не определена.
"""
import glob
import os
import shutil
import subprocess
import tempfile


def _run(cmd, timeout=600):
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def determine_haplogroups(vcf_path, build="GRCh38"):
    out = {"y": None, "mt": None, "y_tool": False, "mt_tool": False, "source": "absent"}
    if not vcf_path or not os.path.exists(vcf_path):
        out["y"] = out["mt"] = "не определена (VCF отсутствует)"
        return out
    tmp = tempfile.mkdtemp(prefix="hg_")   # отдельная папка на запуск (не /tmp общий для всех пользователей)
    try:
        # Y: Yleaf
        if shutil.which("Yleaf"):
            out["y_tool"] = True
            rg = "hg38" if build == "GRCh38" else "hg19"
            try:
                r = _run(["Yleaf", "-vcf", vcf_path, "-o", os.path.join(tmp, "yleaf"), "-rg", rg])
                hits = glob.glob(os.path.join(tmp, "yleaf", "**", "*.hg"), recursive=True)
                if r.returncode == 0 and hits:
                    with open(hits[0], encoding="utf-8") as f:
                        lines = [l.rstrip("\n").split("\t") for l in f if l.strip()]
                    out["y"] = lines[1][1] if len(lines) > 1 and len(lines[1]) > 1 else "не определена (пустой результат)"
                else:
                    out["y"] = f"не определена (Yleaf код {r.returncode}: {r.stderr.strip()[-150:]})"
            except Exception as e:
                out["y"] = f"не определена (ошибка: {e})"
        else:
            out["y"] = "не определена (Yleaf не найден)"
        # mtDNA: HaploGrep3
        if shutil.which("haplogrep3"):
            out["mt_tool"] = True
            res = os.path.join(tmp, "mt.txt")
            try:
                r = _run(["haplogrep3", "classify", "--tree", "phylotree-rcrs@17.2", "--input", vcf_path, "--output", res])
                if r.returncode == 0 and os.path.exists(res):
                    with open(res, encoding="utf-8") as f:
                        rows = [l.rstrip("\n").replace('"', "").split("\t") for l in f if l.strip()]
                    out["mt"] = rows[1][1] if len(rows) > 1 and len(rows[1]) > 1 else "не определена (пустой результат)"
                else:
                    out["mt"] = f"не определена (HaploGrep3 код {r.returncode}: {r.stderr.strip()[-150:]})"
            except Exception as e:
                out["mt"] = f"не определена (ошибка: {e})"
        else:
            out["mt"] = "не определена (HaploGrep3 не найден)"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if out["y_tool"] or out["mt_tool"]:
        out["source"] = "tool"
    return out
