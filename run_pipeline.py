# run_pipeline.py
"""
Запуск из командной строки / из app.py:
    python run_pipeline.py wgs <R1> <R2> <workdir>
    python run_pipeline.py snp <file>    <workdir>
    python run_pipeline.py vcf <file>    <workdir>
В <workdir> должны лежать pheno.json (Phenopacket) и необязательный params.json:
    {"population": "nfe", "sex": "unknown", "known_moi": "unknown", "build": "GRCh38", "regions_bed": null}
Результат: <workdir>/result.json; ход выполнения: <workdir>/status.json (запись атомарная).
"""
import json
import os
import sys
import traceback
from pathlib import Path

import core_pipeline as core
from build_detector import detect_build_from_vcf
from haplogroups import determine_haplogroups
from pipeline_converter import convert_snp_to_vcf
from pipeline_wgs import run_wgs_pipeline
from vcf_utils import vcf_to_table


def update_status(workdir, step, progress, message=""):
    p = Path(workdir) / "status.json"
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps({"step": step, "progress": progress, "message": message}, ensure_ascii=False), encoding="utf-8")
    os.replace(tmp, p)


def df_to_records(df):
    return json.loads(df.to_json(orient="records", force_ascii=False))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    input_type, workdir = sys.argv[1], sys.argv[-1]
    files = sys.argv[2:-1]
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        params = {"population": "nfe", "sex": "unknown", "known_moi": "unknown", "build": None, "regions_bed": None}
        pp = wd / "params.json"
        if pp.exists():
            params.update(json.loads(pp.read_text(encoding="utf-8")))
        update_status(workdir, "start", 0, "Инициализация")

        build = params["build"]
        if input_type == "wgs":
            if len(files) != 2:
                raise ValueError("Для wgs нужны два файла: R1 и R2.")
            update_status(workdir, "wgs", 5, "fastp, выравнивание, MarkDuplicates, вызов вариантов (долго)")
            vcf_path = run_wgs_pipeline(files[0], files[1], "GRCh38", workdir)
            build = "GRCh38"
        elif input_type == "snp":
            update_status(workdir, "convert", 10, "Конвертация SNP-чипа в VCF")
            vcf_path, snp_build = convert_snp_to_vcf(files[0], str(wd / "converted.vcf"))
            build = build or snp_build
        elif input_type == "vcf":
            vcf_path = files[0]
        else:
            raise ValueError(f"Неизвестный тип входа: {input_type}")

        detected = detect_build_from_vcf(vcf_path)
        build = build or detected
        if build not in ("GRCh37", "GRCh38"):
            raise ValueError(f"Не удалось определить сборку генома (найдено: {detected}). "
                             "Укажите GRCh37 или GRCh38 явно - от неё зависит корректность аннотации.")
        if detected and detected != build:
            raise ValueError(f"Сборка в VCF ({detected}) не совпадает с заданной ({build}).")
        update_status(workdir, "build", 25, f"Сборка: {build}")

        df, vstats = vcf_to_table(vcf_path, regions_bed=params.get("regions_bed"), require_pass=True,
                                  min_dp=core.CFG["min_dp"], min_gq=core.CFG["min_gq"],
                                  max_variants=core.CFG["max_variants_total"])
        update_status(workdir, "vep", 30, f"Аннотация {len(df)} вариантов (Ensembl VEP)")
        case = core.prepare_case(df, build=build,
                                 progress=lambda i, n: update_status(workdir, "vep", 30 + int(40 * i / max(n, 1)),
                                                                     f"VEP: {i}/{n}"))
        update_status(workdir, "rank", 75, "Ранжирование и проверка частот gnomAD")
        pheno = core.parse_phenopacket(str(wd / "pheno.json"))
        res = core.rank_case(case, pheno, sex=params["sex"], known_moi=params["known_moi"],
                             population=params["population"])
        res["meta"]["vcf_stats"] = vstats
        update_status(workdir, "haplo", 92, "Гаплогруппы")
        hg = determine_haplogroups(vcf_path, build)

        out = {"status": "done", "build": build, "ranked": df_to_records(res["primary"]),
               "carrier": df_to_records(res["carrier"]), "secondary": df_to_records(res["secondary"]),
               "meta": json.loads(json.dumps(res["meta"], default=str)), "haplogroups": hg}
        (wd / "result.json").write_text(json.dumps(out, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        update_status(workdir, "done", 100, "Готово")
    except Exception as e:
        (wd / "error.log").write_text(traceback.format_exc(), encoding="utf-8")
        update_status(workdir, "error", 0, f"{type(e).__name__}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
