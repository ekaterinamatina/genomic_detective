# pipeline_converter.py
import subprocess
from pathlib import Path


def convert_snp_to_vcf(input_file, out_vcf):
    """
    23andMe / AncestryDNA -> VCF (библиотека snps). Возвращает путь к реально созданному файлу.
    Внимание: чип покрывает ~0.6-0.7 млн в основном частых SNP; для диагностики редких болезней
    непригоден (большинство патогенных вариантов на чипе отсутствует).
    """
    try:
        from snps import SNPs
    except ImportError:
        raise RuntimeError("Для конвертации SNP-чипов нужна библиотека 'snps': pip install snps")
    out = Path(out_vcf).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    s = SNPs(str(input_file), output_dir=str(out.parent))   # snps пишет в output_dir, а не по произвольному пути
    if not s.valid:
        raise ValueError("Файл SNP не распознан библиотекой snps.")
    result = s.to_vcf(out.name)
    produced = Path(result) if result else out
    if not produced.exists():
        produced = out
    if not produced.exists():
        raise RuntimeError("snps не создал VCF-файл.")
    build = {37: "GRCh37", 38: "GRCh38"}.get(getattr(s, "build", None))
    return str(produced), build


def convert_plink_to_vcf(bed_prefix, out_prefix, ref_fasta):
    """PLINK bed/bim/fam -> VCF (plink2). REF-аллели берутся из референса (--ref-from-fa)."""
    subprocess.run(["plink2", "--bfile", bed_prefix, "--fa", ref_fasta, "--ref-from-fa", "force",
                    "--export", "vcf", "bgz", "--output-chr", "MT", "--out", out_prefix], check=True)
    return out_prefix + ".vcf.gz"
