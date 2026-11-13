# pipeline_wgs.py
"""
FASTQ -> VCF: fastp -> bwa-mem2 -> MarkDuplicates (GATK) -> bcftools call -> фильтрация.

Примечания:
  * bcftools call подходит для прототипа; для диагностики индели надёжнее вызывать GATK HaplotypeCaller
    или DeepVariant (Poplin et al., Nat Biotechnol 2018;36:983).
  * Полногеномные данные 30x - сотни гигабайт: загрузка через веб-форму и синхронный запуск в Streamlit
    для них не годятся (нужны пути к файлам на сервере и очередь задач).
"""
import shlex
import shutil
import subprocess
from pathlib import Path

REQUIRED_TOOLS = ["fastp", "bwa-mem2", "samtools", "gatk", "bcftools"]


def check_tools():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        raise RuntimeError("Не найдены программы: " + ", ".join(missing))


def run_wgs_pipeline(r1, r2, build, workdir, threads=8, sample="SAMPLE", ref_root="refs"):
    check_tools()
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    ref = Path(ref_root) / build / "genome.fa"
    for suffix in ("", ".fai"):
        if not Path(str(ref) + suffix).exists():
            raise FileNotFoundError(f"Нет референса/индекса: {ref}{suffix}")
    if not Path(str(ref) + ".bwt.2bit.64").exists():
        raise FileNotFoundError(f"Нет индекса bwa-mem2 для {ref} (bwa-mem2 index {ref})")

    q = shlex.quote
    r1_trim, r2_trim = str(wd / "R1.trim.fastq.gz"), str(wd / "R2.trim.fastq.gz")
    subprocess.run(["fastp", "-i", r1, "-I", r2, "-o", r1_trim, "-O", r2_trim,
                    "--qualified_quality_phred", "20", "--length_required", "50",
                    "--thread", str(threads), "-j", str(wd / "fastp.json"), "-h", str(wd / "fastp.html")], check=True)

    bam = str(wd / "aligned.bam")
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}"   # read group нужен GATK
    subprocess.run(
        f"bwa-mem2 mem -t {threads} -R {q(rg)} {q(str(ref))} {q(r1_trim)} {q(r2_trim)} | "
        f"samtools sort -@ {threads} -o {q(bam)} -", shell=True, check=True, executable="/bin/bash")
    subprocess.run(["samtools", "index", bam], check=True)

    bam_md = str(wd / "aligned.md.bam")
    subprocess.run(["gatk", "MarkDuplicates", "-I", bam, "-O", bam_md, "-M", str(wd / "md.metrics.txt")], check=True)
    subprocess.run(["samtools", "index", bam_md], check=True)

    vcf_raw = str(wd / "raw.vcf.gz")
    subprocess.run(
        f"bcftools mpileup -f {q(str(ref))} -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -Ou {q(bam_md)} | "
        f"bcftools call -mv -Ou | bcftools norm -f {q(str(ref))} -m-any -Oz -o {q(vcf_raw)}",
        shell=True, check=True, executable="/bin/bash")
    subprocess.run(["bcftools", "index", "-t", vcf_raw], check=True)

    # bcftools call пишет FILTER='.', поэтому '-f PASS' удалил бы ВСЕ варианты; нужно '.,PASS'
    vcf_final = str(wd / "final.vcf.gz")
    subprocess.run(["bcftools", "view", "-f", ".,PASS", "-i", "QUAL>=30 && FORMAT/DP>=10",
                    "-Oz", "-o", vcf_final, vcf_raw], check=True)
    subprocess.run(["bcftools", "index", "-t", vcf_final], check=True)
    return vcf_final
