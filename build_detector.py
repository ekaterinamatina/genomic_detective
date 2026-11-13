# build_detector.py
"""Определение сборки генома по длинам контигов (надёжнее, чем поиск строк hg19/hg38)."""
import gzip
import re
import subprocess

# Длина хромосомы 1 в разных сборках
CHR1_LEN = {248956422: "GRCh38", 249250621: "GRCh37", 248387328: "T2T-CHM13"}


def _open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith((".gz", ".bgz")) \
        else open(path, encoding="utf-8", errors="replace")


def build_from_header_lines(lines):
    """lines - строки заголовка VCF или SAM. Возвращает 'GRCh37' | 'GRCh38' | 'T2T-CHM13' | None."""
    ref_hint = None
    for line in lines:
        m = re.search(r"##contig=<.*?ID=(?:chr)?1[,>].*?length=(\d+)", line)
        if m and int(m.group(1)) in CHR1_LEN:
            return CHR1_LEN[int(m.group(1))]
        m = re.match(r"@SQ\s+SN:(?:chr)?1\s+LN:(\d+)", line)
        if m and int(m.group(1)) in CHR1_LEN:
            return CHR1_LEN[int(m.group(1))]
        if line.startswith("##reference"):
            ref_hint = line.lower()
    if ref_hint:
        if any(x in ref_hint for x in ("grch38", "hg38")):
            return "GRCh38"
        if any(x in ref_hint for x in ("grch37", "hg19", "b37", "hs37")):
            return "GRCh37"
    return None


def detect_build_from_vcf(vcf_path):
    try:
        header = []
        with _open(vcf_path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                header.append(line.rstrip("\n"))
        return build_from_header_lines(header)
    except Exception:
        return None


def detect_build_from_bam(bam_path):
    try:
        r = subprocess.run(["samtools", "view", "-H", bam_path], capture_output=True, text=True, timeout=60)
        if r.returncode == 0:
            return build_from_header_lines(r.stdout.splitlines())
    except Exception:
        pass
    return None
