# Tiers

⭢ Link do Google Colab para a visualização das tabelas.

https://colab.research.google.com/drive/1xLRH-LPHSKECMB5pYR3pKcw9iIGSR6n1?usp=sharing

## 1. Clonagem de GitHub

```
! git clone https://github.com/renatopuga/lmabrasil-hg38.git
```

### Output:
```
Cloning into 'lmabrasil-hg38'...
remote: Enumerating objects: 226, done.
remote: Counting objects: 100% (168/168), done.
remote: Compressing objects: 100% (108/108), done.
remote: Total 226 (delta 90), reused 114 (delta 56), pack-reused 58 (from 1)
Receiving objects: 100% (226/226), 8.63 MiB | 22.96 MiB/s, done.
Resolving deltas: 100% (106/106), done.
```

## 2. Instalar bcftools (+split-vep é uma ferramenta dentro de bcftools)

```
%%bash
sudo apt install bcftools
```

```
!bcftools +split-vep
```

```
!ls lmabrasil-hg38/vep_output/*.vep.vcf
```

# WP017

### Etapa 1
```
%%bash
VEP_VCF="lmabrasil-hg38/vep_output/liftOver_WP017_hg19ToHg38.vep.vcf"

# criar o cabeçalho
bcftools +split-vep -l $VEP_VCF | \
cut -f2  | \
tr '\n\r' '\t' | \
awk '{print("CHROM\tPOS\tREF\tALT\t"$0"FILTER\tTumorID\tGT\tDP\tAD\tAF\tNormalID\tNGT\tNDP\tNAD\tNAF")}' \
> liftOver_WP017_hg19ToHg38.vep.tsv

# adicionar as variantes
bcftools +split-vep \
-f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%FILTER\t[%SAMPLE\t%GT\t%DP\t%AD\t%AF\t]\n' \
-i 'FMT/DP>=20 && FMT/AF>=0.1' -d -A tab $VEP_VCF \
-p x  >> liftOver_WP017_hg19ToHg38.vep.tsv
```

### Etapa 2
```
import pandas as pd
import numpy as np
import re
import requests

# =========================
# CONFIG
# =========================
INPUT_TSV  = "/content/liftOver_WP017_hg19ToHg38.vep.tsv"
OUTPUT_TSV = "/content/WP017-tier.tsv"

DRIVER_RAW_URL = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES
# =========================
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def pick_first_existing_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def parse_percent_vaf(x):
    # retorna VAF em % (0-100). Se vier 0-1, converte para 0-100.
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", ".", "NA", "NaN", "nan", "None"}:
        return np.nan
    s = s.replace(",", ".").replace("%", "")
    try:
        v = float(s)
    except:
        return np.nan
    return v * 100.0 if v <= 1.0 else v

def filter_is_pass(x):
    return str(x).strip().upper() == "PASS"

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    # seu arquivo raw está em UMA linha com genes separados por espaço :contentReference[oaicite:1]{index=1}
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)

# =========================
# LOAD
# =========================
driver_genes = load_driver_genes(DRIVER_RAW_URL)
print("Driver genes carregados:", len(driver_genes))

df = smart_read_tsv(INPUT_TSV)

# detectar colunas (ajuste se necessário)
gene_col     = pick_first_existing_col(df, ["SYMBOL", "Gene", "GENE", "HGNC"])
filter_col   = pick_first_existing_col(df, ["FILTER", "Filter"])
dp_col       = pick_first_existing_col(df, ["DP", "DP_tumor", "TUMOR_DP", "DEPTH"])
vaf_col      = pick_first_existing_col(df, ["VAF_tumor", "VAF", "AF", "TUMOR_AF"])
existing_col = pick_first_existing_col(df, ["Existing_variation", "ExistingVariation", "existing_variation"])

missing = [("SYMBOL", gene_col), ("FILTER", filter_col), ("DP", dp_col), ("VAF/AF", vaf_col), ("Existing_variation", existing_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas necessárias: {missing}\nColunas do arquivo: {list(df.columns)}")

# =========================
# TIERS (hierárquico)
# =========================
df["_GENE_UP"] = df[gene_col].astype(str).str.strip().str.upper()
df["_PASS"]    = df[filter_col].apply(filter_is_pass)
df["_DP"]      = pd.to_numeric(df[dp_col].astype(str).str.replace(",", ".", regex=False), errors="coerce").fillna(0).astype(int)
df["_VAF_PCT"] = df[vaf_col].apply(parse_percent_vaf)
df["_EXIST"]   = df[existing_col].astype(str)

base  = (df["_PASS"]) & (df["_DP"] >= 20) & (df["_VAF_PCT"] >= 10)

tier1 = base & (df["_GENE_UP"].isin(driver_genes))
tier2 = base & (~tier1) & (df["_EXIST"].str.contains(r"\bCOSV\b", flags=re.IGNORECASE, regex=True))
tier3 = ~(tier1 | tier2)

df["Tier"] = np.select([tier1, tier2, tier3], ["Tier 1", "Tier 2", "Tier 3"], default="Tier 3")

# (opcional) auditoria rápida
df["Tier_rule_PASS"]    = df["_PASS"]
df["Tier_rule_DP"]      = df["_DP"]
df["Tier_rule_VAF_pct"] = df["_VAF_PCT"]

df.drop(columns=[c for c in df.columns if c.startswith("_")], inplace=True)
df.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("OK ->", OUTPUT_TSV)
print(df["Tier"].value_counts(dropna=False))
```

### Output:
```
Driver genes carregados: 64
OK -> /content/WP017-tier.tsv
Tier
Tier 3    2226
Tier 1       8
Name: count, dtype: int64
```

### Etapa 3
```
# =========================
# IMPORTS
# =========================
import pandas as pd
import numpy as np
import requests

# =========================
# ARQUIVOS
# =========================
input_tsv  = "/content/liftOver_WP017_hg19ToHg38.vep.tsv"
output_tsv = "/content/WP017_tier_table.tsv"

driver_url = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES (Re-definir ou garantir que smart_read_tsv esteja disponível)
# =========================
# A função smart_read_tsv já está definida na célula anterior (uzRA_DALNTe5)
# e pode ser reutilizada aqui. Se esta célula for executada isoladamente,
# a função precisaria ser copiada ou importada.
# Por simplicidade e dado o contexto do notebook, vamos assumir que está disponível.

# Função smart_read_tsv (copiada para garantir disponibilidade se a célula anterior não for executada)
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)


# =========================
# LER VEP
# =========================
df = smart_read_tsv(input_tsv)

# =========================
# GENES DRIVER
# =========================
driver_genes = load_driver_genes(driver_url)

df
```

#### ↳ O Output esperado é uma longa tabela e você encontra ela no Colab. 

### Etapa 4 
```
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import zipfile
import os

# Define the path to the zip file
ZIP_ARCHIVE_PATH = "/content/cgi_somaticmfwp058_results.zip"
# Define the directory where content will be extracted
EXTRACT_TO_DIR = "/content/extracted_cgi_data"

# Ensure the extraction directory exists
os.makedirs(EXTRACT_TO_DIR, exist_ok=True)

extracted_alterations_tsv_path = None
try:
    with zipfile.ZipFile(ZIP_ARCHIVE_PATH, 'r') as zip_ref:
        # Search for 'alterations.tsv' within the zip file.
        # This handles cases where it might be directly at root or nested in a folder.
        for member in zip_ref.namelist():
            if member.endswith('alterations.tsv'):
                zip_ref.extract(member, EXTRACT_TO_DIR)
                extracted_alterations_tsv_path = os.path.join(EXTRACT_TO_DIR, member)
                break

        if not extracted_alterations_tsv_path:
            raise FileNotFoundError(f"Could not find 'alterations.tsv' within {ZIP_ARCHIVE_PATH}")

    # Update ALTERATIONS_TSV to point to the newly extracted file
    ALTERATIONS_TSV = extracted_alterations_tsv_path
    print(f"'{os.path.basename(ALTERATIONS_TSV)}' extracted and path updated to: {ALTERATIONS_TSV}")

except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the zip file exists and contains 'alterations.tsv'.")
    raise
except zipfile.BadZipFile:
    print(f"Error: '{ZIP_ARCHIVE_PATH}' is not a valid zip file.")
    raise
except Exception as e:
    print(f"An unexpected error occurred during zip extraction: {e}")
    raise

TIERS_TSV       = "/content/WP058-tier.tsv"    # seu TSV com coluna Tier (Tier 1/2/3)

def pick_first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def norm_chr(x):
    s = str(x).strip()
    if not s:
        return s
    return s if s.lower().startswith("chr") else "chr" + s

def make_var_id(df, chr_col, pos_col, ref_col, alt_col):
    return (
        df[chr_col].astype(str).map(norm_chr) + ":" +
        df[pos_col].astype(str) + ":" +
        df[ref_col].astype(str) + ":" +
        df[alt_col].astype(str)
    )

# ----------------------------
# 1) DRIVER set (a partir do alterations.tsv)
# ----------------------------
alt = pd.read_csv(ALTERATIONS_TSV, sep="\t", dtype=str)

chr_col = pick_first_existing_col(alt, ["CHR", "CHROM", "#CHROM", "CHROMOSOME"])
pos_col = pick_first_existing_col(alt, ["POS", "POSITION"])
ref_col = pick_first_existing_col(alt, ["REF"])
alt_col = pick_first_existing_col(alt, ["ALT"])
pred_col = pick_first_existing_col(alt, ["CGI-Oncogenic Prediction", "CGI-Oncogenic Summary"])

missing = [("CHR", chr_col), ("POS", pos_col), ("REF", ref_col), ("ALT", alt_col), ("DriverColumn", pred_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no alterations.tsv: {missing}\nColunas: {list(alt.columns)}")

alt["VAR_ID"] = make_var_id(alt, chr_col, pos_col, ref_col, alt_col)

# considera driver se a coluna CGI tiver a palavra "driver"
is_driver = alt[pred_col].fillna("").astype(str).str.contains(r"\bdriver\b", flags=re.IGNORECASE, regex=True)
set_driver = set(alt.loc[is_driver, "VAR_ID"])

# ----------------------------
# 2) T1 set (a partir do tiers TSV)
# ----------------------------
tiers = pd.read_csv(TIERS_TSV, sep="\t", dtype=str)

t_chr = pick_first_existing_col(tiers, ["CHROM", "#CHROM", "CHR", "Chrom"])
t_pos = pick_first_existing_col(tiers, ["POS", "Position", "pos"])
t_ref = pick_first_existing_col(tiers, ["REF", "Ref"])
t_alt = pick_first_existing_col(tiers, ["ALT", "Alt"])

missing = [("CHROM", t_chr), ("POS", t_pos), ("REF", t_ref), ("ALT", t_alt)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no tiers TSV: {missing}\nColunas: {list(tiers.columns)}")

tiers["VAR_ID"] = make_var_id(tiers, t_chr, t_pos, t_ref, t_alt)
set_t1 = set(tiers.loc[tiers["Tier"].astype(str) == "Tier 1", "VAR_ID"])

# ----------------------------
# 3) Venn
# ----------------------------
plt.figure(figsize=(6, 6))
venn2([set_driver, set_t1], set_labels=("Driver (alterations.tsv)", "Tier 1 (tiers.tsv)"))
plt.title("Venn: chr:pos:ref:alt — Driver vs Tier 1")
plt.show()

print("N Driver:", len(set_driver))
print("N Tier 1:", len(set_t1))
print("Interseção:", len(set_driver & set_t1))
print("Driver somente:", len(set_driver - set_t1))
print("T1 somente:", len(set_t1 - set_driver))

# (Opcional) salvar listas
pd.Series(sorted(set_driver & set_t1)).to_csv("/content/venn_intersection.txt", index=False, header=False)
pd.Series(sorted(set_driver - set_t1)).to_csv("/content/venn_driver_only.txt", index=False, header=False)
pd.Series(sorted(set_t1 - set_driver)).to_csv("/content/venn_t1_only.txt", index=False, header=False)
```
<img width="476" height="458" alt="image" src="https://github.com/user-attachments/assets/9140f5e9-4288-45d9-a8d6-06002f86e283" />

# WP019

### Etapa 1
```
%%bash
VEP_VCF="lmabrasil-hg38/vep_output/liftOver_WP019_hg19ToHg38.vep.vcf"

bcftools +split-vep -l $VEP_VCF | \
cut -f2  | \
tr '\n\r' '\t' | \
awk '{print("CHROM\tPOS\tREF\tALT\t"$0"FILTER\tTumorID\tGT\tDP\tAD\tAF\tNormalID\tNGT\tNDP\tNAD\tNAF")}' \
> liftOver_WP019_hg19ToHg38.vep.tsv

bcftools +split-vep \
-f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%FILTER\t[%SAMPLE\t%GT\t%DP\t%AD\t%AF\t]\n' \
-i 'FMT/DP>=20 && FMT/AF>=0.1' -d -A tab $VEP_VCF \
-p x  >> liftOver_WP019_hg19ToHg38.vep.tsv
```

### Etapa 2
```
import pandas as pd
import numpy as np
import re
import requests

# =========================
# CONFIG
# =========================
INPUT_TSV  = "/content/liftOver_WP019_hg19ToHg38.vep.tsv"
OUTPUT_TSV = "/content/WP019-tier.tsv"

DRIVER_RAW_URL = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES
# =========================
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def pick_first_existing_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def parse_percent_vaf(x):
    # retorna VAF em % (0-100). Se vier 0-1, converte para 0-100.
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", ".", "NA", "NaN", "nan", "None"}:
        return np.nan
    s = s.replace(",", ".").replace("%", "")
    try:
        v = float(s)
    except:
        return np.nan
    return v * 100.0 if v <= 1.0 else v

def filter_is_pass(x):
    return str(x).strip().upper() == "PASS"

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    # seu arquivo raw está em UMA linha com genes separados por espaço :contentReference[oaicite:1]{index=1}
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)

# =========================
# LOAD
# =========================
driver_genes = load_driver_genes(DRIVER_RAW_URL)
print("Driver genes carregados:", len(driver_genes))

df = smart_read_tsv(INPUT_TSV)

# detectar colunas (ajuste se necessário)
gene_col     = pick_first_existing_col(df, ["SYMBOL", "Gene", "GENE", "HGNC"])
filter_col   = pick_first_existing_col(df, ["FILTER", "Filter"])
dp_col       = pick_first_existing_col(df, ["DP", "DP_tumor", "TUMOR_DP", "DEPTH"])
vaf_col      = pick_first_existing_col(df, ["VAF_tumor", "VAF", "AF", "TUMOR_AF"])
existing_col = pick_first_existing_col(df, ["Existing_variation", "ExistingVariation", "existing_variation"])

missing = [("SYMBOL", gene_col), ("FILTER", filter_col), ("DP", dp_col), ("VAF/AF", vaf_col), ("Existing_variation", existing_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas necessárias: {missing}\nColunas do arquivo: {list(df.columns)}")

# =========================
# TIERS (hierárquico)
# =========================
df["_GENE_UP"] = df[gene_col].astype(str).str.strip().str.upper()
df["_PASS"]    = df[filter_col].apply(filter_is_pass)
df["_DP"]      = pd.to_numeric(df[dp_col].astype(str).str.replace(",", ".", regex=False), errors="coerce").fillna(0).astype(int)
df["_VAF_PCT"] = df[vaf_col].apply(parse_percent_vaf)
df["_EXIST"]   = df[existing_col].astype(str)

base  = (df["_PASS"]) & (df["_DP"] >= 20) & (df["_VAF_PCT"] >= 10)

tier1 = base & (df["_GENE_UP"].isin(driver_genes))
tier2 = base & (~tier1) & (df["_EXIST"].str.contains(r"\bCOSV\b", flags=re.IGNORECASE, regex=True))
tier3 = ~(tier1 | tier2)

df["Tier"] = np.select([tier1, tier2, tier3], ["Tier 1", "Tier 2", "Tier 3"], default="Tier 3")

# (opcional) auditoria rápida
df["Tier_rule_PASS"]    = df["_PASS"]
df["Tier_rule_DP"]      = df["_DP"]
df["Tier_rule_VAF_pct"] = df["_VAF_PCT"]

df.drop(columns=[c for c in df.columns if c.startswith("_")], inplace=True)
df.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("OK ->", OUTPUT_TSV)
print(df["Tier"].value_counts(dropna=False))
```

### Output:
```
Driver genes carregados: 64
OK -> /content/WP019-tier.tsv
Tier
Tier 3    2186
Tier 1       3
Name: count, dtype: int64
```

### Etapa 3
```
# =========================
# IMPORTS
# =========================
import pandas as pd
import numpy as np
import requests

# =========================
# ARQUIVOS
# =========================
input_tsv  = "/content/liftOver_WP019_hg19ToHg38.vep.tsv"
output_tsv = "/content/WP019_tier_table.tsv"

driver_url = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES (Re-definir ou garantir que smart_read_tsv esteja disponível)
# =========================
# A função smart_read_tsv já está definida na célula anterior (uzRA_DALNTe5)
# e pode ser reutilizada aqui. Se esta célula for executada isoladamente,
# a função precisaria ser copiada ou importada.
# Por simplicidade e dado o contexto do notebook, vamos assumir que está disponível.

# Função smart_read_tsv (copiada para garantir disponibilidade se a célula anterior não for executada)
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)


# =========================
# LER VEP
# =========================
df = smart_read_tsv(input_tsv)

# =========================
# GENES DRIVER
# =========================
driver_genes = load_driver_genes(driver_url)

df
```

#### ↳ O Output esperado é uma longa tabela e você encontra ela no Colab. 

### Etapa 4
```
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import zipfile
import os

# Define the path to the zip file
ZIP_ARCHIVE_PATH = "/content/cgi_somaticmfwp019_results.zip"
# Define the directory where content will be extracted
EXTRACT_TO_DIR = "/content/extracted_cgi_data"

# Ensure the extraction directory exists
os.makedirs(EXTRACT_TO_DIR, exist_ok=True)

extracted_alterations_tsv_path = None
try:
    with zipfile.ZipFile(ZIP_ARCHIVE_PATH, 'r') as zip_ref:
        # Search for 'alterations.tsv' within the zip file.
        # This handles cases where it might be directly at root or nested in a folder.
        for member in zip_ref.namelist():
            if member.endswith('alterations.tsv'):
                zip_ref.extract(member, EXTRACT_TO_DIR)
                extracted_alterations_tsv_path = os.path.join(EXTRACT_TO_DIR, member)
                break

        if not extracted_alterations_tsv_path:
            raise FileNotFoundError(f"Could not find 'alterations.tsv' within {ZIP_ARCHIVE_PATH}")

    # Update ALTERATIONS_TSV to point to the newly extracted file
    ALTERATIONS_TSV = extracted_alterations_tsv_path
    print(f"'{os.path.basename(ALTERATIONS_TSV)}' extracted and path updated to: {ALTERATIONS_TSV}")

except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the zip file exists and contains 'alterations.tsv'.")
    raise
except zipfile.BadZipFile:
    print(f"Error: '{ZIP_ARCHIVE_PATH}' is not a valid zip file.")
    raise
except Exception as e:
    print(f"An unexpected error occurred during zip extraction: {e}")
    raise

TIERS_TSV       = "/content/WP019-tier.tsv"    # seu TSV com coluna Tier (Tier 1/2/3)

def pick_first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def norm_chr(x):
    s = str(x).strip()
    if not s:
        return s
    return s if s.lower().startswith("chr") else "chr" + s

def make_var_id(df, chr_col, pos_col, ref_col, alt_col):
    return (
        df[chr_col].astype(str).map(norm_chr) + ":" +
        df[pos_col].astype(str) + ":" +
        df[ref_col].astype(str) + ":" +
        df[alt_col].astype(str)
    )

# ----------------------------
# 1) DRIVER set (a partir do alterations.tsv)
# ----------------------------
alt = pd.read_csv(ALTERATIONS_TSV, sep="\t", dtype=str)

chr_col = pick_first_existing_col(alt, ["CHR", "CHROM", "#CHROM", "CHROMOSOME"])
pos_col = pick_first_existing_col(alt, ["POS", "POSITION"])
ref_col = pick_first_existing_col(alt, ["REF"])
alt_col = pick_first_existing_col(alt, ["ALT"])
pred_col = pick_first_existing_col(alt, ["CGI-Oncogenic Prediction", "CGI-Oncogenic Summary"])

missing = [("CHR", chr_col), ("POS", pos_col), ("REF", ref_col), ("ALT", alt_col), ("DriverColumn", pred_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no alterations.tsv: {missing}\nColunas: {list(alt.columns)}")

alt["VAR_ID"] = make_var_id(alt, chr_col, pos_col, ref_col, alt_col)

# considera driver se a coluna CGI tiver a palavra "driver"
is_driver = alt[pred_col].fillna("").astype(str).str.contains(r"\bdriver\b", flags=re.IGNORECASE, regex=True)
set_driver = set(alt.loc[is_driver, "VAR_ID"])

# ----------------------------
# 2) T1 set (a partir do tiers TSV)
# ----------------------------
tiers = pd.read_csv(TIERS_TSV, sep="\t", dtype=str)

t_chr = pick_first_existing_col(tiers, ["CHROM", "#CHROM", "CHR", "Chrom"])
t_pos = pick_first_existing_col(tiers, ["POS", "Position", "pos"])
t_ref = pick_first_existing_col(tiers, ["REF", "Ref"])
t_alt = pick_first_existing_col(tiers, ["ALT", "Alt"])

missing = [("CHROM", t_chr), ("POS", t_pos), ("REF", t_ref), ("ALT", t_alt)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no tiers TSV: {missing}\nColunas: {list(tiers.columns)}")

tiers["VAR_ID"] = make_var_id(tiers, t_chr, t_pos, t_ref, t_alt)
set_t1 = set(tiers.loc[tiers["Tier"].astype(str) == "Tier 1", "VAR_ID"])

# ----------------------------
# 3) Venn
# ----------------------------
plt.figure(figsize=(6, 6))
venn2([set_driver, set_t1], set_labels=("Driver (alterations.tsv)", "Tier 1 (tiers.tsv)"))
plt.title("Venn: chr:pos:ref:alt — Driver vs Tier 1")
plt.show()

print("N Driver:", len(set_driver))
print("N Tier 1:", len(set_t1))
print("Interseção:", len(set_driver & set_t1))
print("Driver somente:", len(set_driver - set_t1))
print("T1 somente:", len(set_t1 - set_driver))

# (Opcional) salvar listas
pd.Series(sorted(set_driver & set_t1)).to_csv("/content/venn_intersection.txt", index=False, header=False)
pd.Series(sorted(set_driver - set_t1)).to_csv("/content/venn_driver_only.txt", index=False, header=False)
pd.Series(sorted(set_t1 - set_driver)).to_csv("/content/venn_t1_only.txt", index=False, header=False)
```

<img width="436" height="476" alt="image" src="https://github.com/user-attachments/assets/32becc77-fa70-4219-ab7f-d317aaf77d98" />

# WP058

### Etapa 1
```
%%bash
VEP_VCF="lmabrasil-hg38/vep_output/liftOver_WP058_hg19ToHg38.vep.vcf"

bcftools +split-vep -l $VEP_VCF | \
cut -f2  | \
tr '\n\r' '\t' | \
awk '{print("CHROM\tPOS\tREF\tALT\t"$0"FILTER\tTumorID\tGT\tDP\tAD\tAF\tNormalID\tNGT\tNDP\tNAD\tNAF")}' \
> liftOver_WP058_hg19ToHg38.vep.tsv

bcftools +split-vep \
-f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%FILTER\t[%SAMPLE\t%GT\t%DP\t%AD\t%AF\t]\n' \
-i 'FMT/DP>=20 && FMT/AF>=0.1' -d -A tab $VEP_VCF \
-p x  >> liftOver_WP058_hg19ToHg38.vep.tsv
```

### Etapa 2
```
import pandas as pd
import numpy as np
import re
import requests

# =========================
# CONFIG
# =========================
INPUT_TSV  = "/content/liftOver_WP058_hg19ToHg38.vep.tsv"
OUTPUT_TSV = "/content/WP058-tier.tsv"

DRIVER_RAW_URL = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES
# =========================
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def pick_first_existing_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def parse_percent_vaf(x):
    # retorna VAF em % (0-100). Se vier 0-1, converte para 0-100.
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", ".", "NA", "NaN", "nan", "None"}:
        return np.nan
    s = s.replace(",", ".").replace("%", "")
    try:
        v = float(s)
    except:
        return np.nan
    return v * 100.0 if v <= 1.0 else v

def filter_is_pass(x):
    return str(x).strip().upper() == "PASS"

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    # seu arquivo raw está em UMA linha com genes separados por espaço :contentReference[oaicite:1]{index=1}
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)

# =========================
# LOAD
# =========================
driver_genes = load_driver_genes(DRIVER_RAW_URL)
print("Driver genes carregados:", len(driver_genes))

df = smart_read_tsv(INPUT_TSV)

# detectar colunas (ajuste se necessário)
gene_col     = pick_first_existing_col(df, ["SYMBOL", "Gene", "GENE", "HGNC"])
filter_col   = pick_first_existing_col(df, ["FILTER", "Filter"])
dp_col       = pick_first_existing_col(df, ["DP", "DP_tumor", "TUMOR_DP", "DEPTH"])
vaf_col      = pick_first_existing_col(df, ["VAF_tumor", "VAF", "AF", "TUMOR_AF"])
existing_col = pick_first_existing_col(df, ["Existing_variation", "ExistingVariation", "existing_variation"])

missing = [("SYMBOL", gene_col), ("FILTER", filter_col), ("DP", dp_col), ("VAF/AF", vaf_col), ("Existing_variation", existing_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas necessárias: {missing}\nColunas do arquivo: {list(df.columns)}")

# =========================
# TIERS (hierárquico)
# =========================
df["_GENE_UP"] = df[gene_col].astype(str).str.strip().str.upper()
df["_PASS"]    = df[filter_col].apply(filter_is_pass)
df["_DP"]      = pd.to_numeric(df[dp_col].astype(str).str.replace(",", ".", regex=False), errors="coerce").fillna(0).astype(int)
df["_VAF_PCT"] = df[vaf_col].apply(parse_percent_vaf)
df["_EXIST"]   = df[existing_col].astype(str)

base  = (df["_PASS"]) & (df["_DP"] >= 20) & (df["_VAF_PCT"] >= 10)

tier1 = base & (df["_GENE_UP"].isin(driver_genes))
tier2 = base & (~tier1) & (df["_EXIST"].str.contains(r"\bCOSV\b", flags=re.IGNORECASE, regex=True))
tier3 = ~(tier1 | tier2)

df["Tier"] = np.select([tier1, tier2, tier3], ["Tier 1", "Tier 2", "Tier 3"], default="Tier 3")

# (opcional) auditoria rápida
df["Tier_rule_PASS"]    = df["_PASS"]
df["Tier_rule_DP"]      = df["_DP"]
df["Tier_rule_VAF_pct"] = df["_VAF_PCT"]

df.drop(columns=[c for c in df.columns if c.startswith("_")], inplace=True)
df.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("OK ->", OUTPUT_TSV)
print(df["Tier"].value_counts(dropna=False))
```

### Output:
```
Driver genes carregados: 64
OK -> /content/WP058-tier.tsv
Tier
Tier 3    1331
Name: count, dtype: int64
```

### Etapa 3
```
# =========================
# IMPORTS
# =========================
import pandas as pd
import numpy as np
import requests

# =========================
# ARQUIVOS
# =========================
input_tsv  = "/content/liftOver_WP058_hg19ToHg38.vep.tsv"
output_tsv = "/content/WP058_tier_table.tsv"

driver_url = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES (Re-definir ou garantir que smart_read_tsv esteja disponível)
# =========================
# A função smart_read_tsv já está definida na célula anterior (uzRA_DALNTe5)
# e pode ser reutilizada aqui. Se esta célula for executada isoladamente,
# a função precisaria ser copiada ou importada.
# Por simplicidade e dado o contexto do notebook, vamos assumir que está disponível.

# Função smart_read_tsv (copiada para garantir disponibilidade se a célula anterior não for executada)
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)


# =========================
# LER VEP
# =========================
df = smart_read_tsv(input_tsv)

# =========================
# GENES DRIVER
# =========================
driver_genes = load_driver_genes(driver_url)

df
```

#### ↳ O Output esperado é uma longa tabela e você encontra ela no Colab. 

### Etapa 4
```
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import zipfile
import os

# Define the path to the zip file
ZIP_ARCHIVE_PATH = "/content/cgi_somaticmfwp058_results.zip"
# Define the directory where content will be extracted
EXTRACT_TO_DIR = "/content/extracted_cgi_data"

# Ensure the extraction directory exists
os.makedirs(EXTRACT_TO_DIR, exist_ok=True)

extracted_alterations_tsv_path = None
try:
    with zipfile.ZipFile(ZIP_ARCHIVE_PATH, 'r') as zip_ref:
        # Search for 'alterations.tsv' within the zip file.
        # This handles cases where it might be directly at root or nested in a folder.
        for member in zip_ref.namelist():
            if member.endswith('alterations.tsv'):
                zip_ref.extract(member, EXTRACT_TO_DIR)
                extracted_alterations_tsv_path = os.path.join(EXTRACT_TO_DIR, member)
                break

        if not extracted_alterations_tsv_path:
            raise FileNotFoundError(f"Could not find 'alterations.tsv' within {ZIP_ARCHIVE_PATH}")

    # Update ALTERATIONS_TSV to point to the newly extracted file
    ALTERATIONS_TSV = extracted_alterations_tsv_path
    print(f"'{os.path.basename(ALTERATIONS_TSV)}' extracted and path updated to: {ALTERATIONS_TSV}")

except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the zip file exists and contains 'alterations.tsv'.")
    raise
except zipfile.BadZipFile:
    print(f"Error: '{ZIP_ARCHIVE_PATH}' is not a valid zip file.")
    raise
except Exception as e:
    print(f"An unexpected error occurred during zip extraction: {e}")
    raise

TIERS_TSV       = "/content/WP058-tier.tsv"    # seu TSV com coluna Tier (Tier 1/2/3)

def pick_first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def norm_chr(x):
    s = str(x).strip()
    if not s:
        return s
    return s if s.lower().startswith("chr") else "chr" + s

def make_var_id(df, chr_col, pos_col, ref_col, alt_col):
    return (
        df[chr_col].astype(str).map(norm_chr) + ":" +
        df[pos_col].astype(str) + ":" +
        df[ref_col].astype(str) + ":" +
        df[alt_col].astype(str)
    )

# ----------------------------
# 1) DRIVER set (a partir do alterations.tsv)
# ----------------------------
alt = pd.read_csv(ALTERATIONS_TSV, sep="\t", dtype=str)

chr_col = pick_first_existing_col(alt, ["CHR", "CHROM", "#CHROM", "CHROMOSOME"])
pos_col = pick_first_existing_col(alt, ["POS", "POSITION"])
ref_col = pick_first_existing_col(alt, ["REF"])
alt_col = pick_first_existing_col(alt, ["ALT"])
pred_col = pick_first_existing_col(alt, ["CGI-Oncogenic Prediction", "CGI-Oncogenic Summary"])

missing = [("CHR", chr_col), ("POS", pos_col), ("REF", ref_col), ("ALT", alt_col), ("DriverColumn", pred_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no alterations.tsv: {missing}\nColunas: {list(alt.columns)}")

alt["VAR_ID"] = make_var_id(alt, chr_col, pos_col, ref_col, alt_col)

# considera driver se a coluna CGI tiver a palavra "driver"
is_driver = alt[pred_col].fillna("").astype(str).str.contains(r"\bdriver\b", flags=re.IGNORECASE, regex=True)
set_driver = set(alt.loc[is_driver, "VAR_ID"])

# ----------------------------
# 2) T1 set (a partir do tiers TSV)
# ----------------------------
tiers = pd.read_csv(TIERS_TSV, sep="\t", dtype=str)

t_chr = pick_first_existing_col(tiers, ["CHROM", "#CHROM", "CHR", "Chrom"])
t_pos = pick_first_existing_col(tiers, ["POS", "Position", "pos"])
t_ref = pick_first_existing_col(tiers, ["REF", "Ref"])
t_alt = pick_first_existing_col(tiers, ["ALT", "Alt"])

missing = [("CHROM", t_chr), ("POS", t_pos), ("REF", t_ref), ("ALT", t_alt)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no tiers TSV: {missing}\nColunas: {list(tiers.columns)}")

tiers["VAR_ID"] = make_var_id(tiers, t_chr, t_pos, t_ref, t_alt)
set_t1 = set(tiers.loc[tiers["Tier"].astype(str) == "Tier 1", "VAR_ID"])

# ----------------------------
# 3) Venn
# ----------------------------
plt.figure(figsize=(6, 6))
venn2([set_driver, set_t1], set_labels=("Driver (alterations.tsv)", "Tier 1 (tiers.tsv)"))
plt.title("Venn: chr:pos:ref:alt — Driver vs Tier 1")
plt.show()

print("N Driver:", len(set_driver))
print("N Tier 1:", len(set_t1))
print("Interseção:", len(set_driver & set_t1))
print("Driver somente:", len(set_driver - set_t1))
print("T1 somente:", len(set_t1 - set_driver))

# (Opcional) salvar listas
pd.Series(sorted(set_driver & set_t1)).to_csv("/content/venn_intersection.txt", index=False, header=False)
pd.Series(sorted(set_driver - set_t1)).to_csv("/content/venn_driver_only.txt", index=False, header=False)
pd.Series(sorted(set_t1 - set_driver)).to_csv("/content/venn_t1_only.txt", index=False, header=False)
```

<img width="504" height="410" alt="image" src="https://github.com/user-attachments/assets/a4cc5994-35a7-4fa5-985c-b852570c24d7" />

# WP068

### Etapa 1
```
%%bash
VEP_VCF="lmabrasil-hg38/vep_output/liftOver_WP068_hg19ToHg38.vep.vcf"

bcftools +split-vep -l $VEP_VCF | \
cut -f2  | \
tr '\n\r' '\t' | \
awk '{print("CHROM\tPOS\tREF\tALT\t"$0"FILTER\tTumorID\tGT\tDP\tAD\tAF\tNormalID\tNGT\tNDP\tNAD\tNAF")}' \
> liftOver_WP068_hg19ToHg38.vep.tsv

bcftools +split-vep \
-f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%FILTER\t[%SAMPLE\t%GT\t%DP\t%AD\t%AF\t]\n' \
-i 'FMT/DP>=20 && FMT/AF>=0.1' -d -A tab $VEP_VCF \
-p x  >> liftOver_WP068_hg19ToHg38.vep.tsv
```

### Etapa 2
```
import pandas as pd
import numpy as np
import re
import requests

# =========================
# CONFIG
# =========================
INPUT_TSV  = "/content/liftOver_WP068_hg19ToHg38.vep.tsv"
OUTPUT_TSV = "/content/WP068-tier.tsv"

DRIVER_RAW_URL = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES
# =========================
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def pick_first_existing_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def parse_percent_vaf(x):
    # retorna VAF em % (0-100). Se vier 0-1, converte para 0-100.
    if x is None:
        return np.nan
    s = str(x).strip()
    if s in {"", ".", "NA", "NaN", "nan", "None"}:
        return np.nan
    s = s.replace(",", ".").replace("%", "")
    try:
        v = float(s)
    except:
        return np.nan
    return v * 100.0 if v <= 1.0 else v

def filter_is_pass(x):
    return str(x).strip().upper() == "PASS"

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    # seu arquivo raw está em UMA linha com genes separados por espaço :contentReference[oaicite:1]{index=1}
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)

# =========================
# LOAD
# =========================
driver_genes = load_driver_genes(DRIVER_RAW_URL)
print("Driver genes carregados:", len(driver_genes))

df = smart_read_tsv(INPUT_TSV)

# detectar colunas (ajuste se necessário)
gene_col     = pick_first_existing_col(df, ["SYMBOL", "Gene", "GENE", "HGNC"])
filter_col   = pick_first_existing_col(df, ["FILTER", "Filter"])
dp_col       = pick_first_existing_col(df, ["DP", "DP_tumor", "TUMOR_DP", "DEPTH"])
vaf_col      = pick_first_existing_col(df, ["VAF_tumor", "VAF", "AF", "TUMOR_AF"])
existing_col = pick_first_existing_col(df, ["Existing_variation", "ExistingVariation", "existing_variation"])

missing = [("SYMBOL", gene_col), ("FILTER", filter_col), ("DP", dp_col), ("VAF/AF", vaf_col), ("Existing_variation", existing_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas necessárias: {missing}\nColunas do arquivo: {list(df.columns)}")

# =========================
# TIERS (hierárquico)
# =========================
df["_GENE_UP"] = df[gene_col].astype(str).str.strip().str.upper()
df["_PASS"]    = df[filter_col].apply(filter_is_pass)
df["_DP"]      = pd.to_numeric(df[dp_col].astype(str).str.replace(",", ".", regex=False), errors="coerce").fillna(0).astype(int)
df["_VAF_PCT"] = df[vaf_col].apply(parse_percent_vaf)
df["_EXIST"]   = df[existing_col].astype(str)

base  = (df["_PASS"]) & (df["_DP"] >= 20) & (df["_VAF_PCT"] >= 10)

tier1 = base & (df["_GENE_UP"].isin(driver_genes))
tier2 = base & (~tier1) & (df["_EXIST"].str.contains(r"\bCOSV\b", flags=re.IGNORECASE, regex=True))
tier3 = ~(tier1 | tier2)

df["Tier"] = np.select([tier1, tier2, tier3], ["Tier 1", "Tier 2", "Tier 3"], default="Tier 3")

# (opcional) auditoria rápida
df["Tier_rule_PASS"]    = df["_PASS"]
df["Tier_rule_DP"]      = df["_DP"]
df["Tier_rule_VAF_pct"] = df["_VAF_PCT"]

df.drop(columns=[c for c in df.columns if c.startswith("_")], inplace=True)
df.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("OK ->", OUTPUT_TSV)
print(df["Tier"].value_counts(dropna=False))
```

### Output:
```
Driver genes carregados: 64
OK -> /content/WP068-tier.tsv
Tier
Tier 3    1222
Tier 1       1
Name: count, dtype: int64
```

### Etapa 3
```
# =========================
# IMPORTS
# =========================
import pandas as pd
import numpy as np
import requests

# =========================
# ARQUIVOS
# =========================
input_tsv  = "/content/liftOver_WP068_hg19ToHg38.vep.tsv"
output_tsv = "/content/WP068_tier_table.tsv"

driver_url = "https://raw.githubusercontent.com/renatopuga/lmabrasil-hg38/refs/heads/main/hpo/Clonal_Hematopoiesis_driver_genes.txt"

# =========================
# FUNÇÕES (Re-definir ou garantir que smart_read_tsv esteja disponível)
# =========================
# A função smart_read_tsv já está definida na célula anterior (uzRA_DALNTe5)
# e pode ser reutilizada aqui. Se esta célula for executada isoladamente,
# a função precisaria ser copiada ou importada.
# Por simplicidade e dado o contexto do notebook, vamos assumir que está disponível.

# Função smart_read_tsv (copiada para garantir disponibilidade se a célula anterior não for executada)
def smart_read_tsv(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        n = len(cols)

        rows = []
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) > n:
                parts = parts[:n-1] + ["\t".join(parts[n-1:])]
            elif len(parts) < n:
                parts = parts + [""] * (n - len(parts))
            rows.append(parts)

    return pd.DataFrame(rows, columns=cols)

def load_driver_genes(url: str) -> set:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    tokens = re.split(r"[\s,;]+", r.text.strip())
    genes = [t.strip().upper() for t in tokens if t.strip() and not t.startswith("#")]
    return set(genes)


# =========================
# LER VEP
# =========================
df = smart_read_tsv(input_tsv)

# =========================
# GENES DRIVER
# =========================
driver_genes = load_driver_genes(driver_url)

df
```

#### ↳ O Output esperado é uma longa tabela e você encontra ela no Colab. 

### Etapa 4
```
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import zipfile
import os

# Define the path to the zip file
ZIP_ARCHIVE_PATH = "/content/cgi_somaticmfwp068_results.zip"
# Define the directory where content will be extracted
EXTRACT_TO_DIR = "/content/extracted_cgi_data"

# Ensure the extraction directory exists
os.makedirs(EXTRACT_TO_DIR, exist_ok=True)

extracted_alterations_tsv_path = None
try:
    with zipfile.ZipFile(ZIP_ARCHIVE_PATH, 'r') as zip_ref:
        # Search for 'alterations.tsv' within the zip file.
        # This handles cases where it might be directly at root or nested in a folder.
        for member in zip_ref.namelist():
            if member.endswith('alterations.tsv'):
                zip_ref.extract(member, EXTRACT_TO_DIR)
                extracted_alterations_tsv_path = os.path.join(EXTRACT_TO_DIR, member)
                break

        if not extracted_alterations_tsv_path:
            raise FileNotFoundError(f"Could not find 'alterations.tsv' within {ZIP_ARCHIVE_PATH}")

    # Update ALTERATIONS_TSV to point to the newly extracted file
    ALTERATIONS_TSV = extracted_alterations_tsv_path
    print(f"'{os.path.basename(ALTERATIONS_TSV)}' extracted and path updated to: {ALTERATIONS_TSV}")

except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the zip file exists and contains 'alterations.tsv'.")
    raise
except zipfile.BadZipFile:
    print(f"Error: '{ZIP_ARCHIVE_PATH}' is not a valid zip file.")
    raise
except Exception as e:
    print(f"An unexpected error occurred during zip extraction: {e}")
    raise

TIERS_TSV       = "/content/WP068-tier.tsv"    # seu TSV com coluna Tier (Tier 1/2/3)

def pick_first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def norm_chr(x):
    s = str(x).strip()
    if not s:
        return s
    return s if s.lower().startswith("chr") else "chr" + s

def make_var_id(df, chr_col, pos_col, ref_col, alt_col):
    return (
        df[chr_col].astype(str).map(norm_chr) + ":" +
        df[pos_col].astype(str) + ":" +
        df[ref_col].astype(str) + ":" +
        df[alt_col].astype(str)
    )

# ----------------------------
# 1) DRIVER set (a partir do alterations.tsv)
# ----------------------------
alt = pd.read_csv(ALTERATIONS_TSV, sep="\t", dtype=str)

chr_col = pick_first_existing_col(alt, ["CHR", "CHROM", "#CHROM", "CHROMOSOME"])
pos_col = pick_first_existing_col(alt, ["POS", "POSITION"])
ref_col = pick_first_existing_col(alt, ["REF"])
alt_col = pick_first_existing_col(alt, ["ALT"])
pred_col = pick_first_existing_col(alt, ["CGI-Oncogenic Prediction", "CGI-Oncogenic Summary"])

missing = [("CHR", chr_col), ("POS", pos_col), ("REF", ref_col), ("ALT", alt_col), ("DriverColumn", pred_col)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no alterations.tsv: {missing}\nColunas: {list(alt.columns)}")

alt["VAR_ID"] = make_var_id(alt, chr_col, pos_col, ref_col, alt_col)

# considera driver se a coluna CGI tiver a palavra "driver"
is_driver = alt[pred_col].fillna("").astype(str).str.contains(r"\bdriver\b", flags=re.IGNORECASE, regex=True)
set_driver = set(alt.loc[is_driver, "VAR_ID"])

# ----------------------------
# 2) T1 set (a partir do tiers TSV)
# ----------------------------
tiers = pd.read_csv(TIERS_TSV, sep="\t", dtype=str)

t_chr = pick_first_existing_col(tiers, ["CHROM", "#CHROM", "CHR", "Chrom"])
t_pos = pick_first_existing_col(tiers, ["POS", "Position", "pos"])
t_ref = pick_first_existing_col(tiers, ["REF", "Ref"])
t_alt = pick_first_existing_col(tiers, ["ALT", "Alt"])

missing = [("CHROM", t_chr), ("POS", t_pos), ("REF", t_ref), ("ALT", t_alt)]
missing = [name for name, col in missing if col is None]
if missing:
    raise ValueError(f"Faltando colunas no tiers TSV: {missing}\nColunas: {list(tiers.columns)}")

tiers["VAR_ID"] = make_var_id(tiers, t_chr, t_pos, t_ref, t_alt)
set_t1 = set(tiers.loc[tiers["Tier"].astype(str) == "Tier 1", "VAR_ID"])

# ----------------------------
# 3) Venn
# ----------------------------
plt.figure(figsize=(6, 6))
venn2([set_driver, set_t1], set_labels=("Driver (alterations.tsv)", "Tier 1 (tiers.tsv)"))
plt.title("Venn: chr:pos:ref:alt — Driver vs Tier 1")
plt.show()

print("N Driver:", len(set_driver))
print("N Tier 1:", len(set_t1))
print("Interseção:", len(set_driver & set_t1))
print("Driver somente:", len(set_driver - set_t1))
print("T1 somente:", len(set_t1 - set_driver))

# (Opcional) salvar listas
pd.Series(sorted(set_driver & set_t1)).to_csv("/content/venn_intersection.txt", index=False, header=False)
pd.Series(sorted(set_driver - set_t1)).to_csv("/content/venn_driver_only.txt", index=False, header=False)
pd.Series(sorted(set_t1 - set_driver)).to_csv("/content/venn_t1_only.txt", index=False, header=False)
```

<img width="510" height="276" alt="image" src="https://github.com/user-attachments/assets/9ce8e3e9-a2d2-4916-9c0a-67dff74d5107" />
