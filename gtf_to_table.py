import os
import csv
from collections import defaultdict

import pandas as pd
import mysql.connector


# ============================================================
# ======================= CONFIG START =======================
# ============================================================

# --- MySQL connection ---
DB_CONFIG = {
    "host": "localhost",
    "user": "root",
    "password": "",
    "database": "as_db",
}

# --- Paths (EDIT THESE) ---
CANONICAL_GTF_FILE = "/Users/galno/tmp_for_local_run/Arabidopsis_thaliana.TAIR10.60.gtf"
BAM_DIRECTORY = "/Volumes/bam"
WORKING_DIRECTORY = "/Users/galno/as_2"

# --- External tool ---
STRINGTIE_BIN = "stringtie"  # if not in PATH, set full path e.g. "/usr/local/bin/stringtie"

# --- Experiment design ---
GENOTYPES = ["Col-0", "tfiis", "upf1", "upf3", "tfiis-upf1", "tfiis-upf3"]
TREATMENTS = ["NT", "1h", "1d"]
N_REPLICATES = 3

# --- Output / intermediate filenames ---
CANONICAL_CSV_FILENAME = os.path.join(WORKING_DIRECTORY, "from_canonical_gtf.csv")

# --- MySQL table ---
MYSQL_GTF_TABLE = "gtf"

# ============================================================
# ======================== CONFIG END ========================
# ============================================================


def gtf_to_csv(canonical_gtf_file: str, out_csv: str) -> str:
    gene_ids = []
    starts = []
    ends = []
    strands = []

    with open(canonical_gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            if parts[2] == "transcript":
                attributes = parts[8]

                attr_dict = {}
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    try:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.replace('"', "")
                    except ValueError:
                        continue

                gene_id = attr_dict.get("gene_id")
                if gene_id is None:
                    continue

                start_position = int(parts[3])
                end_position = int(parts[4])
                strand = parts[6]

                gene_ids.append(gene_id)
                starts.append(start_position)
                ends.append(end_position)
                strands.append(strand)

    df = pd.DataFrame(
        {"gene": gene_ids, "start_position": starts, "end_position": ends, "strand": strands}
    )
    df.to_csv(out_csv, index=False)
    return out_csv


def load_gene_intervals(csv_path: str):
    genes = []
    with open(csv_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # NOTE: this assumes something like "AT1G01010" and uses row['gene'][2] -> "1"
            genes.append(
                {
                    "gene_id": row["gene"],
                    "chromosome": row["gene"][2],
                    "start": int(row["start_position"]),
                    "end": int(row["end_position"]),
                    "strand": row["strand"],
                }
            )
    return genes


def find_matching_gene(gtf_start, gtf_end, gtf_strand, chromosome, genes):
    gtf_start = int(gtf_start)
    gtf_end = int(gtf_end)

    # 1) Exact match first
    for gene in genes:
        if gene["strand"] != gtf_strand:
            continue
        if gene["chromosome"] != chromosome:
            continue

        gene_start = min(gene["start"], gene["end"])
        gene_end = max(gene["start"], gene["end"])

        if gene_start == gtf_start and gene_end == gtf_end:
            return gene["gene_id"]

    # 2) Otherwise first overlapping hit
    for gene in genes:
        if gene["strand"] != gtf_strand:
            continue
        if gene["chromosome"] != chromosome:
            continue

        gene_start = min(gene["start"], gene["end"])
        gene_end = max(gene["start"], gene["end"])

        if not (gtf_end < gene_start or gtf_start > gene_end):
            return gene["gene_id"]

    return None


def gene_ids_from_csv(gtf_file: str, genes, csv_file: str) -> None:
    gene_ids = []
    strg_ids = []
    strg_transcripts = []

    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            if parts[2] != "transcript":
                continue

            attributes = parts[8]
            attr_dict = {}
            for attr in attributes.split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                try:
                    key, value = attr.split(" ", 1)
                    attr_dict[key] = value.replace('"', "")
                except ValueError:
                    continue

            chromosome = parts[0]
            start_position = parts[3]
            end_position = parts[4]
            strand = parts[6]

            strg_id = attr_dict.get("gene_id")
            strg_transcript_id = attr_dict.get("transcript_id")

            if strand == ".":
                continue
            if not strg_transcript_id or not strg_id:
                continue

            matched_gene_id = find_matching_gene(
                start_position, end_position, strand, chromosome, genes
            )
            if matched_gene_id:
                gene_ids.append(matched_gene_id)
                strg_ids.append(strg_id)
                strg_transcripts.append(strg_transcript_id)

    df = pd.DataFrame(
        {"gene": gene_ids, "strg_gene_id": strg_ids, "strg_transcript_id": strg_transcripts}
    )
    df.to_csv(csv_file, index=False)


def modify_gtf(gtf_file: str, csv_file: str, modified_gtf_file: str) -> None:
    df = pd.read_csv(csv_file)
    if df.empty:
        # still write original (or empty) file to make pipeline consistent
        with open(gtf_file, "r") as src, open(modified_gtf_file, "w") as dst:
            dst.writelines(src.readlines())
        return

    mapping = dict(zip(df["strg_transcript_id"], df["gene"]))

    new_lines = []
    with open(gtf_file, "r") as f:
        for line in f:
            if "#" not in line and "gene_id" in line:
                for strg_transcript_id, gene in mapping.items():
                    if f'transcript_id "{strg_transcript_id}"' in line:
                        strg_gene = strg_transcript_id.rsplit(".", 1)[0]
                        line = line.replace(f'gene_id "{strg_gene}"', f'gene_id "{gene}"')
                        break
            new_lines.append(line)

    with open(modified_gtf_file, "w") as f:
        f.writelines(new_lines)


def parse_gtf_exons(gtf_file: str):
    transcripts = defaultdict(list)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            attr_dict = {}
            for attr in attributes.split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                key, value = attr.split(" ", 1)
                attr_dict[key] = value.replace('"', "")

            transcript_id = attr_dict.get("transcript_id")
            fpkm = attr_dict.get("FPKM")
            gene_id = attr_dict.get("gene_id")

            if feature_type == "transcript" and transcript_id:
                transcripts[transcript_id] = {
                    "gene_id": gene_id,
                    "strand": strand,
                    "fpkm": float(fpkm) if fpkm else None,
                    "exons": [],
                }

            elif feature_type == "exon" and transcript_id:
                # only append if transcript key exists (stringtie gtf should have transcript before exon)
                if transcript_id in transcripts:
                    transcripts[transcript_id]["exons"].append((start, end))

    return transcripts


def insert_gtf_records(sample: str, gtf_file: str):
    gtf_records = []

    transcripts = parse_gtf_exons(gtf_file)
    for transcript_id, data in transcripts.items():
        exon_positions = ";".join([f"{s}-{e}" for s, e in data["exons"]])
        gene_id = data["gene_id"]
        uid = f"{sample}_{transcript_id}"

        gtf_records.append(
            (
                uid,
                gene_id,
                sample,
                data["strand"],
                data["fpkm"],
                exon_positions,
            )
        )

    return gtf_records


def build_bam_list(genotype: str, treatment: str) -> str:
    bams = [
        os.path.join(BAM_DIRECTORY, f"{genotype}_{treatment}_{i+1}_sorted.bam")
        for i in range(N_REPLICATES)
    ]
    return " ".join(bams)


def run_stringtie(sample: str, gtf_out: str) -> None:
    bam_list = build_bam_list(sample.split("_")[0], sample.split("_")[1])
    cmd = f'{STRINGTIE_BIN} {bam_list} -G "{CANONICAL_GTF_FILE}" -o "{gtf_out}"'
    os.system(cmd)


def main():
    os.makedirs(WORKING_DIRECTORY, exist_ok=True)

    # --- DB connect ---
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()

    # --- canonical gtf -> csv -> genes ---
    if not os.path.exists(CANONICAL_CSV_FILENAME):
        gtf_to_csv(CANONICAL_GTF_FILE, CANONICAL_CSV_FILENAME)

    genes = load_gene_intervals(CANONICAL_CSV_FILENAME)

    for genotype in GENOTYPES:
        for treatment in TREATMENTS:
            sample = f"{genotype}_{treatment}"

            # 1) Create gtf from BAMs
            gtf_file = os.path.join(WORKING_DIRECTORY, f"{sample}.gtf")
            if not os.path.exists(gtf_file):
                print(f"generate gtf for {sample}")
                bam_list = " ".join(
                    [
                        os.path.join(BAM_DIRECTORY, f"{genotype}_{treatment}_{n+1}_sorted.bam")
                        for n in range(N_REPLICATES)
                    ]
                )
                os.system(
                    f'{STRINGTIE_BIN} {bam_list} -G "{CANONICAL_GTF_FILE}" -o "{gtf_file}"'
                )
            else:
                print(f"skipping {sample}")

            # 2) Create mapping csv
            csv_file = os.path.join(WORKING_DIRECTORY, f"{sample}.csv")
            gene_ids_from_csv(gtf_file, genes, csv_file)

            # 3) Modify gtf
            print(f"modify gtf file: {gtf_file}")
            modified_gtf_file = os.path.join(WORKING_DIRECTORY, f"{sample}_modified.gtf")
            modify_gtf(gtf_file, csv_file, modified_gtf_file)

            # 4) Insert to mysql
            print("insert into mysql...")
            gtf_records = insert_gtf_records(sample, modified_gtf_file)

            insert_query = f"""
                INSERT IGNORE INTO {MYSQL_GTF_TABLE}
                (id, gene_id, sample, strand, fpkm, exon_positions)
                VALUES (%s, %s, %s, %s, %s, %s)
            """

            for record in gtf_records:
                cursor.execute(insert_query, record)

            conn.commit()
            print(f"{sample} is done")

    cursor.close()
    conn.close()


if __name__ == "__main__":
    main()
