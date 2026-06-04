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


def _parse_gtf_attributes(raw: str) -> dict:
    attr_dict = {}
    for attr in raw.split(";"):
        attr = attr.strip()
        if not attr:
            continue
        try:
            key, value = attr.split(" ", 1)
            attr_dict[key] = value.replace('"', "")
        except ValueError:
            continue
    return attr_dict


class GeneIndex:
    """Coordinate-based lookup index built from the canonical TAIR10 GTF.

    Genes are pre-bucketed by (chromosome, strand) so each lookup only
    searches the relevant subset instead of the full gene list.
    """

    def __init__(self, canonical_gtf_file: str, cache_csv: str) -> None:
        if not os.path.exists(cache_csv):
            self._parse_and_cache(canonical_gtf_file, cache_csv)
        self._index: dict[tuple, list] = defaultdict(list)
        self._load(cache_csv)

    def _parse_and_cache(self, gtf_file: str, out_csv: str) -> None:
        gene_ids, starts, ends, strands = [], [], [], []

        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "transcript":
                    continue

                attrs = _parse_gtf_attributes(parts[8])
                gene_id = attrs.get("gene_id")
                if gene_id is None:
                    continue

                gene_ids.append(gene_id)
                starts.append(int(parts[3]))
                ends.append(int(parts[4]))
                strands.append(parts[6])

        pd.DataFrame(
            {"gene": gene_ids, "start_position": starts, "end_position": ends, "strand": strands}
        ).to_csv(out_csv, index=False)

    def _load(self, csv_path: str) -> None:
        with open(csv_path, newline="") as f:
            for row in csv.DictReader(f):
                gene_id = row["gene"]
                # AT1G01010 → chromosome is the character at index 2
                key = (gene_id[2], row["strand"])
                self._index[key].append(
                    {
                        "gene_id": gene_id,
                        "start": int(row["start_position"]),
                        "end": int(row["end_position"]),
                    }
                )

    def find(self, start: int, end: int, strand: str, chromosome: str) -> str | None:
        start, end = int(start), int(end)
        candidates = self._index.get((chromosome, strand), [])

        # 1) Exact match
        for gene in candidates:
            g_start = min(gene["start"], gene["end"])
            g_end = max(gene["start"], gene["end"])
            if g_start == start and g_end == end:
                return gene["gene_id"]

        # 2) First overlapping hit
        for gene in candidates:
            g_start = min(gene["start"], gene["end"])
            g_end = max(gene["start"], gene["end"])
            if not (end < g_start or start > g_end):
                return gene["gene_id"]

        return None



class GTFProcessor:
    """Handles all per-sample GTF operations: StringTie → mapping → modify → DB records."""

    def __init__(self, genotype: str, treatment: str, gene_index: GeneIndex) -> None:
        self.sample = f"{genotype}_{treatment}"
        self.genotype = genotype
        self.treatment = treatment
        self.gene_index = gene_index

        self.gtf_file = os.path.join(WORKING_DIRECTORY, f"{self.sample}.gtf")
        self.csv_file = os.path.join(WORKING_DIRECTORY, f"{self.sample}.csv")
        self.modified_gtf_file = os.path.join(WORKING_DIRECTORY, f"{self.sample}_modified.gtf")

    def run_stringtie(self) -> None:
        if os.path.exists(self.gtf_file):
            print(f"skipping {self.sample}")
            return

        print(f"generate gtf for {self.sample}")
        bams = " ".join(
            os.path.join(BAM_DIRECTORY, f"{self.genotype}_{self.treatment}_{n+1}_sorted.bam")
            for n in range(N_REPLICATES)
        )
        os.system(f'{STRINGTIE_BIN} {bams} -G "{CANONICAL_GTF_FILE}" -o "{self.gtf_file}"')

    def build_mapping_csv(self) -> None:
        gene_ids, strg_ids, strg_transcripts = [], [], []

        with open(self.gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "transcript":
                    continue

                attrs = _parse_gtf_attributes(parts[8])
                strand = parts[6]
                strg_id = attrs.get("gene_id")
                strg_transcript_id = attrs.get("transcript_id")

                if strand == "." or not strg_transcript_id or not strg_id:
                    continue

                matched = self.gene_index.find(parts[3], parts[4], strand, parts[0])
                if matched:
                    gene_ids.append(matched)
                    strg_ids.append(strg_id)
                    strg_transcripts.append(strg_transcript_id)

        pd.DataFrame(
            {"gene": gene_ids, "strg_gene_id": strg_ids, "strg_transcript_id": strg_transcripts}
        ).to_csv(self.csv_file, index=False)

    def write_modified_gtf(self) -> None:
        df = pd.read_csv(self.csv_file)
        if df.empty:
            with open(self.gtf_file) as src, open(self.modified_gtf_file, "w") as dst:
                dst.writelines(src.readlines())
            return

        mapping = dict(zip(df["strg_transcript_id"], df["gene"]))

        new_lines = []
        with open(self.gtf_file) as f:
            for line in f:
                if "#" not in line and "gene_id" in line:
                    for strg_transcript_id, gene in mapping.items():
                        if f'transcript_id "{strg_transcript_id}"' in line:
                            strg_gene = strg_transcript_id.rsplit(".", 1)[0]
                            line = line.replace(f'gene_id "{strg_gene}"', f'gene_id "{gene}"')
                            break
                new_lines.append(line)

        with open(self.modified_gtf_file, "w") as f:
            f.writelines(new_lines)

    def to_db_records(self) -> list[tuple]:
        transcripts: dict = defaultdict(dict)

        with open(self.modified_gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                attrs = _parse_gtf_attributes(parts[8])
                transcript_id = attrs.get("transcript_id")
                feature_type = parts[2]

                if feature_type == "transcript" and transcript_id:
                    fpkm = attrs.get("FPKM")
                    transcripts[transcript_id] = {
                        "gene_id": attrs.get("gene_id"),
                        "strand": parts[6],
                        "fpkm": float(fpkm) if fpkm else None,
                        "exons": [],
                    }
                elif feature_type == "exon" and transcript_id in transcripts:
                    transcripts[transcript_id]["exons"].append((int(parts[3]), int(parts[4])))

        records = []
        for transcript_id, data in transcripts.items():
            exon_positions = ";".join(f"{s}-{e}" for s, e in data["exons"])
            records.append((
                f"{self.sample}_{transcript_id}",
                data["gene_id"],
                self.sample,
                data["strand"],
                data["fpkm"],
                exon_positions,
            ))
        return records


class Pipeline:
    """Orchestrates the full analysis: StringTie → GTF mapping → MySQL insert."""

    INSERT_QUERY = f"""
        INSERT IGNORE INTO {MYSQL_GTF_TABLE}
        (id, gene_id, sample, strand, fpkm, exon_positions)
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    def __init__(self) -> None:
        os.makedirs(WORKING_DIRECTORY, exist_ok=True)
        self.gene_index = GeneIndex(CANONICAL_GTF_FILE, CANONICAL_CSV_FILENAME)
        self.conn = mysql.connector.connect(**DB_CONFIG)
        self.cursor = self.conn.cursor()

    def run(self) -> None:
        try:
            for genotype in GENOTYPES:
                for treatment in TREATMENTS:
                    self._process_sample(genotype, treatment)
        finally:
            self.cursor.close()
            self.conn.close()

    def _process_sample(self, genotype: str, treatment: str) -> None:
        proc = GTFProcessor(genotype, treatment, self.gene_index)

        proc.run_stringtie()
        proc.build_mapping_csv()

        print(f"modify gtf: {proc.gtf_file}")
        proc.write_modified_gtf()

        print("insert into mysql...")
        for record in proc.to_db_records():
            self.cursor.execute(self.INSERT_QUERY, record)

        self.conn.commit()
        print(f"{proc.sample} done")


if __name__ == "__main__":
    Pipeline().run()
