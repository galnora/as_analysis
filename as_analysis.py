import os
import pandas as pd
import mysql.connector
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


class Config:
    def __init__(self):
        self.canonical_gtf_file = '/Users/galno/tmp_for_local_run/Arabidopsis_thaliana.TAIR10.60.gtf'
        self.bam_directory = '/Volumes/bam'
        self.working_directory = '/Users/galno/as_2'
        self.genotypes = ['Col-0', 'tfiis', 'upf1', 'upf3', 'tfiis-upf1', 'tfiis-upf3']
        self.treatments = ['NT', '1h', '1d']
        self.db_config = {
            'host': 'localhost',
            'user': 'root',
            'password': '',
            'database': 'as_db'
        }


class GTFProcessor:
    def __init__(self, config: Config):
        self.config = config
        self.conn = mysql.connector.connect(**config.db_config)
        self.cursor = self.conn.cursor()
        self.genes = self.load_genes()

    def gtf_to_csv(self, gtf_path, out_path):
        gene_ids, starts, ends, strands = [], [], [], []
        with open(gtf_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split("\t")
                    if parts[2] == 'transcript':
                        attributes = parts[8]
                        attr_dict = {}
                        for attr in attributes.split(";"):
                            attr = attr.strip()
                            if attr:
                                try:
                                    key, value = attr.split(" ", 1)
                                    attr_dict[key] = value.replace('"', '')
                                except ValueError:
                                    continue
                        gene_id = attr_dict.get("gene_id")
                        start_position = int(parts[3])
                        end_position = int(parts[4])
                        strand = parts[6]

                        gene_ids.append(gene_id)
                        starts.append(start_position)
                        ends.append(end_position)
                        strands.append(strand)

        df = pd.DataFrame({'gene': gene_ids, 'start_position': starts, 'end_position': ends, 'strand': strands})
        df.to_csv(out_path, index=False)

    def load_gene_intervals(self, csv_path):
        genes = []
        with open(csv_path) as f:
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                genes.append({
                    'gene_id': row['gene'],
                    'chromosome': row['gene'][2],
                    'start': int(row['start_position']),
                    'end': int(row['end_position']),
                    'strand': row['strand']
                })
        return genes

    def load_genes(self):
        csv_path = os.path.join(self.config.working_directory, 'from_canonical_gtf.csv')
        if not os.path.exists(csv_path):
            self.gtf_to_csv(self.config.canonical_gtf_file, csv_path)
        return self.load_gene_intervals(csv_path)

    def find_matching_gene(self, gtf_start, gtf_end, gtf_strand, chromosome):
        gtf_start = int(gtf_start)
        gtf_end = int(gtf_end)
        for gene in self.genes:
            if gene['strand'] != gtf_strand or gene['chromosome'] != chromosome:
                continue
            gene_start = min(gene['start'], gene['end'])
            gene_end = max(gene['start'], gene['end'])
            if gene_start == gtf_start and gene_end == gtf_end:
                return gene['gene_id']
        for gene in self.genes:
            if gene['strand'] != gtf_strand or gene['chromosome'] != chromosome:
                continue
            gene_start = min(gene['start'], gene['end'])
            gene_end = max(gene['start'], gene['end'])
            if not (gtf_end < gene_start or gtf_start > gene_end):
                return gene['gene_id']
        return None

    def gene_ids_from_csv(self, gtf_file, csv_file):
        gene_ids, strg_ids, strg_transcripts = [], [], []
        with open(gtf_file) as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split("\t")
                    if parts[2] == 'transcript':
                        attributes = parts[8]
                        attr_dict = {}
                        for attr in attributes.split(";"):
                            attr = attr.strip()
                            if attr:
                                try:
                                    key, value = attr.split(" ", 1)
                                    attr_dict[key] = value.replace('"', '')
                                except ValueError:
                                    continue
                        chromosome = parts[0]
                        start_position = parts[3]
                        end_position = parts[4]
                        strg_id = attr_dict.get('gene_id')
                        strg_transcript_id = attr_dict.get('transcript_id')
                        strand = parts[6]
                        if strand != '.':
                            matched_gene_id = self.find_matching_gene(start_position, end_position, strand, chromosome)
                            if matched_gene_id:
                                gene_ids.append(matched_gene_id)
                                strg_ids.append(strg_id)
                                strg_transcripts.append(strg_transcript_id)
        df = pd.DataFrame({'gene': gene_ids, 'strg_gene_id': strg_ids, 'strg_transcript_id': strg_transcripts})
        df.to_csv(csv_file, index=False)

    def modify_gtf(self, gtf_file, csv_file, modified_gtf_file):
        df = pd.read_csv(csv_file)
        mapping = dict(zip(df['strg_transcript_id'], df['gene']))
        new_lines = []
        with open(gtf_file) as f:
            for line in f:
                if not '#' in line and 'gene_id' in line:
                    for tid, gene in mapping.items():
                        if f'transcript_id "{tid}"' in line:
                            strg_gene = tid.rsplit('.', 1)[0]
                            line = line.replace(f'gene_id "{strg_gene}"', f'gene_id "{gene}"')
                            break
                new_lines.append(line)
        with open(modified_gtf_file, 'w') as f:
            f.writelines(new_lines)

    def parse_gtf_exons(self, gtf_file):
        transcripts = defaultdict(lambda: {"exons": []})
        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                start, end = int(parts[3]), int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                attr_dict = {}
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if attr:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.replace('"', '')
                transcript_id = attr_dict.get("transcript_id")
                fpkm = attr_dict.get("FPKM")
                gene_id = attr_dict.get("gene_id")
                if transcript_id:
                    transcripts[transcript_id]["strand"] = strand
                    transcripts[transcript_id]["gene_id"] = gene_id
                    transcripts[transcript_id]["fpkm"] = float(fpkm) if fpkm else None
                    if feature_type == "exon":
                        transcripts[transcript_id]["exons"].append((start, end))
        return transcripts

    def insert_gtf_records(self, sample, gtf_file):
        records = []
        transcripts = self.parse_gtf_exons(gtf_file)
        for tid, data in transcripts.items():
            uid = f"{sample}_{tid}"
            exons = ";".join([f"{s}-{e}" for s, e in data["exons"]])
            records.append((uid, data["gene_id"], sample, data["strand"], data["fpkm"], exons))
        for record in records:
            self.cursor.execute('''
                INSERT IGNORE INTO gtf (id, gene_id, sample, strand, fpkm, exon_positions)
                VALUES (%s, %s, %s, %s, %s, %s)
            ''', record)
        self.conn.commit()

    def process_sample(self, genotype, treatment):
        sample = f"{genotype}_{treatment}"
        gtf_file = os.path.join(self.config.working_directory, f"{sample}.gtf")
        bam_list = ' '.join([f"{self.config.bam_directory}/{sample}_{i+1}_sorted.bam" for i in range(3)])

        if not os.path.exists(gtf_file):
            print(f"[GTF] Generating for {sample}")
            os.system(f"stringtie {bam_list} -G {self.config.canonical_gtf_file} -o {gtf_file}")
        else:
            print(f"[GTF] Skipping {sample}, already exists")

        csv_file = os.path.join(self.config.working_directory, f"{sample}.csv")
        self.gene_ids_from_csv(gtf_file, csv_file)

        modified_gtf_file = os.path.join(self.config.working_directory, f"{sample}_modified.gtf")
        self.modify_gtf(gtf_file, csv_file, modified_gtf_file)

        self.insert_gtf_records(sample, modified_gtf_file)

    def close(self):
        self.cursor.close()
        self.conn.close()


class ASProcessor:
    def __init__(self, config: Config):
        self.config = config
        self.conn = mysql.connector.connect(**config.db_config)
        self.cursor = self.conn.cursor()

    def parse_gtf_exons(self, gtf_file):
        transcripts = defaultdict(lambda: {"exons": []})
        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                start, end = int(parts[3]), int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                attr_dict = {}
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if attr:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.replace('"', '')
                transcript_id = attr_dict.get("transcript_id")
                gene_id = attr_dict.get("gene_id")
                if transcript_id:
                    transcripts[transcript_id]["strand"] = strand
                    transcripts[transcript_id]["gene_id"] = gene_id
                    if feature_type == "exon":
                        transcripts[transcript_id]["exons"].append((start, end))
        return transcripts

    def extract_transcript_sequence(self, sample, transcripts):
        fasta_file = os.path.join(self.config.working_directory, f"{sample}.fa")
        fasta_index = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        sequences = {}
        for tid, data in transcripts.items():
            gene_id = data["gene_id"]
            chrom = gene_id[2]  # Arabidopsis: 3rd char is chromosome
            exons = sorted(data["exons"], key=lambda x: x[0])
            full_seq = ""
            for start, end in exons:
                exon_seq = str(fasta_index[chrom].seq[start - 1:end])
                full_seq += exon_seq
            if data["strand"] == "-":
                full_seq = str(Seq(full_seq).reverse_complement())
            sequences[tid] = full_seq
        return sequences

    def predict_protein_sequence(self, transcript_seqs):
        return {tid: str(Seq(seq).translate(to_stop=True)) for tid, seq in transcript_seqs.items()}

    def predict_cds_sequence(self, protein_seqs, transcript_seqs):
        return {tid: transcript_seqs[tid][:len(prot) * 3] for tid, prot in protein_seqs.items()}

    def insert_sequences_to_db(self, sample, transcripts, transcript_seqs, protein_seqs, cds_seqs):
        for tid in transcripts:
            uid = f"{sample}_{tid}"
            gene_id = transcripts[tid]['gene_id']
            strand = transcripts[tid]['strand']
            transcript = transcript_seqs.get(tid)
            protein = protein_seqs.get(tid)
            cds = cds_seqs.get(tid)
            if transcript and protein and cds:
                self.cursor.execute('''
                    INSERT IGNORE INTO sequences (id, gene_id, sample, strand, nuc_seq, aa_seq, cds_seq)
                    VALUES (%s, %s, %s, %s, %s, %s, %s)
                ''', (uid, gene_id, sample, strand, transcript, protein, cds))
        self.conn.commit()

    def process_sample(self, genotype, treatment):
        sample = f"{genotype}_{treatment}"
        gtf_file = os.path.join(self.config.working_directory, f"{sample}_modified.gtf")
        transcripts = self.parse_gtf_exons(gtf_file)
        transcript_seqs = self.extract_transcript_sequence(sample, transcripts)
        protein_seqs = self.predict_protein_sequence(transcript_seqs)
        cds_seqs = self.predict_cds_sequence(protein_seqs, transcript_seqs)
        self.insert_sequences_to_db(sample, transcripts, transcript_seqs, protein_seqs, cds_seqs)

    def close(self):
        self.cursor.close()
        self.conn.close()


def main():
    config = Config()

    print("[STEP 1] GTF feldolgozás")
    gtf_proc = GTFProcessor(config)
    for genotype in config.genotypes:
        for treatment in config.treatments:
            gtf_proc.process_sample(genotype, treatment)
    gtf_proc.close()

    print("[STEP 2] Szekvencia becslés")
    as_proc = ASProcessor(config)
    for genotype in config.genotypes:
        for treatment in config.treatments:
            as_proc.process_sample(genotype, treatment)
    as_proc.close()


if __name__ == "__main__":
    main()
