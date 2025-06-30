import pandas as pd
import csv
import mysql.connector
import regex as re
import os 
from collections import defaultdict


def gtf_to_csv(canonical_gtf_file):

    gene_ids = []
    starts = []
    ends = []
    strands = []

    canonical_gtf_file = open(canonical_gtf_file, 'r')

    for line in canonical_gtf_file.readlines():
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


    my_dict = {'gene': gene_ids, 'start_position': starts, 'end_position': ends, 'strand': strands}

    df = pd.DataFrame(my_dict)
    filename = 'from_canonical_gtf.csv'
    df.to_csv(filename)
    return filename

def load_gene_intervals(csv_path):
    genes = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            genes.append({
                'gene_id': row['gene'],
                'chromosome': row['gene'][2],
                'start': int(row['start_position']),
                'end': int(row['end_position']),
                'strand': row['strand']
            })
    return genes


def find_matching_gene(gtf_start, gtf_end, gtf_strand, chromosome, genes):
    gtf_start = int(gtf_start)
    gtf_end = int(gtf_end)
    gtf_half = (gtf_start + gtf_end) // 2

    for gene in genes:
        if gene['strand'] != gtf_strand:
            continue

        if gene['chromosome'] != chromosome:
            continue

        gene_start = min(gene['start'], gene['end'])
        gene_end = max(gene['start'], gene['end'])

        if gene_start == gtf_start and gene_end == gtf_end:
            return gene['gene_id']
        
    for gene in genes:
        if gene['strand'] != gtf_strand:
            continue

        if gene['chromosome'] != chromosome:
            continue

        gene_start = min(gene['start'], gene['end'])
        gene_end = max(gene['start'], gene['end'])
        if not (gtf_end < gene_start or gtf_start > gene_end):
            return gene['gene_id']


    return None


def gene_ids_from_csv(gtf_file, genes, csv_file):

    gene_ids = []
    starts = []
    ends = []
    strg_ids = []
    strg_transcripts = []

    gtf_file = open(gtf_file, 'r')

    for line in gtf_file.readlines():
        
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
                fpkm = attr_dict.get("FPKM")
                ref_gene_id = attr_dict.get("ref_gene_id")
                strg_id = attr_dict.get('gene_id')
                strg_transcript_id = attr_dict.get('transcript_id')
                strand = parts[6]

                if strand != '.':

                    matched_gene_id = find_matching_gene(start_position, end_position, strand, chromosome, genes)

                    if matched_gene_id:
                        gene_ids.append(matched_gene_id)
                        strg_ids.append(strg_id)
                        strg_transcripts.append(strg_transcript_id)

    my_dict = {'gene': gene_ids, 'strg_gene_id': strg_ids, 'strg_transcript_id': strg_transcripts}
    df = pd.DataFrame(my_dict)
    df.to_csv(csv_file)


def modify_gtf(gtf_file, csv_file, modified_gtf_file):

    gtf_file = open(gtf_file,'r')

    df = pd.read_csv(csv_file)

    mapping = dict(zip(df['strg_transcript_id'], df['gene']))

    new_lines = []

    for line in gtf_file.readlines():
        if not '#' in line:
            if 'gene_id' in line:
                for strg_transcript_id, gene in mapping.items():
                    if f'transcript_id "{strg_transcript_id}"' in line:
                        strg_gene = strg_transcript_id.rsplit('.', 1)[0]
                        line = line.replace(f'gene_id "{strg_gene}"', f'gene_id "{gene}"')
                        break
        new_lines.append(line)

    with open(modified_gtf_file, "w") as f:
        f.writelines(new_lines)

def parse_gtf_exons(gtf_file):
    transcripts = defaultdict(list)
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                print(parts)
                break

            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
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

            if feature_type == "transcript":
                transcripts[transcript_id] = {
                    "gene_id": gene_id,
                    "strand": strand,
                    "fpkm": float(fpkm) if fpkm else None,
                    "exons": []
                }

            elif feature_type == "exon" and transcript_id:
                transcripts[transcript_id]["exons"].append((start, end))

    return transcripts

def insert_gtf_records(sample, gtf_file):
    gtf_records = []

    transcripts = parse_gtf_exons(gtf_file)
    for transcript_id, data in transcripts.items():
        exon_positions = ";".join([f"{start}-{end}" for start, end in data["exons"]])

        gene_id = data["gene_id"]
        uid = f"{sample}_{transcript_id}"

        gtf_records.append((
            uid,
            gene_id,
            sample,
            data["strand"],
            data["fpkm"],
            exon_positions
        ))

    return gtf_records


def main():
    #mysql connection
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    #define paths
    canonical_gtf_file = '/Users/galno/tmp_for_local_run/Arabidopsis_thaliana.TAIR10.60.gtf'
    bam_directory = '/Volumes/bam'
    working_directory = '/Users/galno/as_2'

    #canonical gtf file to csv and to the genes variable
    canonical_csv_filename = f'{working_directory}/from_canonical_gtf.csv'
    if not os.path.exists(canonical_csv_filename):
        gtf_to_csv(canonical_gtf_file)
    genes = load_gene_intervals(canonical_csv_filename)

    genotypes = ['Col-0', 'tfiis', 'upf1', 'upf3', 'tfiis-upf1', 'tfiis-upf3']
    treatments = ['NT', '1h', '1d']

    # genotypes = ['Col-0']
    # treatments = ['NT']

    for genotype in genotypes:
        for treatment in treatments:
            #create gtf
            sample = f'{genotype}_{treatment}'
            
            gtf_file = f'{working_directory}/{genotype}_{treatment}.gtf'
            bam_list = ' '.join([f'{bam_directory}/{genotype}_{treatment}_{n+1}_sorted.bam' for n in range(3)])
            if not os.path.exists(gtf_file):
                print(f'generate gtf for {sample}')
                os.system(f'stringtie {bam_list} -G {canonical_gtf_file} -o {gtf_file}')
            else:
                print(f'skipping {sample}')

            #create csv
            csv_file = f'{working_directory}/{genotype}_{treatment}.csv'
            gene_ids_from_csv(gtf_file, genes, csv_file)

            #modify gtf
            print(f'modify gtf file: {gtf_file}')
            modified_gtf_file = f'{working_directory}/{genotype}_{treatment}_modified.gtf'
            modify_gtf(gtf_file, csv_file, modified_gtf_file)

            #insert to mysq
            print('insert into mysql...')
            gtf_records = insert_gtf_records(sample, modified_gtf_file)
            for record in gtf_records:
                insert_query = '''
                INSERT IGNORE INTO gtf (id, gene_id, sample, strand, fpkm, exon_positions)
                VALUES (%s, %s, %s, %s, %s, %s)
                '''
                cursor.execute(insert_query, record)

            conn.commit()

            print(f'{sample} is done')
    
    cursor.close()
    conn.close()



if __name__ == "__main__":
    main()




