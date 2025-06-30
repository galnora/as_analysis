import mysql.connector
import pandas as pd


def collect_genes(sheet, treatment):
    df = pd.read_excel('/Users/galno/as/meta.xlsx', sheet_name=sheet)
    filtered_df = df[df['treatment'] == treatment]
    gene_ids = list(set(filtered_df['gene'].tolist()))
    print(f"  Treatment: {treatment} — {len(gene_ids)} genes to process")
    return gene_ids


def main():
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    # genotypes = ['tfiis', 'upf1', 'upf3', 'tfiis-upf1', 'tfiis-upf3']
    genotypes = ['tfiis-upf1', 'tfiis-upf3']
    treatments = ['NT', '1h', '1d']

    result_rows = []

    for genotype in genotypes:
        print(f'Processing genotype: {genotype}')
        sheet = f'tfiis_{genotype}'

        for treatment in treatments:
            gene_ids = collect_genes(sheet, treatment)

            if not gene_ids:
                print(f"    No genes for {genotype} {treatment}, skipping")
                continue

            samples = [f'{genotype}_{treatment}', f'tfiis_{treatment}']

            gene_placeholders = ','.join(['%s'] * len(gene_ids))
            sample_placeholders = ','.join(['%s'] * len(samples))

            query = f"""
                SELECT as_seq
                FROM seq
                WHERE gene_id IN ({gene_placeholders})
                AND sample IN ({sample_placeholders})
            """

            cursor.execute(query, gene_ids + samples)
            rows = cursor.fetchall()

            unique_as_seqs = set(row[0] for row in rows)
            unique_count = len(unique_as_seqs)

            label = f'tfiis_{genotype}_{treatment}'
            print(f"    {label}: {unique_count} unique transcripts")

            result_rows.append({
                'sample': label,
                'unique_transcript_count': unique_count
            })

    result_df = pd.DataFrame(result_rows)
    result_df.to_excel('/Users/galno/as/summary_transcript_counts_by_treatment_tfiis.xlsx', index=False)

    cursor.close()
    conn.close()


main()
