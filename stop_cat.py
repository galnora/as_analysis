import mysql.connector


def main():
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    # 🔥 Lekérjük az összes canonical szekvenciát
    cursor.execute('''
        SELECT gene_id, as_seq 
        FROM filtered 
        WHERE canonical = TRUE
    ''')
    canonical_rows = cursor.fetchall()
    canonical_dict = {gene_id: as_seq.strip('*') for gene_id, as_seq in canonical_rows}

    print(f'{len(canonical_dict)} canonical sequences loaded')

    # 🔥 Lekérjük az összes nem-canonical sort
    cursor.execute('''
        SELECT 
            id, 
            gene_id, 
            as_seq, 
            short 
        FROM 
            filtered
        WHERE 
            canonical = FALSE AND as_seq IS NOT NULL
    ''')

    rows = cursor.fetchall()
    print(f'{len(rows)} rows loaded')

    update_query = "UPDATE filtered SET short_group = %s WHERE id = %s"

    for count, row in enumerate(rows, 1):
        id, gene_id, as_seq, short = row

        canonical_seq = canonical_dict.get(gene_id)

        if canonical_seq is None:
            # 🔍 Plusz check: gene_id szerepel-e bármelyik Col-0 mintában?
            cursor.execute('''
                SELECT COUNT(*) 
                FROM filtered 
                WHERE gene_id = %s AND sample LIKE "Col-0%%"
            ''', (gene_id,))
            col0_count = cursor.fetchone()[0]

            if col0_count > 0:
                print(f'WARNING: No canonical for {gene_id}, but present in {col0_count} Col-0 samples')
            # else:
            #     print(f'No canonical for {gene_id} and not found in any Col-0 sample')

            continue

        as_seq = as_seq.strip('*')

        # Utolsó 10 nukleotid összehasonlítása
        if as_seq[-10:] == canonical_seq[-10:]:
            short_group = 'short-original' if short else 'long-original'
        else:
            short_group = 'short-neo' if short else 'long-neo'

        cursor.execute(update_query, (short_group, id))

        if count % 1000 == 0:
            print(f'{count} rows processed')

    conn.commit()
    cursor.close()
    conn.close()


if __name__ == "__main__":
    main()
