import mysql.connector


def short():
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    # Canonical transzkriptek betöltése memóriába
    cursor.execute('SELECT gene_id, as_seq FROM filtered WHERE canonical = TRUE')
    canonical_rows = cursor.fetchall()
    canonical_dict = {gene_id: as_seq for gene_id, as_seq in canonical_rows}

    print(f"Canonical transzkriptek betöltve: {len(canonical_dict)}")

    # Nem-canonical sorok betöltése
    cursor.execute('SELECT id, gene_id, as_seq FROM filtered WHERE canonical = FALSE')
    rows = cursor.fetchall()

    print(f"Összes nem-canonical sor: {len(rows)}")

    updates = []

    for count, (id, gene_id, as_seq) in enumerate(rows, start=1):
        if count % 1000 == 0:
            print(f"{count} sor feldolgozva")

        canonical_as_seq = canonical_dict.get(gene_id)

        if as_seq and canonical_as_seq:
            if len(as_seq) < len(canonical_as_seq):
                updates.append((id,))

    print(f"Rövidebb isoformák száma: {len(updates)}")

    # Batch update
    if updates:
        cursor.executemany("UPDATE filtered SET short = TRUE WHERE id = %s", updates)
        conn.commit()
        print("UPDATE kész")

    cursor.close()
    conn.close()


short()
