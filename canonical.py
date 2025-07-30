import mysql.connector
from concurrent.futures import ThreadPoolExecutor, as_completed


DB_CONFIG = {
    'host': 'localhost',
    'user': 'root',
    'password': '',
    'database': 'as_db'
}

N_CORES = 10


# 🔄 Reset canonical flag
def reset_canonical_flags():
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    cursor.execute("UPDATE filtered SET canonical = FALSE")
    conn.commit()
    cursor.close()
    conn.close()
    print("Canonical flags reset to FALSE.")


# 🔍 Helper — Get gene IDs missing canonical
def get_missing_genes(collected_genes):
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()
    cursor.execute("SELECT DISTINCT gene_id FROM filtered")
    all_genes = {row[0] for row in cursor.fetchall()}
    missing = all_genes - collected_genes
    cursor.close()
    conn.close()
    return missing


# 🧠 Canonical selection for a given sample set
def set_canonical_for_sample(sample_list, exclude_genes):
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()

    placeholders = ','.join(['%s'] * len(sample_list))

    cursor.execute(f"""
        SELECT DISTINCT gene_id
        FROM filtered
        WHERE sample IN ({placeholders})
    """, tuple(sample_list))
    gene_list = {row[0] for row in cursor.fetchall()}

    target_genes = gene_list - exclude_genes

    def process_gene(gene_id):
        local_conn = mysql.connector.connect(**DB_CONFIG)
        local_cursor = local_conn.cursor()

        local_cursor.execute(f"""
            SELECT id, fpkm
            FROM filtered
            WHERE gene_id = %s AND sample IN ({placeholders})
        """, (gene_id, *sample_list))
        rows = local_cursor.fetchall()

        if not rows:
            local_cursor.close()
            local_conn.close()
            return None

        # Select row with maximum FPKM
        best_id = max(rows, key=lambda x: x[1])[0]

        local_cursor.execute("""
            UPDATE filtered SET canonical = TRUE WHERE id = %s
        """, (best_id,))
        local_conn.commit()

        local_cursor.close()
        local_conn.close()
        return gene_id

    # Run parallel
    with ThreadPoolExecutor(max_workers=N_CORES) as executor:
        futures = {executor.submit(process_gene, gid): gid for gid in target_genes}
        results = []
        for future in as_completed(futures):
            gene_id = futures[future]
            try:
                result = future.result()
                if result:
                    results.append(result)
            except Exception as e:
                print(f"Error processing gene {gene_id}: {e}")

    cursor.close()
    conn.close()

    return set(results)


# 🔗 Canonical propagate by sequence
def canonical_same():
    conn = mysql.connector.connect(**DB_CONFIG)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT DISTINCT nuc_seq
        FROM filtered
        WHERE canonical = TRUE AND sample IN ('Col-0_NT', 'Col-0_1h', 'Col-0_1d')
    """)
    canonical_seqs = [row[0] for row in cursor.fetchall()]

    if not canonical_seqs:
        print("No canonical sequences found in Col-0.")
        cursor.close()
        conn.close()
        return

    print(f"Found {len(canonical_seqs)} canonical sequences to propagate.")

    def process_seq(nuc_seq):
        local_conn = mysql.connector.connect(**DB_CONFIG)
        local_cursor = local_conn.cursor()

        local_cursor.execute("""
            UPDATE filtered
            SET canonical = TRUE
            WHERE nuc_seq = %s
        """, (nuc_seq,))
        affected = local_cursor.rowcount
        local_conn.commit()

        local_cursor.close()
        local_conn.close()
        return affected

    # Run parallel
    with ThreadPoolExecutor(max_workers=N_CORES) as executor:
        futures = [executor.submit(process_seq, seq) for seq in canonical_seqs]
        total = 0
        for future in as_completed(futures):
            try:
                count = future.result()
                total += count
            except Exception as e:
                print(f"Error updating canonical for a sequence: {e}")

    print(f"Canonical flags propagated to {total} rows.")

    cursor.close()
    conn.close()


# 🏃 Main pipeline
def main():
    print("Resetting canonical flags...")
    reset_canonical_flags()

    collected_genes = set()

    print("Setting canonical for Col-0_NT...")
    genes_NT = set_canonical_for_sample(['Col-0_NT'], set())
    collected_genes |= genes_NT
    print(f"Canonical genes from NT: {len(genes_NT)}")

    print("Setting canonical for Col-0_1h...")
    genes_1h = set_canonical_for_sample(['Col-0_1h'], collected_genes)
    collected_genes |= genes_1h
    print(f"Canonical genes from 1h: {len(genes_1h)}")

    print("Setting canonical for Col-0_1d...")
    genes_1d = set_canonical_for_sample(['Col-0_1d'], collected_genes)
    collected_genes |= genes_1d
    print(f"Canonical genes from 1d: {len(genes_1d)}")

    print("Propagating canonical flags based on sequence identity...")
    canonical_same()

    print("All done.")


if __name__ == "__main__":
    main()
