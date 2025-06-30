import mysql.connector

def col_NT():
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    a = 0

    cursor.execute("""
    SELECT g1.id, g1.gene_id, g1.fpkm
    FROM gtf g1
    JOIN (
        SELECT gene_id, MAX(fpkm) AS max_fpkm
        FROM gtf
        WHERE sample = 'Col-0_NT'
        GROUP BY gene_id
    ) g2
    ON g1.gene_id = g2.gene_id AND g1.fpkm = g2.max_fpkm
    WHERE g1.sample = 'Col-0_NT'
    """)
    rows = cursor.fetchall()

    for row in rows:
        id, gene_id, fpkm = row
        a += 1
        cursor.execute("""
            UPDATE gtf SET canonical = TRUE
            WHERE id = %s AND sample = 'Col-0_NT'""", (id,))

    conn.commit()
    cursor.close()
    conn.close()

    print(a)

def col_1h():

    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    a = 0

    cursor.execute("""
        SELECT DISTINCT g1.gene_id
        FROM gtf g1
        WHERE g1.sample = 'Col-0_1h'
          AND g1.gene_id NOT IN (
              SELECT DISTINCT gene_id
              FROM gtf
              WHERE sample = 'Col-0_NT'
                AND canonical = TRUE
          )
    """)
    exclusive_gene_ids = [row[0] for row in cursor.fetchall()]

    # print(f"Genes only in Col-0_1d: {len(exclusive_gene_ids)} found.")

    for gene_id in exclusive_gene_ids:
        if not 'STRG' in gene_id:
            cursor.execute("""
                SELECT id, fpkm FROM gtf
                WHERE sample = 'Col-0_1h' AND gene_id = %s
                ORDER BY fpkm DESC LIMIT 1""", (gene_id,))
            result = cursor.fetchone()
            if result:
                id, fpkm = result
                a += 1
                # print(f'Canonical transcript for {gene_id} (only in Col-0_1d): {id} (FPKM={fpkm})')
                cursor.execute("""
                    UPDATE gtf SET canonical = TRUE
                    WHERE id = %s AND sample = 'Col-0_1h'""", (id,))

    conn.commit()
    cursor.close()
    conn.close()

    print(a)



# def col_1d():

#     conn = mysql.connector.connect(
#         host='localhost',
#         user='root',
#         password='',
#         database='as_db'
#     )
#     cursor = conn.cursor()

#     a = 0

#     cursor.execute("""
#         SELECT DISTINCT g1.gene_id
#         FROM gtf g1
#         WHERE g1.sample = 'Col-0_1d'
#           AND g1.gene_id NOT IN (
#               SELECT DISTINCT gene_id
#               FROM gtf
#               WHERE sample = 'Col-0_NT'
#                 AND canonical = TRUE
#           )
#     """)
#     exclusive_gene_ids = [row[0] for row in cursor.fetchall()]

#     # print(f"Genes only in Col-0_1d: {len(exclusive_gene_ids)} found.")

#     for gene_id in exclusive_gene_ids:
#         if not 'STRG' in gene_id:
#             cursor.execute("""
#                 SELECT id, fpkm FROM gtf
#                 WHERE sample = 'Col-0_1d' AND gene_id = %s
#                 ORDER BY fpkm DESC LIMIT 1""", (gene_id,))
#             result = cursor.fetchone()
#             if result:
#                 id, fpkm = result
#                 a += 1
#                 # print(f'Canonical transcript for {gene_id} (only in Col-0_1d): {id} (FPKM={fpkm})')
#                 cursor.execute("""
#                     UPDATE gtf SET canonical = TRUE
#                     WHERE id = %s AND sample = 'Col-0_1d'""", (id,))

#     conn.commit()
#     cursor.close()
#     conn.close()

#     print(a)


def col_1d():
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    cursor.execute("""
        SELECT DISTINCT g1.gene_id
        FROM gtf g1
        WHERE sample = 'Col-0_1d'
        AND g1.gene_id NOT IN (
            SELECT DISTINCT gene_id
            FROM gtf
            WHERE (sample = 'Col-0_NT' OR sample = 'Col-0_1h')
            AND canonical = TRUE
        )
    """)

    exclusive_gene_ids = [row[0] for row in cursor.fetchall()]

    a = 0

    for gene_id in exclusive_gene_ids:
        if not 'STRG' in gene_id:
            # print(gene_id)
            cursor.execute("""
                SELECT id, fpkm FROM gtf
                WHERE sample = 'Col-0_1d' AND gene_id = %s
                ORDER BY fpkm DESC LIMIT 1""", (gene_id,))
            result = cursor.fetchone()
            if result:
                a += 1
                id, fpkm = result
                # print(f'Canonical transcript for {gene_id} (only in Col-0_1d): {id} (FPKM={fpkm})')
                cursor.execute("""
                    UPDATE gtf SET canonical = TRUE
                    WHERE id = %s AND sample = 'Col-0_1d'""", (id,))
    conn.commit()
    cursor.close()
    conn.close()

    print(a)

def rest():
    conn = mysql.connector.connect(
        host='localhost',
        user='root',
        password='',
        database='as_db'
    )
    cursor = conn.cursor()

    cursor.execute("""
        SELECT DISTINCT g1.gene_id FROM gtf g1 WHERE g1.gene_id NOT IN (SELECT DISTINCT gene_id FROM gtf WHERE canonical = True) AND gene_id NOT LIKE 'STRG%';
    """)
    exclusive_gene_ids = [row[0] for row in cursor.fetchall()]

    # print(len(exclusive_gene_ids))

    # a = 0

    for gene_id in exclusive_gene_ids:
        cursor.execute("""
            SELECT COUNT(*) FROM gtf where gene_id = %s""", (gene_id,))
        
        myresult = cursor.fetchall()
        if myresult[-1][-1] == 1:
            print(gene_id)
            cursor.execute("""
                UPDATE gtf SET canonical = TRUE WHERE gene_id = %s""", (gene_id,))
    conn.commit()
    cursor.close()
    conn.close()


# col_NT()

# col_1h()

col_1d()

# upf3_NT()

# rest()