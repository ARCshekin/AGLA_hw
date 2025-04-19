def orthogonal_complement(matrix):
    rref_matrix, pivot_columns = rref(matrix)

    rank = len(pivot_columns)
    n = len(matrix[0])

    basis = []
    for j in range(n):
        if j not in pivot_columns:
            vec = [0] * n
            vec[j] = 1
            for i in range(rank): vec[pivot_columns[i]] = -rref_matrix[i][j]
            basis.append(vec)

    return basis


def rref(matrix):
    m = len(matrix)
    n = len(matrix[0]) if m > 0 else 0

    rref = [row.copy() for row in matrix]
    pivot_columns = []

    for r in range(m):
        pivot = r
        while pivot < m and rref[pivot][r] == 0: pivot += 1

        if pivot >= m:
            continue

        if pivot != r:
            rref[r], rref[pivot] = rref[pivot], rref[r]

        pivot_columns.append(r)

        pivot_val = rref[r][r]
        if pivot_val != 0:
            for c in range(n):
                rref[r][c] /= pivot_val

        for i in range(m):
            if i != r and rref[i][r] != 0:
                factor = rref[i][r]
                for c in range(n):
                    rref[i][c] -= factor * rref[r][c]

    return rref, pivot_columns


if __name__ == "__main__":
    A = [
        [1, 2, 3],
        [4, 5, 6]
    ]

    print("Original matrix:")
    for row in A:
        print(row)

    orth_comp = orthogonal_complement(A)
    print("\nOrthogonal complement basis vectors:")
    for vec in orth_comp:
        print(vec)

    print("\nVerification (dot products should be zero):")
    for row in A:
        for vec in orth_comp:
            dot = sum(row[i] * vec[i] for i in range(len(row)))
            print(f"Dot product of {row} and {vec}: {dot}")
