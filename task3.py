TOLERANCE = 10 ** -6


class Matrix:
    def __init__(self, matrix=None, m=0, n=0):
        if matrix is not None:
            self.m = len(matrix)
            self.n = len(matrix[0]) if self.m > 0 else 0
            self.matrix = [[elem for elem in row] for row in matrix]
        else:
            self.m = m
            self.n = n
            self.matrix = [[0] * n for _ in range(m)]

    def copy(self):
        return Matrix([[el for el in row] for row in self.matrix])

    def transposed(self):
        return Matrix([[self.matrix[row][col] for row in range(self.m)]
                       for col in range(self.n)])

    def rref(self):
        rref = self.copy()
        rank = 0

        for col in range(rref.n):
            pivot = -1
            for row in range(rank, rref.m):
                if abs(rref.matrix[row][col]) > TOLERANCE:
                    pivot = row
                    break

            if pivot == -1:
                continue

            rref.swap_rows(rank, pivot)
            rref.multiply_row(rank, 1 / rref.matrix[rank][col])

            for row in range(rref.m):
                if row != rank and abs(rref.matrix[row][col]) > TOLERANCE:
                    factor = rref.matrix[row][col]
                    rref.subtract_rows(rank, row, factor)

            rank += 1

        return rref, rank

    def column_space(self):
        rref, rank = self.rref()
        space = []
        pivot_cols = set()

        for row in range(rank):
            for col in range(rref.n):
                if abs(rref.matrix[row][col] - 1) < TOLERANCE:
                    pivot_cols.add(col)
                    break

        for col in pivot_cols:
            space.append([self.matrix[row][col] for row in range(self.m)])

        return space

    def row_space(self):
        rref, rank = self.rref()
        return [row for row in rref.matrix[:rank]]

    def null_space(self):
        rref, rank = self.rref()
        space = []
        pivot_cols = []

        for row in range(rank):
            for col in range(rref.n):
                if abs(rref.matrix[row][col] - 1) < TOLERANCE:
                    pivot_cols.append(col)
                    break

        free_cols = sorted(set(range(rref.n)) - set(pivot_cols))

        for free in free_cols:
            vec = [0] * rref.n
            vec[free] = 1

            for row in reversed(range(rank)):
                pivot_col = pivot_cols[row]
                vec[pivot_col] = -sum(rref.matrix[row][col] * vec[col]
                                      for col in range(pivot_col + 1, rref.n))

            space.append(vec)

        return space

    def left_null_space(self):
        return self.transposed().null_space()

    def subtract_rows(self, r1, r2, k):
        for i in range(self.n):
            self.matrix[r2][i] -= self.matrix[r1][i] * k

    def swap_rows(self, r1, r2):
        self.matrix[r1], self.matrix[r2] = self.matrix[r2], self.matrix[r1]

    def multiply_row(self, row, scalar):
        for i in range(self.n):
            self.matrix[row][i] *= scalar

    def __mul__(self, other):
        if self.n != other.m:
            raise "Matrix dimensions don't match"

        result = [[0] * other.n for _ in range(self.m)]
        for i in range(self.m):
            for j in range(other.n):
                result[i][j] = sum(self.matrix[i][k] * other.matrix[k][j]
                                   for k in range(self.n))
        return Matrix(result)

    def output(self):
        for row in self.matrix:
            print(*[round(el, 2) for el in row])


A = Matrix([
    [1, 0, -1],
    [2, 1, 1],
    [3, 1, 0]
])

"""

A = Matrix([
    [1, 2, 3, 4],
    [2, 4, 6, 8],
    [3, 6, 9, 12]
])

"""

A.output()

print("Rank:", A.rref()[1])

print("Column Space:")
for col in A.column_space():
    print(*col)

print("Row Space:")
for row in A.row_space():
    print(*row)

print("Null Space:")
null_space = A.null_space()
if null_space:
    for vec in null_space: print(*vec)
else:
    print("Null space is trivial (only zero vector)")

print("Left Null Space:")
left_null_space = A.left_null_space()
if left_null_space:
    for vec in left_null_space: print(*vec)
else: print("Left null space is trivial (only zero vector)")
