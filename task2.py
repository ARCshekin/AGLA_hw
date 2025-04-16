TOLERANCE = 10 ** -6

class Matrix:
    def __init__(self, matrix):
        self.n = len(matrix[0])
        self.m = len(matrix)
        self.matrix = matrix

    def subtract_rows(self, r1, r2, k):
        # R2 = R2 - k*R1
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

    def gauss_jordan_inverse(self):
        if self.n != self.m:
            raise "Matrix is not square"

        augmented = [row.copy() + [1 if i == j else 0 for i in range(self.n)]
                     for j, row in enumerate(self.matrix)]

        for col in range(self.n):
            max_row = max(range(col, self.n), key=lambda r: abs(augmented[r][col]))
            if col != max_row:
                augmented[col], augmented[max_row] = augmented[max_row], augmented[col]

            if abs(augmented[col][col]) < TOLERANCE:
                raise "Matrix is singular"

            # Normalize current row
            pivot = augmented[col][col]
            for j in range(col, 2 * self.n):
                augmented[col][j] /= pivot

            for row in range(self.n):
                if row != col and augmented[row][col] != 0:
                    factor = augmented[row][col]
                    for j in range(col, 2 * self.n):
                        augmented[row][j] -= augmented[col][j] * factor

        inverse_matrix = [row[self.n:] for row in augmented]
        return Matrix(inverse_matrix)

    def output(self):
        for row in self.matrix:
            print(*[round(el, 2) for el in row])

# Create a matrix
m = Matrix([
    [4, 7],
    [2, 6]
])

# Calculate its inverse
inv = m.gauss_jordan_inverse()

inv.output()

# Check if A * A-1 = I
verif = m * inv

verif.output()
