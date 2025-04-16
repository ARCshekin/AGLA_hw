class Matrix:
    def __init__(self, matrix):
        self.n = len(matrix[0])
        self.m = len(matrix)
        self.matrix = matrix

    def subtract_rows(self, r1, r2, k):
        # R2 = R2 - k*R1
        for i in range(self.n):
            self.matrix[r2][i] -= self.matrix[r1][i] * k

    def gaussian_elimination(self, answers):
        # Augment matrix with answers
        for i in range(self.m):
            self.matrix[i].append(answers[i])
        self.n += 1

        # Forward pass: create upper triangular matrix
        for col in range(self.n - 1):
            max_row = max(range(col, self.m), key=lambda r: abs(self.matrix[r][col]))
            self.matrix[col], self.matrix[max_row] = self.matrix[max_row], self.matrix[col]

            if abs(self.matrix[col][col]) == 0:
                raise "Non-singular matrix"

            for row in range(col + 1, self.m):
                self.subtract_rows(col, row, self.matrix[row][col] / self.matrix[col][col])

        # Back substitution: solve from bottom up
        solutions = [0] * (self.n - 1)
        for row in reversed(range(self.m)):
            if all(abs(x) < 1e-10 for x in self.matrix[row][:-1]):
                continue
            solutions[row] = (self.matrix[row][-1] -
                              sum(self.matrix[row][col] * solutions[col]
                                  for col in range(row + 1, self.n - 1))) / self.matrix[row][row]

        for i, el in enumerate(solutions):
            solutions[i] = round(el, 3)
        return solutions
'''
---- Example system: ----
 2x + y - z = 8
 -3x - y + 2z = -11
 -2x + y + 2z = -3
 
 Solution: x=2, y=3, z=-1
'''
coefficients = [
    [2, 1, -1],
    [-3, -1, 2],
    [-2, 1, 2]
]
answers = [8, -11, -3]

m = Matrix(coefficients)
solutions = m.gaussian_elimination(answers)
print(solutions)
