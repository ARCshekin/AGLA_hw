TOLERANCE = 10 ** -6

class Matrix:
    def __init__(self, matrix):
        self.n = len(matrix[0])
        self.m = len(matrix)
        self.matrix = matrix

    def subtract_rows(self, r1, r2, k):
        for i in range(self.n):
            self.matrix[r2][i] -= self.matrix[r1][i] * k

    def swap_rows(self, r1, r2):
        self.matrix[r1], self.matrix[r2] = self.matrix[r2], self.matrix[r1]

    def gaussian_elimination(self, answers):
        augmented = [row.copy() for row in self.matrix]
        for i in range(self.m):
            augmented[i].append(answers[i])
        num_rows = self.m
        num_cols = self.n + 1
        rank = 0

        for col in range(num_cols - 1):  # Skip augmented column
            pivot_row = -1
            max_val = 0
            for r in range(rank, num_rows):
                current_val = abs(augmented[r][col])
                if current_val > max_val:
                    max_val = current_val
                    pivot_row = r

            if pivot_row == -1:
                continue

            augmented[rank], augmented[pivot_row] = augmented[pivot_row], augmented[rank]

            pivot = augmented[rank][col]
            if abs(pivot) < TOLERANCE:
                continue

            for c in range(col, num_cols):
                augmented[rank][c] /= pivot

            for r in range(rank + 1, num_rows):
                factor = augmented[r][col]
                for c in range(col, num_cols):
                    augmented[r][c] -= factor * augmented[rank][c]

            rank += 1

        for r in range(rank, num_rows):
            if all(abs(x) < TOLERANCE for x in augmented[r][:-1]) and abs(augmented[r][-1]) > TOLERANCE:
                raise "No solution"

        if rank < self.n:
            raise "Infinite solutions"

        solutions = [0.0] * self.n
        for r in reversed(range(rank)):
            pivot_col = next((c for c in range(self.n) if abs(augmented[r][c]) > TOLERANCE), None)
            if pivot_col is None:
                continue

            rhs = augmented[r][-1]
            for c in range(pivot_col + 1, self.n):
                rhs -= augmented[r][c] * solutions[c]
            solutions[pivot_col] = rhs / augmented[r][pivot_col]

        solutions = [round(s, 3) for s in solutions]
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
