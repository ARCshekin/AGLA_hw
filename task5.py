from math import sqrt

TOLERANCE = 10 ** -6

def gram_schmidt(vectors, normalize=True):
    if not vectors:
        return []

    basis = [[float(x) for x in v] for v in vectors]
    orthogonal_basis = []

    for i in range(len(basis)):
        v = basis[i].copy()

        for j in range(i):
            projection = project(v, orthogonal_basis[j])
            v = subtract_vectors(v, projection)

        if all(abs(x) < TOLERANCE for x in v):
            continue

        if normalize:
            norm = vector_norm(v)
            if norm > TOLERANCE:
                v = [x / norm for x in v]
            else:
                continue

        orthogonal_basis.append(v)

    return orthogonal_basis


def dot_product(v1, v2):
    return sum(x * y for x, y in zip(v1, v2))


def vector_norm(v):
    return sqrt(dot_product(v, v))


def project(v, u):
    scale = dot_product(v, u) / dot_product(u, u)
    return [scale * x for x in u]


def subtract_vectors(v1, v2):
    return [x - y for x, y in zip(v1, v2)]


def output(vectors):
    for i, v in enumerate(vectors, 1):
        print(f"v_{i}: {[round(x, 4) for x in v]}")


# Example with 4 vectors
if __name__ == "__main__":
    vectors = [
        [1, 1, 1, 1],
        [1, 2, 3, 4],
        [2, 0, 1, 1],
        [0, 1, 0, 1]
    ]

    print("Original vectors:")
    output(vectors)

    orthogonal_basis = gram_schmidt(vectors, normalize=False)
    print("\nOrthogonal basis:")
    output(orthogonal_basis)

    print("\nVerification:")
    for i in range(len(orthogonal_basis)):
        for j in range(i + 1, len(orthogonal_basis)):
            dp = dot_product(orthogonal_basis[i], orthogonal_basis[j])
            print(f"v_{i + 1}·v_{j + 1} = {dp:.2e}")
