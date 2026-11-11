from multiprocessing import Value
import math
import matplotlib.pyplot as plt


vertices = []

for x in range(2):
    for y in range(2):
        for z in range(2):
            vertices.append([x, y, z])

edges = []

for i in range(len(vertices)):
    for j in range(i + 1, len(vertices)):
        differences = sum(
            vertices[i][k] != vertices[j][k]
            for k in range(3)
        )

        if differences == 1:
            edges.append((i, j))

fig = plt.figure()
ax = fig.add_subplot(projection = "3d")

ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.plot([-5, 5], [0, 0], [0, 0], color='black', linewidth=1.5)  
ax.plot([0, 0], [-5, 5], [0, 0], color='black', linewidth=1.5)
ax.plot([0, 0], [0, 0], [-5, 5], color='black', linewidth=1.5) 

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

for start, end in edges:
    ax.plot(
        [vertices[start][0], vertices[end][0]],
        [vertices[start][1], vertices[end][1]],
        [vertices[start][2], vertices[end][2]],
        color="blue"
    )


plt.show()


def matrix(matrice):
    matrice = matrice.strip().strip("()[]")
    rows = matrice.split(",")
    input_matrix = [[float(x) for x in row.split()] for row in rows]
    size = len(input_matrix)
    if size not in (2, 3):
        raise ValueError("Invalid matrix size")
    if any(len(row) != size for row in input_matrix):
        raise ValueError("Invalid dimensions")
    return input_matrix

def format_matrix(M):
    if isinstance(M, str):
        return M
    cells = [[f"{round(v, 4) + 0:g}" for v in row] for row in M]
    width = max(len(c) for row in cells for c in row)

    lines = []
    for i, row in enumerate(cells):
        body = " ".join(c.rjust(width) for c in row)
        left = "( " if i == 0 else " "
        right = " )" if i == len(cells) - 1 else ""
        lines.append(left + body + right)
    return "\n".join(lines)

def determinant(A):
    size = len(A)
    if size == 1:
        return A[0][0]

    if size == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]

    det = 0

    for col in range(size):
        minor = [
            row[:col] + row[col+1:]
            for row in A[1:]
        ]
        det += ((-1)** col) * A[0][col] * determinant(minor)

    return det

def adjugate(A):
    cofactors = []
    size = len(A)
    for row in range(size):
        cofactor_row = []
        for col in range(size):
            minor = [
                r[:col] + r[col+1:]
                for r in A[:row] + A[row+1:]
            ]
            cofactor = ((-1) ** (row + col)) * determinant(minor)
            cofactor_row.append(cofactor)

        cofactors.append(cofactor_row)

    adjugate_matrix = [
        [cofactors[row][col] for row in range(size)]
        for col in range(size)
    ]
    return adjugate_matrix

def inverse_matrix(A):
    adj = adjugate(A)
    det = determinant(A)

    if det == 0:
        return "Matrix is singular; inverse does not exist"

    inverse_matrice = [
        [value / det for value in row]
        for row in adj
    ]
    return inverse_matrice

operations = {
    "+": lambda a,b: a+b,
    "-": lambda a,b: a-b 
}

def matrix_additive_operations(A, B, operator):
    if len(A) != len(B) or any(len(a) != len(b) for a,b in zip(A, B)):
        raise ValueError("Matrices must have matching dimensions")
    result = []

    for row in range(len(A)):
        new_row = []
        
        for col in range(len(A[row])):
            value = operations[operator](A[row][col], B[row][col])
            new_row.append(value)

        result.append(new_row)
    return result

def matrix_multiplication(A, B):
    if len(A[0]) != len(B):
        raise ValueError("Incompatible matrix dimensions")
    result = []
    
    for row in range(len(A)):
        new_row = []
        
        for col in range(len(B[0])):
            total = 0
            for k in range(len(A[0])):
                total += A[row][k] * B[k][col]

            new_row.append(total)

        result.append(new_row)
    return result

matrice = input("Enter matrix in form (a b, c d): ")
A = matrix(matrice)
print(f"""
Matrix : \n{format_matrix(A)}
Determinant: \n{determinant(A)}
Adjugate: \n{format_matrix(adjugate(A))}
Inverse: \n {format_matrix(inverse_matrix(A))}
""")
matrice_two = input("Enter second matrix: ")
B = matrix(matrice_two)
print(f"""
Addition: \n{format_matrix(matrix_additive_operations(A, B, "+"))}
Subtraction: \n{format_matrix(matrix_additive_operations(A, B, "-"))}
Multiplication: \n{format_matrix(matrix_multiplication(A, B))}
""")

Identity_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

def discriminant(a, b, c):
    return b ** 2 - 4 * a * c
def eigenvalues(A):
    trace = A[0][0] + A[1][1]
    det = determinant(A)

    disc = discriminant(1, -trace, det)
    if disc < 0:
        return []
    if disc == 0:
        return [trace/2]
    
    root = math.sqrt(disc)
    return [(trace + root) / 2, (trace - root) / 2]

def eigenvector(A, lam):
    p, q = A[0][0] - lam, A[0][1]
    if abs(p) < 1e-9 and abs(q) < 1e-9:
        p, q = A[1][0], A[1][1] - lam
    if abs(p) < 1e-9 and abs(q) < 1e-9:
        return None
    v = [-q, p]
    first = v[0] if abs (v[0]) > 1e-9 else v[1]
    return [x / first + 0 for x in v]

def eigen_report(A):
    if len(A) != 2:
        return "Eigenvalues are only supported for 2x2 matrices"
    values = eigenvalues(A)
    if not values:
        return "No real eigenvalues"

    blocks = []
    for lam in values:
        v = eigenvector(A, lam)
        if v is None:
            blocks.append(f"Eigen value {lam:g}: any vector works")
            continue

        column = [[v[0]], [v[1]]]

        Av = matrix_multiplication(A, column)
        ok = all(abs(Av[i][0] - lam * v[i]) < 1e-6 for i in range(2))

        blocks.append(
            f"Eigenvalue: {round(lam, 4) + 0:g}\n"
            f"Eigenvector: {'✓' if ok else '✗ check failed'}\n"
            f"{format_matrix(column)}"
        )
    return "\n\n".join(blocks)


print(eigen_report(A))

plt.figure()
plt.axhline(0)
plt.axvline(0)

plt.xlim(-5, 5)
plt.ylim(-5, 5)

plt.grid()
plt.quiver(
    0, 0, 1, 1,
    angles = "xy",
    scale_units = "xy",
    scale = 1
)
plt.quiver(
    0, 0, 3, 3,
    angles = "xy",
    scale_units = "xy",
    scale = 1
)
plt.show()
