import sympy
from symcirc.analysis import AnalyseCircuit
from symcirc import utils
import time

def cramer_ddd_solve(matrix):
    """
    Determines system solution using Cramer's rule and determinant decision diagram (DDD)
    matrix: sympy.Matrix containing a system equation [A | b].
    Returns list of solutions [x0, x1, ...] as sympy expressions.
    """
    # build coefficient matrix A and rhs vector b
    if matrix.cols < 2:
        raise ValueError("Matrix must be augmented [A|b] with at least 2 columns.")

    # Solve determinant of A (main Cramer determinant)
    mat_a = matrix.copy()
    mat_a.col_del(-1)   # remove last column (b) to get A
    det_a = determinant_by_ddd(mat_a)
    if det_a == 0:
        raise ValueError("Coefficient matrix is singular (determinant == 0).")

    # Solve other Cramer determinants
    results = []
    for col in range(mat_a.cols):
        mat_ax = mat_a.copy()
        for r in range(mat_a.rows):
            mat_ax[r, col] = matrix[r, -1]
        det_ax = determinant_by_ddd(mat_ax)
        results.append(sympy.cancel(det_ax/det_a))
    return results


def determinant_by_ddd(mat):
    """
    Compute determinant of square sympy.Matrix `mat` using recursive Laplace expansion
    along the first row, with memoization. Works with symbolic entries.
    """
    if mat.rows != mat.cols:
        raise ValueError("Determinant requires a square matrix.")
    memo = {}
    return _det_recursive(mat, memo)


def _matrix_key(mat):
    """
    Create an immutable key for memoization from a sympy Matrix.
    Basicaly 'unravel' the matrix into a tuple of tuples, representing
    its structure in an immutable form.
    """
    return tuple(tuple(mat[i, j] for j in range(mat.cols)) for i in range(mat.rows))


def _det_recursive(mat, memo):
    n = mat.rows
    key = _matrix_key(mat)

    # Check memo for cached results
    if key in memo:
        return memo[key]

    # 1-size matrix
    if n == 1:
        val = mat[0]
        memo[key] = val
        return val

    # 2-size matrix
    '''if n == 2:
        val = sympy.simplify(mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0])
        memo[key] = val
        return val'''

    total = 0
    # expand along row 0: determinant = sum_j (-1)^(0+j) * a[0,j] * det(minor(0,j))
    for j in range(n):
        a = mat[0, j]
        if a == 0:
            # skip zero entries
            continue
        # build minor matrix by deleting row 0 and column j
        minor = mat.copy()
        minor.row_del(0)
        minor.col_del(j)
        cofactor_sign = (-1) ** j  # since row index is 0 => (-1)^(0+j) = (-1)^j
        subdet = _det_recursive(minor, memo)
        term = cofactor_sign * a * subdet
        total += term

    total = total
    memo[key] = total
    return total


if __name__ == '__main__':
    netlist = utils.load_file("..\\..\\tests\\netlists\\AC5.txt")
    t1 = time.time()
    circuit = AnalyseCircuit(netlist, "TF", symbolic=True, precision=6, method="two_graph_node", sympy_ilt=True)
    print(circuit.solved_dict)
    M = circuit.eqn_matrix
    print(f"Gauss solve time: {time.time()-t1}")
    t2 = time.time()
    sol = cramer_ddd_solve(M)
    print(f"Cramer+DDD solve time: {time.time()-t2}")
    print(sol)
