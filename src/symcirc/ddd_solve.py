"""
LED-DDD: Layered Expansion Diagram based Determinant Decision Diagram
with subgraph reuse Cramer's Rule solver

Based on: Shi, G. (2010). "A Simple Implementation of Determinant Decision Diagram."

Cramer's rule subgraph reuse:
    Solving A*x = b via Cramer's rule requires det(A_k) for each k, where
    A_k is A with column k replaced by b. A_k shares (n-1)/n of its
    sub-minors with A. CramerBuilder reuses the hash table from the
    original det(A) build.
    This avoids re-expanding most of the DDD subgraph.

Hashing notes:
    Minors: Each minor of A is uniquely identified by its row and column index sets.
    Bitmask keys: encode the row/column index sets as Python ints (bitmasks), one bit per index.
    Python int arithmetic and hashing are implemented in C and much faster than constructing sorted tuples.
    For n <= 64 each bitmask fits in a word. For larger n Python uses arbitrary-precision ints automatically.

Optimizations:
Precomputed sparsity mask: frozenset of (r,c) pairs replaces repeated A[r,c]!=0 SymPy __eq__ calls
Bitmask minor keys: (row_bits, col_bits) int pair replaces sorted-tuple construction
Shared evaluate() cache: det(A) and all det(A_k) share one memoisation dict - shared subgraphs are evaluated once
Bisect+slice tuple ops: replace sorted(rows).index() and genexpr filtering with bisect_left + slice
"""

from __future__ import annotations

from bisect import bisect_left
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import sympy as sp

# Terminal node constants for DDD
TERMINAL_ONE  = 1
TERMINAL_ZERO = 0

_node_id_counter = 0

def _get_next_node_id():
    global _node_id_counter
    _node_id_counter += 1
    return _node_id_counter

@dataclass
class DDDNode:
    """
    A non-terminal (not 0 or 1) DDD vertex.

    symbol   : the matrix element a[abs_row, abs_col]  (SymPy expr or number)
    abs_row  : absolute row index in the matrix
    abs_col  : absolute column index in the matrix
    sign     : +/-1  -- cofactor sign (-1)**(local_row + local_col) where local_row / local_col are 1-based
               (starts at 1 not 0) positions inside the current minor at expansion time (Shi 2000, Section III)
    one      : 1-edge child. Points to the DDD representing the cofactor (minor).
    zero     : 0-edge child. Points to the next sibling in the expansion chain (the remainder).
    """
    symbol: Any
    abs_row: int
    abs_col: int
    sign: int
    one: Any = None
    zero: Any = None

    _uid: int = field(default_factory=_get_next_node_id)
    _cached: Any = field(default=None)

    @property
    def id(self):
        return self._uid

    def __repr__(self) -> str:
        s = "+" if self.sign >= 0 else "-"
        return f"DDD[{s}a[{self.abs_row},{self.abs_col}]]"


# Two-element tuple representing the minor key by a bitmask
Minor = Tuple[int, int]   # (row_bitmask, col_bitmask)


def _minor_key(rows: tuple, cols: tuple) -> Minor:
    """
    Encodes the row and column into bitmasks to optimize hashing speed and equality checks.
    """
    r_bits = 0
    c_bits = 0
    for r in rows:
        r_bits |= 1 << r
    for c in cols:
        c_bits |= 1 << c
    return r_bits, c_bits


class MinorHashTable:
    """
    Maps {tuple(row_bitmask, col_bitmask): node (DDD subgraph root)}.
    """
    def __init__(self) -> None:
        self._table: Dict[Minor, Any] = {}

    def lookup(self, minor_key: Minor) -> Optional[Any]:
        return self._table.get(minor_key)

    def store(self, minor_key: Minor, node: Any) -> None:
        self._table[minor_key] = node

    def __len__(self) -> int:
        return len(self._table)


def _row_degrees(nonzero: frozenset, rows: tuple, cols: tuple) -> Dict[int, int]:
    """
    Count nonzero entries per row using the precomputed sparsity set.
    """
    t0 = time.time()
    ret = {}
    for r in rows:
        nonzero_count = 0
        for c in cols:
            if (r, c) in nonzero:
                nonzero_count += 1
        ret[r] = nonzero_count
    return ret


def _select_expansion_row(nonzero: frozenset, rows: tuple, cols: tuple, strategy: str) -> int:
    """
    Picks the expansion row for the determinant.
    "min_degree": row with fewest nonzeros in the current minor. Best for sparse circuit matrices.
    "first": first remaining row in original order.
    """

    if strategy == "min_degree":
        deg = _row_degrees(nonzero, rows, cols)
        # primary heuristic is degree of row, secondary is row index, ensures canonicity of the DDD in case of ties
        # (as opposed to the second heuristic being random choice)
        return min(rows, key=lambda r: (deg[r], r))
    else: # default to "first" row strategy
        return rows[0]


def _minor_is_singular(nonzero: frozenset, rows: tuple, cols: tuple) -> bool:
    """
    Singularity check using the precomputed nonzero set.
    Return True if any row or column in the sub-minor is structurally all-zero.
    """
    for r in rows:
        if not any((r, c) in nonzero for c in cols):
            return True

    for c in cols:
        if not any((r, c) in nonzero for r in rows):
            return True
    return False


def _minor_Ak_is_singular(nonzero_A: frozenset, b_vec: list, k: int, rows: tuple, cols: tuple) -> bool:
    """
    Singularity check for A_k (column k replaced by b_vec).
    Modification of the more general _minor_is_singular().

    For column k it additionally checks b_vec[r] != 0.
    """
    for r in rows:
        if all((b_vec[r] == 0 if c == k else (r, c) not in nonzero_A) for c in cols):
            return True
    for c in cols:
        if c == k: # b_vec check for col k
            if all(b_vec[r] == 0 for r in rows):
                return True
        else: # standard logic
            if not any((r, c) in nonzero_A for r in rows):
                return True
    return False


def evaluate(node: Any, cache: Optional[dict] = None) -> Any:
    """
    Evaluate a DDD to its value.
    """
    if cache is None:
        cache = {}

    if node is TERMINAL_ZERO:
        return sp.S.Zero
    if node is TERMINAL_ONE:
        return sp.S.One
    if not isinstance(node, DDDNode):
        return sp.sympify(node)

    nid = node.id
    if nid in cache:
        return cache[nid]

    zero_val = evaluate(node.zero, cache)
    one_val  = evaluate(node.one, cache)

    result = zero_val + node.sign * node.symbol * one_val
    cache[nid] = result
    return result


def evaluate_nested(node: Any, cache: Optional[dict] = None) -> Any:
    """
    Evaluate a DDD to its value in the nested factored form.
    The factored form is exponentially more compact for sparse matrices.
    Example:
        Standard eval: a*(b+c) + d*(e+f) → ab + ac + de + df (expanded)
        Nested eval:   a*(b+c) + d*(e+f) → a*(b+c) + d*(e+f) (factored)
    """
    if cache is None:
        cache = {}

    # Trivial cases
    if node is TERMINAL_ZERO: # Node is zero
        return sp.S.Zero
    if node is TERMINAL_ONE: # Node is one
        return sp.S.One
    if not isinstance(node, DDDNode): # Node is a leaf
        return sp.sympify(node)

    # Check cache for already solved subdeterminants
    nid = node.id
    if nid in cache: # Cache hit
        return cache[nid]

    # Recursively evaluate children
    zero_val = evaluate_nested(node.zero, cache)
    one_val = evaluate_nested(node.one, cache)

    # Build cofactor without expansion using Mul and Add
    # (use of simple '*', '+' causes automatic simplifications)
    node_is_zero = (node.symbol == TERMINAL_ZERO)
    one_edge_is_zero = (one_val == TERMINAL_ZERO)
    node_is_one = (node.symbol == TERMINAL_ONE)
    one_edge_is_one = (one_val == TERMINAL_ONE)

    if node_is_zero or one_edge_is_zero: # to avoid explicit mult by 0 (e.g. 0*a)
        cofactor = TERMINAL_ZERO
    elif node_is_one: # to avoid explicit mult by 1 (e.g. 1*a)
        cofactor = one_val
    elif one_edge_is_one: # to avoid explicit mult by 1 (e.g. a*1)
        cofactor = node.symbol
    else: # general case
        if node.sign == 1:
            cofactor = sp.Mul(node.symbol, one_val, evaluate=False)
        else:
            cofactor = sp.Mul(-1, node.symbol, one_val, evaluate=False)
    # Combine with remainder, avoiding explicit zero add (e.g. a+0)
    if zero_val == 0:
        result = cofactor
    elif cofactor == 0:
        result = zero_val
    else:
        result = sp.Add(zero_val, cofactor, evaluate=False)

    cache[nid] = result
    return result


def ddd_size(root: Any) -> int:
    """
    Count non-terminal DDD vertices reachable from root.
    """
    visited: set = set()

    def _visit(node: Any) -> None:
        if node is TERMINAL_ZERO or node is TERMINAL_ONE:
            return
        if not isinstance(node, DDDNode):
            return
        nid = node.id
        if nid in visited:
            return
        visited.add(nid)
        _visit(node.one)
        _visit(node.zero)

    _visit(root)
    return len(visited)


class DDDBuilder:
    """
    Builder class for constructing LED-DDDs of det(A) and det(A_k) for Cramer's.

    If k is None, it builds the standard determinant DDD.
    If k is an integer, it builds the numerator for Cramer's rule DDD, reusing the expansion order and hash table
    from a reference builder.
    """

    def __init__(self, A: sp.Matrix, strategy: str = "min_degree",
                 b_vec: Optional[List[Any]] = None, k: Optional[int] = None,
                 reference: Optional[DDDBuilder] = None) -> None:
        self.A = A
        self.n = A.shape[0]
        self.strategy = strategy
        self.b_vec = b_vec
        self.k = k

        if reference: # Cramer Ak det building
            self._nonzero_A = reference._nonzero_A
            self._expansion_rows = reference._expansion_rows
            self._hash_shared = reference._hash_local  # Reuse hash from det(A)
        else: # Standard A det building
            self._nonzero_A = self._build_nonzero_set(A)
            self._expansion_rows = {}
            self._hash_shared = None

        self._hash_local = MinorHashTable()
        self.node_count = 0

    @staticmethod
    def _build_nonzero_set(matrix: sp.Matrix) -> frozenset:
        """
        Build a frozenset of (r, c) pairs where A[r, c] != 0.
        Note: float(0) != sympy.S.Zero
        """
        zero = sp.S.Zero
        n = matrix.shape[0]
        return frozenset((r, c) for r in range(n) for c in range(n) if matrix[r, c] != zero)

    def build(self) -> Any:
        """
        Builds the LED-DDD.
        """
        rows = tuple(range(self.n))
        cols = tuple(range(self.n))
        return self._expand(rows, cols)

    def _expand(self, rows: tuple, cols: tuple) -> Any:
        """
        Recursively expand det(A[rows, cols]) into a LED-DDD.
        """
        # 1. Try empty minor check
        if not rows:
            return TERMINAL_ONE

        # 2. Try hash table lookup
        key = _minor_key(rows, cols)
        # 2a. Try local cache (Ak)
        cached = self._hash_local.lookup(key)
        if cached is not None:
            return cached
        # 2b. Try shared cache (A) if it exists
        if self._hash_shared:
            if self.k not in cols:
                shared = self._hash_shared.lookup(key)
                if shared is not None:
                    return shared

        # 3. Try singularity check
        if self.k in cols:
            is_singular = _minor_Ak_is_singular(self._nonzero_A, self.b_vec, self.k, rows, cols)
        else:
            is_singular = _minor_is_singular(self._nonzero_A, rows, cols)
        if is_singular:
            self._hash_local.store(key, TERMINAL_ZERO)
            return TERMINAL_ZERO

        # 4. Try leaf (1x1 minor) check
        if len(rows) == 1:
            r, c = rows[0], cols[0]
            if c == self.k: # For Cramer det(Ak)
                val = self.b_vec[r]
            else: # For standard det(A)
                val = self.A[r, c]
            if val == 0:
                self._hash_local.store(key, TERMINAL_ZERO)
                return TERMINAL_ZERO
            node = DDDNode(val, r, c, 1, TERMINAL_ONE, TERMINAL_ZERO)
            self.node_count += 1
            self._hash_local.store(key, node)
            return node

        # 5. General case: select, record and remove expansion row
        expansion_row = self._expansion_rows.get(key) # For det(Ak)
        if expansion_row is None: # For det(A)
            expansion_row = _select_expansion_row(self._nonzero_A, rows, cols, self.strategy)
            self._expansion_rows[key] = expansion_row

        # bisect_left on already-sorted rows gives position without re-sorting.
        local_row_pos = bisect_left(rows, expansion_row) + 1
        rem_rows = rows[:local_row_pos - 1] + rows[local_row_pos:]

        leader_node = None
        prev_node = None
        for local_col_pos, c in enumerate(cols, start=1):
            if c == self.k: # For det(Ak)
                val = self.b_vec[expansion_row]
            else:   # For det(A)
                val = self.A[expansion_row, c]
            if val == 0:
                continue

            # Compute subdeterminant sign
            if (local_row_pos + local_col_pos) & 1 == 0:
                sign = 1
            else:
                sign = -1

            # New node
            node = DDDNode(val, expansion_row, c, sign)
            self.node_count += 1

            # remove column c from sorted cols using bisect+slice
            c_idx = bisect_left(cols, c)
            rem_cols = cols[:c_idx] + cols[c_idx + 1:]

            # Core recursion
            node.one = self._expand(rem_rows, rem_cols)

            if not leader_node:
                leader_node = node
            else:
                prev_node.zero = node
            prev_node = node

        if not leader_node:
            self._hash_local.store(key, TERMINAL_ZERO)
            return TERMINAL_ZERO

        prev_node.zero = TERMINAL_ZERO
        self._hash_local.store(key, leader_node)
        return leader_node


class DDDSolver:
    """
    Solver for standard determinant and Cramer's Rule subdeterminants (with reference cache sharing).
    """
    def __init__(self, A: sp.Matrix, strategy: str = "min_degree"):
        self.A = A
        self.strategy = strategy
        self.denominator_builder = None
        self._shared_eval_cache = {}

    def solve(self, b: Optional[Sequence[Any]] = None) -> Any:
        # Build Denominator det(A)
        self.denominator_builder = DDDBuilder(self.A, self.strategy)
        root_A = self.denominator_builder.build()
        det_A = evaluate(root_A, self._shared_eval_cache)

        if b is None: # Simple determinant (not Cramer's rule)
            return det_A

        else:
            # Build Numerators det(Ak)
            b_list = [sp.sympify(x) for x in b]
            results = []

            for k in range(self.A.shape[0]):
                num_builder = DDDBuilder(self.A, b_vec=b_list, k=k, reference=self.denominator_builder)
                root_Ak = num_builder.build()
                det_Ak = evaluate(root_Ak, self._shared_eval_cache)
                results.append(det_Ak / det_A)

            return results


def solve_cramer_ddd(Ab: sp.Matrix, symbols: List[sp.Symbol]) -> Dict[sp.Symbol, Any]:
    """
    Convenience matrix equation system solver using LED-DDD and Cramer's rule.
    """
    A = Ab[:, :-1]
    b = Ab[:, -1]
    solver = DDDSolver(A)
    solution = solver.solve(b)
    ret = dict(zip(symbols, solution))
    return ret
