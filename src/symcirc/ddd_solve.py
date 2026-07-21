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
    Python int arithmetic and hashing are much faster than constructing sorted tuples.
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
from typing import Any, Dict, List, Optional, Sequence, Tuple, Set, Union

import sympy as sp
from sympy.polys.polyerrors import PolynomialError

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
        string = "+" if self.sign >= 0 else "-"
        return f"DDD[{string}a[{self.abs_row},{self.abs_col}]]"


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
    elif node.sign == 1:
        if node_is_one:
            cofactor = one_val
        elif one_edge_is_one:
            cofactor = node.symbol
        else:
            cofactor = sp.Mul(node.symbol, one_val, evaluate=False)
    else:
        if node_is_one:
            cofactor = -one_val
        elif one_edge_is_one:
            cofactor = -node.symbol  # same
        else:
            cofactor = sp.Mul(-node.symbol, one_val, evaluate=False)
    # Combine with remainder, avoiding explicit zero add (e.g. a+0)
    if zero_val == 0:
        result = cofactor
    elif cofactor == 0:
        result = zero_val
    else:
        result = sp.Add(zero_val, cofactor, evaluate=False)

    cache[nid] = result
    return result

_NUMERIC_TYPES = (int, float, sp.Float, sp.Integer, sp.Rational)


def _resolve_symbol_value(symbol: Any,
                          value_map: Dict[sp.Symbol, float],
                          keep_symbolic: Set[sp.Symbol],
                          resolve_cache: Dict[int, Tuple[Any, bool]]) -> Tuple[Any, bool]:
    """
    Substitutes known numeric values into a single matrix-entry symbol/expression,
    returning (substituted_value, is_purely_numeric).
    """
    key = id(symbol)
    cached = resolve_cache.get(key)
    if cached is not None:
        return cached

    if isinstance(symbol, sp.Symbol):
        if symbol in keep_symbolic:
            val = symbol
        elif symbol in value_map:
            val = sp.Float(value_map[symbol])
        else:
            val = symbol
    elif isinstance(symbol, _NUMERIC_TYPES):
        val = symbol
    else:
        free = symbol.free_symbols if hasattr(symbol, 'free_symbols') else set()
        if free & keep_symbolic:
            val = symbol
        elif free:
            subs_dict = {symb: value_map[symb] for symb in free
                         if symb in value_map and symb not in keep_symbolic}
            if subs_dict:
                val = symbol.subs(subs_dict)
                if not val.free_symbols:
                    val = sp.Float(float(val))
            else:
                val = symbol
        else:
            val = symbol

    is_numeric = isinstance(val, _NUMERIC_TYPES) or (hasattr(val, 'is_number') and bool(val.is_number))
    result = (val, is_numeric)
    resolve_cache[key] = result
    return result


def _evaluate_semi_symbolic(node: Any,
                            value_map: Dict[sp.Symbol, float],
                            keep_symbolic: Set[sp.Symbol],
                            cache: Dict[int, Tuple[Any, bool]],
                            resolve_cache: Dict[int, Tuple[Any, bool]]) -> Tuple[Any, bool]:
    """
    Core evaluator. Returns (value, is_purely_numeric).
    is_purely_numeric is tracked incrementally.
    """
    if node is TERMINAL_ZERO:
        return sp.S.Zero, True
    if node is TERMINAL_ONE:
        return sp.S.One, True

    if not isinstance(node, DDDNode):
        if isinstance(node, sp.Symbol):
            if node in keep_symbolic:
                return node, False
            elif node in value_map:
                return sp.Float(value_map[node]), True
            else:
                return node, False
        return node, True

    nid = node.id
    hit = cache.get(nid)
    if hit is not None:
        return hit

    zero_val, zero_numeric = _evaluate_semi_symbolic(node.zero, value_map, keep_symbolic, cache, resolve_cache)
    one_val, one_numeric = _evaluate_semi_symbolic(node.one, value_map, keep_symbolic, cache, resolve_cache)
    symbol_val, symbol_numeric = _resolve_symbol_value(node.symbol, value_map, keep_symbolic, resolve_cache)

    # Avoid constructing a Mul with an explicit zero factor.
    if symbol_val == 0 or one_val == 0:
        cofactor = sp.S.Zero
        cofactor_numeric = True
    else:
        cofactor_numeric = symbol_numeric and one_numeric
        # Numeric branches can be safely auto-evaluated (cheap, keeps them collapsed
        # to a single Float); symbolic branches keep evaluate=False to preserve the
        # compact factored form instead of expanding it.
        cofactor = sp.Mul(node.sign, symbol_val, one_val, evaluate=cofactor_numeric)

    if zero_val == 0:
        result = cofactor
        result_numeric = cofactor_numeric
    elif cofactor == 0:
        result = zero_val
        result_numeric = zero_numeric
    else:
        result_numeric = cofactor_numeric and zero_numeric
        result = sp.Add(zero_val, cofactor, evaluate=result_numeric)

    if result_numeric and isinstance(result, (sp.Add, sp.Mul)):
        result = sp.Float(float(result))

    ret = (result, result_numeric)
    cache[nid] = ret
    return ret


def evaluate_semi_symbolic(node: Any,
                           value_map: Dict[sp.Symbol, float],
                           keep_symbolic: Set[sp.Symbol] = None,
                           cache: Dict = None,
                           resolve_cache: Dict = None) -> Union[float, sp.Expr]:
    """
    Evaluate a DDD to its value, substituting `value_map` into every symbol except
    those in `keep_symbolic`. Thin backwards-compatible wrapper around the
    (value, is_numeric) core evaluator - see `_evaluate_semi_symbolic`.
    """
    if cache is None:
        cache = {}
    if keep_symbolic is None:
        keep_symbolic = set()
    if resolve_cache is None:
        resolve_cache = {}
    value, _ = _evaluate_semi_symbolic(node, value_map, keep_symbolic, cache, resolve_cache)
    return value


def evaluate_numeric_direct(node: Any, value_map: Dict[sp.Symbol, float],
                            cache: Dict = None, expr_cache: Dict = None) -> float:
    """
    Evaluate DDD directly to numeric value without creating symbolic expression.

    This bypasses ALL symbolic manipulation and is the fastest approach.

    Args:
        node: DDD root node
        value_map: {symbol: numeric_value}
        cache: Evaluation cache (uses node._uid as key)
        expr_cache: Cache for evaluated expressions

    Returns:
        Numeric determinant value

    Time complexity: O(DDD nodes) ≈ O(n²) for sparse n×n matrix
    """
    if cache is None:
        cache = {}
    if expr_cache is None:
        expr_cache = {}

    # Terminal cases
    if node is TERMINAL_ZERO:
        return 0.0
    if node is TERMINAL_ONE:
        return 1.0
    if not isinstance(node, DDDNode):
        # Constant
        if isinstance(node, sp.Symbol):
            return float(value_map.get(node, 0.0))
        return float(node)
    nid = node.id
    # Check cache
    if nid in cache:
        return cache[nid]

    # Recursively evaluate children
    zero_val = evaluate_numeric_direct(node.zero, value_map, cache, expr_cache)
    one_val = evaluate_numeric_direct(node.one, value_map, cache, expr_cache)

    # Get numeric value of symbol
    symbol = node.symbol

    # Check expression cache first (keyed by identity: the same entry object
    # recurs across many DDDNode instances, so this avoids re-substituting it
    # and avoids the cost of str() it just to build a cache key)
    symbol_key = id(symbol)
    if symbol_key in expr_cache:
        symbol_val = expr_cache[symbol_key]
    else:
        if isinstance(symbol, sp.Symbol):
            symbol_val = value_map[symbol]
        elif isinstance(symbol, (int, float)):
            symbol_val = float(symbol)
        else:
            # SymPy expression - substitute and cache
            symbol_val = float(symbol.subs(value_map))
        expr_cache[symbol_key] = symbol_val

    result = zero_val + node.sign * symbol_val * one_val

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


class _SparsePoly:
    """
    Custom multivariate polynomial: a dict mapping an exponent tuple
    (one exponent per gen, in a fixed shared order) to a float coefficient.

    This exists to replace sympy Add/Mul/Poly during DDD evaluation. It is
    only valid for building up polynomials - it supports
    nothing else (no division, no non-polynomial terms) - which is
    the structure semisymbolic circuit-matrix entries in LTI have.
    Because coefficients are combined by plain dict
    accumulation, a _SparsePoly is always already in fully collected
    canonical form, so there is no need to simplify.
    """
    __slots__ = ("c",)

    def __init__(self, c: Dict[Tuple[int, ...], float]):
        self.c = c

    def add(self, other: Union["_SparsePoly", float], ngens: int) -> "_SparsePoly":
        if isinstance(other, _SparsePoly):
            other_c = other.c
        else:
            if other == 0:
                return self
            other_c = {(0,) * ngens: float(other)}
        result = dict(self.c)
        for k, v in other_c.items():
            nv = result.get(k, 0.0) + v
            if nv == 0.0:
                result.pop(k, None)
            else:
                result[k] = nv
        return _SparsePoly(result)

    def mul(self, other: Union["_SparsePoly", float], ngens: int) -> "_SparsePoly":
        if isinstance(other, _SparsePoly):
            if not self.c or not other.c:
                return _SparsePoly({})
            result: Dict[Tuple[int, ...], float] = {}
            for k1, v1 in self.c.items():
                for k2, v2 in other.c.items():
                    k = tuple(a + b for a, b in zip(k1, k2))
                    result[k] = result.get(k, 0.0) + v1 * v2
            return _SparsePoly(result)
        else:
            if other == 0 or not self.c:
                return _SparsePoly({})
            fo = float(other)
            return _SparsePoly({k: v * fo for k, v in self.c.items()})

    def to_expr(self, gens: Tuple[sp.Symbol, ...]) -> sp.Expr:
        expr = sp.Integer(0)
        for exps, coeff in self.c.items():
            term = sp.Float(coeff)
            for g, e in zip(gens, exps):
                if e:
                    term = term * g ** e
            expr = expr + term
        return expr

    def to_sympy_poly(self, gens: Tuple[sp.Symbol, ...]) -> sp.Poly:
        return sp.Poly.from_dict({k: sp.Float(v) for k, v in self.c.items()}, *gens)


def _collect_poly_gens(A: sp.Matrix, b: Sequence[Any],
                       value_map: Dict[sp.Symbol, float],
                       keep_symbolic: Set[sp.Symbol]) -> Tuple[sp.Symbol, ...]:
    """
    Determines the fixed set of generators (typically just `s`, plus any
    independent-source symbols the "ddd" solver's construction hack leaves in
    `b`) by looking at the small, bounded set of *distinct* matrix/vector
    entries - not at the (much larger) set of DDD nodes built from them.
    """
    seen_ids: Set[int] = set()
    gens: Set[sp.Symbol] = set()
    for entry in list(A) + list(b):
        if id(entry) in seen_ids:
            continue
        seen_ids.add(id(entry))
        free = entry.free_symbols if hasattr(entry, 'free_symbols') else set()
        for sym in free:
            if sym not in value_map or sym in keep_symbolic:
                gens.add(sym)
    return tuple(sorted(gens, key=str))


def _resolve_symbol_poly(symbol: Any,
                         value_map: Dict[sp.Symbol, float],
                         keep_symbolic: Set[sp.Symbol],
                         gens: Tuple[sp.Symbol, ...],
                         resolve_cache: Dict[int, Tuple[Any, bool]]) -> Tuple[Any, bool]:
    """
    Like `_resolve_symbol_value`, but returns a `_SparsePoly` instead of a
    sympy Expr for the symbolic case. Raises `sp.PolynomialError` if an entry
    turns out not to actually be polynomial in `gens` (e.g. a distributed
    element contributing exp(-s*tau)) - callers should catch this once, at
    the top of a solve, and fall back to the generic Expr-based evaluator.
    """
    key = id(symbol)
    cached = resolve_cache.get(key)
    if cached is not None:
        return cached

    free = symbol.free_symbols if hasattr(symbol, 'free_symbols') else set()
    relevant = {sy for sy in free if sy not in value_map or sy in keep_symbolic}

    if not relevant:
        subs_dict = {sy: value_map[sy] for sy in free if sy in value_map}
        val = float(symbol.subs(subs_dict)) if free else float(symbol)
        result = (val, True)
    else:
        subs_dict = {sy: value_map[sy] for sy in free if sy in value_map and sy not in keep_symbolic}
        val = symbol.subs(subs_dict) if subs_dict else symbol
        poly = sp.Poly(val, *gens)  # raises PolynomialError if not polynomial in gens
        result = (_SparsePoly({k: float(v) for k, v in poly.as_dict().items()}), False)

    resolve_cache[key] = result
    return result


def _evaluate_semi_symbolic_poly(node: Any,
                                 value_map: Dict[sp.Symbol, float],
                                 keep_symbolic: Set[sp.Symbol],
                                 gens: Tuple[sp.Symbol, ...],
                                 cache: Dict[int, Tuple[Any, bool]],
                                 resolve_cache: Dict[int, Tuple[Any, bool]]) -> Tuple[Any, bool]:
    """
    Same recursion/caching structure as `_evaluate_semi_symbolic`, but
    accumulates plain Python floats / `_SparsePoly` objects instead of sympy
    Add/Mul(evaluate=False) trees. No sympy object is created anywhere in the
    recursive path - only plain dict/float arithmetic - and the result is
    already in canonical collected form, so no expand/collect/cancel pass is
    needed afterwards.
    """
    ngens = len(gens)

    if node is TERMINAL_ZERO:
        return 0.0, True
    if node is TERMINAL_ONE:
        return 1.0, True

    if not isinstance(node, DDDNode):
        if isinstance(node, sp.Symbol):
            if node in keep_symbolic:
                exps = tuple(1 if g == node else 0 for g in gens)
                return _SparsePoly({exps: 1.0}), False
            elif node in value_map:
                return float(value_map[node]), True
            else:
                exps = tuple(1 if g == node else 0 for g in gens)
                return _SparsePoly({exps: 1.0}), False
        return float(node), True

    nid = node.id
    hit = cache.get(nid)
    if hit is not None:
        return hit

    zero_val, zero_numeric = _evaluate_semi_symbolic_poly(node.zero, value_map, keep_symbolic, gens, cache, resolve_cache)
    one_val, one_numeric = _evaluate_semi_symbolic_poly(node.one, value_map, keep_symbolic, gens, cache, resolve_cache)
    symbol_val, symbol_numeric = _resolve_symbol_poly(node.symbol, value_map, keep_symbolic, gens, resolve_cache)

    if (symbol_numeric and symbol_val == 0) or (one_numeric and one_val == 0):
        cofactor, cofactor_numeric = 0.0, True
    else:
        signed_symbol = symbol_val * node.sign if symbol_numeric else symbol_val.mul(node.sign, ngens)
        if symbol_numeric and one_numeric:
            cofactor, cofactor_numeric = signed_symbol * one_val, True
        elif symbol_numeric:
            cofactor, cofactor_numeric = one_val.mul(signed_symbol, ngens), False
        else:
            cofactor, cofactor_numeric = signed_symbol.mul(one_val, ngens), False

    if zero_numeric and zero_val == 0:
        result, result_numeric = cofactor, cofactor_numeric
    elif cofactor_numeric and cofactor == 0:
        result, result_numeric = zero_val, zero_numeric
    elif zero_numeric and cofactor_numeric:
        result, result_numeric = zero_val + cofactor, True
    elif zero_numeric:
        result, result_numeric = cofactor.add(zero_val, ngens), False
    elif cofactor_numeric:
        result, result_numeric = zero_val.add(cofactor, ngens), False
    else:
        result, result_numeric = zero_val.add(cofactor, ngens), False

    ret = (result, result_numeric)
    cache[nid] = ret
    return ret


class SemiSymSolver:
    """
    Solver for Cramer's-rule based semi-symbolic solving of a linear system
    (component values substituted in numerically, chosen symbols such as `s`
    left symbolic).

    Cramer's rule computes each unknown independently as det(A_k)/det(A), so
    `solve()` accepts an optional `targets` subset of unknowns to build/evaluate -
    asking for a single output skips all the unrelated numerator work entirely
    instead of paying for every unknown in the system.

    For LTI circuits, every matrix entry stamped by an R/L/C/linearized
    controlled source is a low-degree polynomial of 's'.
    `solve()` exploits that by accumulating determinants as non-sympy
    polynomials (`_SparsePoly`, plain dict-of-floats arithmetic)
    instead of sympy expression trees, which is both much faster than building
    an unevaluated Add/Mul tree and already fully collected.
    If any entry turns out not to be polynomial in the detected generators,
    it falls back to the general symbolic evaluator.
    """
    def __init__(self, A: sp.Matrix, strategy: str = "min_degree"):
        self.A = A
        self.strategy = strategy
        self.denominator_builder = None
        self._shared_eval_cache: Dict[int, Tuple[Any, bool]] = {}
        self._shared_resolve_cache: Dict[int, Tuple[Any, bool]] = {}

    def solve(self, b: Sequence[Any], value_map: Dict[sp.Symbol, float], keep_symbolic: Set[sp.Symbol],
              targets: Optional[Sequence[int]] = None, simplify: bool = False,
              fast_poly: bool = True) -> List[Any]:
        """
        :param targets: optional subset of unknown indices (0-indexed, matching
            the columns of A) to actually solve for. If None (default), solves
            for every unknown, matching the original behavior. Entries not in
            `targets` are left as None in the returned list.
        :param simplify: request a canonical (fully reduced) result. When the
            fast polynomial path is used (the default), the result is already
            fully collected as soon as it's built, so this only adds a cheap
            exact GCD-based reduction to strip a common factor
            between numerator and denominator. Matters for exact pole/zero extraction.
            If the generic fallback path is used instead, this runs sympy.cancel(),
            which is far more expensive since it has to expand the unevaluated
            expression tree before it can search for a common factor.
        :param fast_poly: use the custom polynomial evaluator.
            Good whenever every matrix/b entry is a polynomial in `s`.
            Automatically falls back to the general symbolic evaluator otherwise,
            so turning this off is only useful for debugging/comparison.
        """
        if keep_symbolic is None:
            keep_symbolic = set()
        n = self.A.shape[0]
        target_indices = range(n) if targets is None else list(targets)
        b_list = [sp.sympify(x) for x in b]
        results: List[Any] = [None] * n

        if fast_poly:
            try:
                gens = _collect_poly_gens(self.A, b_list, value_map, keep_symbolic)
                poly_eval_cache: Dict[int, Tuple[Any, bool]] = {}
                poly_resolve_cache: Dict[int, Tuple[Any, bool]] = {}

                self.denominator_builder = DDDBuilder(self.A, self.strategy)
                root_A = self.denominator_builder.build()
                det_A, det_A_numeric = _evaluate_semi_symbolic_poly(
                    root_A, value_map, keep_symbolic, gens, poly_eval_cache, poly_resolve_cache)

                for k in target_indices:
                    numerator_builder = DDDBuilder(self.A, b_vec=b_list, k=k, reference=self.denominator_builder)
                    root_Ak = numerator_builder.build()
                    det_Ak, det_Ak_numeric = _evaluate_semi_symbolic_poly(
                        root_Ak, value_map, keep_symbolic, gens, poly_eval_cache, poly_resolve_cache)

                    if simplify and not det_Ak_numeric and not det_A_numeric:
                        # Both sides are already-collected native polys - build
                        # sympy Poly objects directly from the coefficient dicts
                        # (cheap: no expand needed) and cancel a genuine common
                        # factor, if any, via GCD.
                        num_p = det_Ak.to_sympy_poly(gens)
                        den_p = det_A.to_sympy_poly(gens)
                        g = sp.gcd(num_p, den_p)
                        if not g.is_one:
                            num_p = num_p.quo(g)
                            den_p = den_p.quo(g)
                        ratio = num_p.as_expr() / den_p.as_expr()
                    else:
                        num_expr = det_Ak if det_Ak_numeric else det_Ak.to_expr(gens)
                        den_expr = det_A if det_A_numeric else det_A.to_expr(gens)
                        ratio = num_expr / den_expr
                        if simplify:
                            ratio = sp.cancel(ratio)
                    results[k] = ratio
                return results
            except (PolynomialError, TypeError, ValueError, AttributeError):
                # Some entry isn't polynomial in the detected generators. Fall
                # through to the general evaluator below - none of the
                # fast-path caches/results are reused past this point.
                results = [None] * n

        # General fallback: unevaluated symbolic Add/Mul tree, safe for any
        # entry regardless of whether it's polynomial in s.
        self.denominator_builder = DDDBuilder(self.A, self.strategy)
        root_A = self.denominator_builder.build()
        det_A, _ = _evaluate_semi_symbolic(root_A, value_map, keep_symbolic,
                                           self._shared_eval_cache, self._shared_resolve_cache)

        for k in target_indices:
            numerator_builder = DDDBuilder(self.A, b_vec=b_list, k=k, reference=self.denominator_builder)
            root_Ak = numerator_builder.build()
            det_Ak, _ = _evaluate_semi_symbolic(root_Ak, value_map, keep_symbolic,
                                                self._shared_eval_cache, self._shared_resolve_cache)
            ratio = det_Ak / det_A
            if simplify:
                ratio = sp.cancel(ratio)
            results[k] = ratio
        return results


class DDDSolver:
    """
    Solver for standard determinant and Cramer's Rule subdeterminants (with reference cache sharing).
    """
    def __init__(self, A: sp.Matrix, strategy: str = "min_degree"):
        self.A = A
        self.strategy = strategy
        self.denominator_builder = None
        self._shared_eval_cache = {}
        self._shared_expr_cache = {}

    def solve_num(self, b, value_dict: Dict[sp.Symbol, float]) -> Any:
        self.denominator_builder = DDDBuilder(self.A, self.strategy)
        root_A = self.denominator_builder.build()

        det_A = evaluate_numeric_direct(root_A, value_dict, self._shared_eval_cache, self._shared_expr_cache)
        # Build Numerators det(Ak)
        b_list = [sp.sympify(x) for x in b]
        results = []

        for k in range(self.A.shape[0]):
            numerator_builder = DDDBuilder(self.A, b_vec=b_list, k=k, reference=self.denominator_builder)
            root_Ak = numerator_builder.build()
            det_Ak = evaluate_numeric_direct(root_Ak, value_dict, self._shared_eval_cache, self._shared_expr_cache)
            results.append(det_Ak / det_A)

        return results

    def solve(self, b: Optional[Sequence[Any]] = None, nested: bool = True,
              targets: Optional[Sequence[int]] = None) -> Any:
        """
        :param targets: optional subset of row indices to solve for. Cramer's
            rule computes each unknown independently, so restricting to the
            ones you need skips the unrelated numerator work - for a
            single target on an n-unknown system this is roughly an n-times
            speedup. If None (default), solves for every unknown.
        """
        # Build Denominator det(A)
        self.denominator_builder = DDDBuilder(self.A, self.strategy)
        root_A = self.denominator_builder.build()
        if nested: det_A = evaluate_nested(root_A, self._shared_eval_cache)
        else: det_A = evaluate(root_A, self._shared_eval_cache)

        # Build only the requested Numerators det(Ak)
        b_list = [sp.sympify(x) for x in b]
        n = self.A.shape[0]
        target_indices = range(n) if targets is None else list(targets)
        results: List[Any] = [None] * n

        for k in target_indices:
            numerator_builder = DDDBuilder(self.A, b_vec=b_list, k=k, reference=self.denominator_builder)
            root_Ak = numerator_builder.build()
            if nested: det_Ak = evaluate_nested(root_Ak, self._shared_eval_cache)
            else: det_Ak = evaluate(root_Ak, self._shared_eval_cache)
            # Plain `det_Ak / det_A` goes through Expr.__truediv__ with
            # evaluate=True, which runs sympy's full assumption system
            # (is_zero/is_infinite/... deduction) over these large unevaluated
            # trees before it can even build the Mul/Pow - construcing the
            # unevaluated ratio directly skips all of that.
            results[k] = sp.Mul(det_Ak, sp.Pow(det_A, -1, evaluate=False), evaluate=False)
        return results

    def solve_node(self, b: Optional[Sequence[Any]] = None, nested: bool=True):
        pass

    def det(self, nested: bool=True) -> Any:
        self.denominator_builder = DDDBuilder(self.A, self.strategy)
        root_A = self.denominator_builder.build()
        if nested:
            det_A = evaluate_nested(root_A, self._shared_eval_cache)
        else:
            det_A = evaluate(root_A, self._shared_eval_cache)
        return det_A


def solve_cramer_ddd(Ab: sp.Matrix, symbols: List[sp.Symbol], nested: bool = True,
                     targets: Optional[Sequence[sp.Symbol]] = None) -> Dict[sp.Symbol, Any]:
    """
    Convenience matrix equation system solver using LED-DDD and Cramer's rule.

    :param targets: optional subset of `symbols` to actually solve for. See
        `DDDSolver.solve`.
    """
    A = Ab[:, :-1]
    b = Ab[:, -1]
    solver = DDDSolver(A)

    if targets is not None:
        target_indices = [symbols.index(sym) for sym in targets]
    else:
        target_indices = None

    solution = solver.solve(b, nested=nested, targets=target_indices)

    if targets is not None:
        ret = {sym: solution[idx] for sym, idx in zip(targets, target_indices)}
    else:
        ret = dict(zip(symbols, solution))
    return ret

def solve_cramer_ddd_semi(Ab: sp.Matrix, symbols: List[sp.Symbol], value_map: Dict[sp.Symbol, float],
                          keep_symbolic: Optional[Set[sp.Symbol]] = None,
                          targets: Optional[Sequence[sp.Symbol]] = None,
                          simplify: bool = False,
                          fast_poly: bool = True) -> Dict[sp.Symbol, Any]:
    """
    :param targets: optional subset of unknown indices (0-indexed, matching
        the columns of A) to actually solve for. If None (default), solves
        for every unknown, matching the original behavior. Entries not in
        `targets` are left as None in the returned list.
    :param simplify: request a canonical (fully reduced) result. When the
        fast polynomial path is used (the default), the result is already
        fully collected as soon as it's built, so this only adds a cheap
        exact GCD-based reduction to strip a common factor
        between numerator and denominator. Matters for exact pole/zero extraction.
        If the generic fallback path is used instead, this runs sympy.cancel(),
        which is far more expensive since it has to expand the unevaluated
        expression tree before it can search for a common factor.
    :param fast_poly: use the custom polynomial evaluator.
        Good whenever every matrix/b entry is a polynomial in `s`.
        Automatically falls back to the general symbolic evaluator otherwise,
        so turning this off is only useful for debugging/comparison.
    """
    #for key in value_map:
    #    print(f"{key}: {value_map[key]}")
    A = Ab[:, :-1]
    b = Ab[:, -1]
    solver = SemiSymSolver(A)

    if targets is not None:
        target_indices = [symbols.index(sym) for sym in targets]
    else:
        target_indices = None

    solution = solver.solve(b, value_map, keep_symbolic, targets=target_indices,
                            simplify=simplify, fast_poly=fast_poly)

    if targets is not None:
        ret = {sym: solution[idx] for sym, idx in zip(targets, target_indices)}
    else:
        ret = dict(zip(symbols, solution))
    return ret

def solve_cramer_ddd_numeric(Ab: sp.Matrix, symbols: List[sp.Symbol],
                             value_dict: Dict[sp.Symbol, float]=True) -> Dict[sp.Symbol, Any]:
    A = Ab[:, :-1]
    b = Ab[:, -1]
    solver = DDDSolver(A)
    solution = solver.solve_num(b, value_dict)
    ret = dict(zip(symbols, solution))
    return ret