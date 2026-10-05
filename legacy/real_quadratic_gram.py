"""
Indecomposables in real quadratic fields and the largest binary (rank <= 2)
totally positive semidefinite Gram matrix carrying distinct square classes of
indecomposables on its diagonal.

This is the corrected version of the experimental code. Two issues from the
original are fixed:

(1) RANK / PSD SOUNDNESS.
    The original accepted a candidate matrix as "rank <= 2, totally positive
    semidefinite" using only: diagonal totally positive, all 2x2 principal
    minors totally >= 0, and all 3x3 PRINCIPAL minors = 0.
    Those conditions are NECESSARY but not SUFFICIENT: rank <= 2 requires ALL
    3x3 minors to vanish (not just principal ones), and principal minors alone
    certify PSD only once PSD is already known. Concretely, the 4x4 correlation
    matrix (3/2)I - (1/2)J (ones on the diagonal, -1/2 off-diagonal) passes all
    three tests yet is rank 4 and indefinite (spectrum {-1/2, 3/2, 3/2, 3/2}).
    So false positives are possible at size >= 4 -- exactly where the
    interesting outputs live.

    FIX: validate each assembled matrix with an EXACT rank check over K,
    rank_K(M) <= 2. This is sufficient for total positive semidefiniteness here:
    valid_off_diagonal already forces every 2x2 principal minor totally >= 0 and
    the diagonal is totally positive, so once rank <= 2 is verified over K every
    principal minor of size >= 3 is identically 0, hence ALL principal minors are
    totally >= 0 at every embedding => totally PSD. No eigenvalues/tolerances
    needed. The cheap 3x3-principal-minor pruning is kept as a (sound) necessary
    condition to keep the search fast; it just no longer serves as acceptance.

(2) SQUARE-CLASS COVERAGE.
    The diagonal candidates were one representative per totally-positive-unit
    orbit. When the fundamental unit has norm +1 (e.g. D = 3, 7, 11, 15, 21, ...)
    each orbit contains TWO distinct square classes of indecomposables,
    [alpha] and [alpha * eta] with eta the totally positive fundamental unit and
    eta not a square. The original only saw [alpha], so it could under-count for
    those fields.

    FIX: square_class_reps() seeds the pool with both alpha and alpha*eta for
    every indecomposable, then dedupes by square class. In norm -1 fields
    eta = epsilon^2 is a square, so alpha*eta is in the same square class as
    alpha and the dedup simply merges it -- harmless. Set expand_square_classes
    =False on largest_gram_matrices to recover the original (restricted) behaviour.

Note on interpretation: largest_gram_matrices computes the maximal size of an
ABSTRACT totally-positive-semidefinite rank-<=2 Gram matrix with these diagonals.
That is a clean UPPER BOUND on how many distinct square classes of indecomposables
a single binary form can represent (a form representing them forces such a matrix
to exist). It does not, by itself, realize the matrix as V^T G V for an honest
binary form G over O_K with integral V.

Indecomposables are computed by the Dress-Scharlau continued-fraction method.

Reference:
    Dress, A.; Scharlau, W. (1975).
    "Indecomposable integral representations of finite groups over real quadratic
    number fields", Advances in Mathematics, 17(3), 231-273.
    Blomer, V.; Kala, V. "On the rank of universal quadratic forms over real
    quadratic fields" (count formula, p. 16).
"""

from sage.all import (
    QuadraticField, RealField, log, Integer, ZZ, continued_fraction,
    matrix, prod
)


class RealQuadraticField:
    """
    Specialized class for computing indecomposables in real quadratic fields and
    the largest binary totally positive semidefinite Gram matrices of distinct
    square classes of indecomposables.
    """

    def __init__(self, D, precision=100):
        """
        Initialize a real quadratic field Q(sqrt(D)).

        Args:
            D: Squarefree positive integer (the radicand)
            precision: Precision for RealField (in bits)
        """
        if not Integer(D).is_squarefree():
            raise ValueError(f"D = {D} must be squarefree")
        if D < 2:
            raise ValueError(f"D = {D} must be greater than 1")

        self.D = D
        self.precision = precision
        self.R = RealField(precision)

        # Number field
        self.K = QuadraticField(D, names='a')
        self.a = self.K.gen()
        self.OK = self.K.ring_of_integers()

        # Field invariants
        self.discriminant = self.K.discriminant()
        self.regulator = self.K.regulator()
        self.class_number = self.K.class_number()

        # Caches
        self._cf = None
        self._delta = None
        self._s = None
        self._fund_unit = None
        self._tp_unit = None
        self._indecomposables = None
        self._alpha_i = None

    # ------------------------------------------------------------------ units

    @property
    def fund_unit(self):
        """Fundamental unit u > 1 of the field."""
        if self._fund_unit is not None:
            return self._fund_unit

        UK = self.K.unit_group()
        u = UK.fundamental_units()[0]

        if abs(u) < 1:
            u = u.trace() - u  # conjugate
        if u < 0:
            u = -u

        assert u > 1, f"Fundamental unit not > 1: {u}"
        assert u.norm() in [1, -1], f"Unit has norm {u.norm()}"

        self._fund_unit = u
        return u

    @property
    def tp_unit(self):
        """Generator eta > 1 of the totally positive units."""
        if self._tp_unit is not None:
            return self._tp_unit

        u = self.fund_unit
        utp = u if u.norm() > 0 else u ** 2  # norm -1 => use u^2

        assert utp.norm() == 1, f"TP unit has norm {utp.norm()}"
        assert utp.trace() > 0, f"TP unit not totally positive: trace = {utp.trace()}"
        assert utp > 1, f"TP unit not > 1: {utp}"

        self._tp_unit = utp
        return utp

    # ----------------------------------------------------- continued fraction

    @property
    def cf_data(self):
        """
        Continued fraction data (cf, delta, s):
            delta = (sqrt(D)+1)/2 if D = 1 (mod 4), else sqrt(D)
            s     = period length
        """
        if self._cf is not None:
            return (self._cf, self._delta, self._s)

        if self.D % 4 == 1:
            cf = continued_fraction((self.a - 1) / 2)
            delta = (self.a + 1) / 2
        else:
            cf = continued_fraction(self.a)
            delta = self.a

        s = len(cf.period())
        self._cf, self._delta, self._s = cf, delta, s
        return (cf, delta, s)

    def _compute_alpha_sequence(self, verbose=False):
        """Sequence alpha_i = p_i + q_i*delta from the convergents p_i/q_i."""
        if self._alpha_i is not None:
            return self._alpha_i

        cf, delta, s = self.cf_data
        num_terms = 2 * s + 10
        alpha_i = [1]  # alpha_0 = 1

        for i in range(num_terms):
            pq = cf.convergent(i)
            p_i, q_i = pq.numer(), pq.denom()
            alpha = self.K(p_i + q_i * delta)
            alpha_i.append(alpha)

            if self.D <= 2000 and len(alpha_i) >= 3:
                expected = cf[i] * alpha_i[-2] + alpha_i[-3]
                assert alpha == expected, \
                    f"Recurrence failed at i={i}: {alpha} != {expected}"

            if self.D <= 2000 and i in [s - 1, 2 * s - 1]:
                assert self.OK(alpha).is_unit(), \
                    f"alpha_{i} = {alpha} is not a unit"
                if i == 2 * s - 1:
                    assert alpha.norm() > 0 and alpha.trace() > 0, \
                        f"alpha_{2*s-1} not totally positive"

        self._alpha_i = alpha_i
        if verbose:
            print("Computing convergents alpha sequence as:", alpha_i)
        return alpha_i

    # -------------------------------------------------------- indecomposables

    def compute_indecomposables_dress_scharlau(self, verbose=True):
        """Indecomposables via Dress-Scharlau (up to totally positive units)."""
        if self._indecomposables is not None:
            return self._indecomposables

        cf, delta, s = self.cf_data
        alpha_i = self._compute_alpha_sequence(verbose=verbose)

        if verbose:
            print(f"Computing indecomposables for Q(sqrt({self.D}))")
            print(f"Discriminant: {self.discriminant}")
            print(f"Period length s = {s}")

        all_indecomposables = []
        for i in range(0, 2 * s + 10, 2):
            cf_i_plus_1 = cf[i + 1]
            for t in range(cf_i_plus_1 + 1):
                alpha_it = self.K(alpha_i[i] + t * alpha_i[i + 1])
                if self.D <= 2000:
                    assert alpha_it.norm() > 0, f"alpha_{i,t} has norm <= 0"
                    assert alpha_it.trace() > 0, f"alpha_{i,t} not totally positive"
                all_indecomposables.append(alpha_it)

        if verbose:
            print(f"Generated {len(all_indecomposables)} candidates")

        indecomposables_final = []
        tp_unit = self.tp_unit
        log_tp = self.R(log(tp_unit).n(max(self.precision, 20 * len(str(tp_unit.trace())))))

        for ind in all_indecomposables:
            is_new = True
            for tmp in indecomposables_final:
                ratio = ind / tmp
                if ratio in self.OK and self.OK(ratio).is_unit():
                    is_new = False
                    break

            if is_new:
                log_ind = self.R(log(ind).n(max(self.precision, 20 * len(str(ind.trace())))))
                l = self.R(((log(ind.norm()) / 2) - log_ind) / log_tp).round()
                good_rep = ind * (tp_unit ** l)
                for cand in [good_rep * tp_unit, good_rep / tp_unit]:
                    if len(str(cand)) < len(str(good_rep)):
                        good_rep = cand
                indecomposables_final.append(good_rep)
                if verbose:
                    print(f"  Found indecomposable: {good_rep} (norm {good_rep.norm()})")

        indecomposables_final.sort(
            key=lambda x: (self.K(x).norm(), self.K(x).trace(), len(str(x)))
        )
        self._indecomposables = indecomposables_final

        if verbose:
            print(f"\nTotal indecomposables: {len(indecomposables_final)}")
            if indecomposables_final:
                max_norm = self.K(indecomposables_final[-1]).norm()
                print(f"Max norm: {max_norm}")
                assert max_norm <= abs(self.discriminant) / 4, \
                    f"Max norm {max_norm} exceeds bound {abs(self.discriminant)/4}"

        return indecomposables_final

    def compute_num_indecomposables(self):
        """Number of indecomposables modulo totally positive units (Blomer-Kala)."""
        if self.D % 4 == 1:
            cf = continued_fraction((1 + self.a) / 2)
        else:
            cf = continued_fraction(self.a)

        cp = cf.period()
        s = cf.period_length()

        if (s % 2) == 0:
            ans = sum(cp[i] for i in range(0, s - 1, 2))
        elif (self.D % 4) == 1:
            ans = 2 * cf[0] + sum(cp[i] for i in range(s - 1)) - 1
        else:
            ans = 2 * cf[0] + sum(cp[i] for i in range(s - 1))
        return ans

    def indecomposables(self):
        """List of indecomposables (one representative per totally positive unit orbit)."""
        if self._indecomposables is None:
            self.compute_indecomposables_dress_scharlau(verbose=False)
        return self._indecomposables

    # --------------------------------------------------------- square classes

    @staticmethod
    def _same_square_class(x, y):
        """True iff x/y is a square in K (equivalently x*y is a square)."""
        is_sq = (x * y).is_square()
        if isinstance(is_sq, tuple):  # some Sage versions return (bool, root)
            is_sq = is_sq[0]
        return bool(is_sq)

    def _dedup_by_square_class(self, pool):
        classes = []
        for ind in pool:
            if not any(self._same_square_class(ind, rep) for rep in classes):
                classes.append(ind)
        return classes

    def get_square_classes(self):
        """
        Square classes among the orbit representatives only (original behaviour).
        Kept for backward compatibility / comparison.
        """
        return self._dedup_by_square_class(self.indecomposables())

    def square_class_reps(self):
        """
        Complete set of square classes of indecomposables.

        For every orbit representative alpha we seed BOTH alpha and alpha*eta
        (eta = tp_unit). The two exhaust the square classes inside the orbit:
        alpha*eta^k depends on k only mod 2. In norm -1 fields eta is a square,
        so alpha*eta merges back into [alpha] under the dedup (no double counting);
        in norm +1 fields it is the genuinely second square class.
        """
        eta = self.tp_unit
        pool = []
        for a in self.indecomposables():
            pool.append(self.K(a))
            pool.append(self.K(a * eta))
        return self._dedup_by_square_class(pool)

    # ------------------------------------------------------- off-diagonals

    def valid_off_diagonal(self, a, b, classical=False):
        """
        All c (in O_K if classical else (1/2)O_K) with a*b - c^2 totally >= 0,
        i.e. the 2x2 principal minor [[a, c], [c, b]] is totally non-negative.

        Floats are used only to bound the integer search box; acceptance is exact.
        """
        K, w = self.K, self.a
        if self.D % 4 == 1:
            w = (1 + self.a) / 2

        emb = K.embeddings(self.R)
        sigma1, sigma2 = emb[0], emb[1]
        a1, a2 = sigma1(a), sigma2(a)
        b1, b2 = sigma1(b), sigma2(b)
        B1 = max(self.R(0), a1 * b1).sqrt()
        B2 = max(self.R(0), a2 * b2).sqrt()
        w1, w2 = sigma1(w), sigma2(w)
        diff_w = abs(w1 - w2)

        D_f = 1 if classical else 2  # denominator: O_K vs (1/2)O_K
        y_max = D_f * (B1 + B2) / diff_w
        y_min = -y_max

        valid_c = []
        for y in range((y_min - 1).floor(), (y_max + 1).ceil() + 1):
            x1_min, x1_max = -D_f * B1 - y * w1,  D_f * B1 - y * w1
            x2_min, x2_max = -D_f * B2 - y * w2,  D_f * B2 - y * w2
            x_min = max(x1_min, x2_min)
            x_max = min(x1_max, x2_max)
            for x in range((x_min - 1).ceil(), (x_max + 1).floor() + 1):
                c_val = (x + y * w) / D_f
                diff = a * b - c_val ** 2
                # totally non-negative: 0, or (trace > 0 and norm >= 0)
                if diff == 0 or (diff.trace() > 0 and diff.norm() >= 0):
                    valid_c.append(c_val)
        return valid_c

    # ---------------------------------------------------------- validation

    def gram_rank(self, M):
        """Exact rank over K of the symmetric matrix given as a list of lists."""
        return matrix(self.K, M).rank()

    def verify_gram_matrix(self, M, tol=1e-9):
        """
        Independent (slower) verification for spot-checking, returning a dict:
          rank_over_K          : exact rank over K (want <= 2)
          min_eig_per_embedding: list of smallest eigenvalues at each real
                                  embedding (want all >= -tol for totally PSD)
        Note: in the main search the exact rank check alone already certifies
        totally PSD (see module docstring); this routine is purely a cross-check.
        """
        rk = self.gram_rank(M)
        min_eigs = []
        for sigma in self.K.embeddings(self.R):
            S = matrix(self.R, [[sigma(e) for e in row] for row in M])
            min_eigs.append(min(S.eigenvalues()))
        return {
            "rank_over_K": rk,
            "min_eig_per_embedding": min_eigs,
            "is_tot_psd_rank_le_2": rk <= 2 and all(e >= -tol for e in min_eigs),
        }

    # ----------------------------------------------------- main computation

    def largest_gram_matrices(self, classical=False, expand_square_classes=True,
                              verbose=True):
        """
        Largest totally positive semidefinite Gram matrices of rank <= 2 with
        distinct square classes of indecomposables on the diagonal.

        Args:
            classical:  off-diagonal entries in O_K (True) or (1/2)O_K (False).
            expand_square_classes:
                        True  -> complete square-class candidate set (the fix);
                        False -> orbit representatives only (original behaviour).
            verbose:    progress output.

        The 3x3-principal-minor (det = 0) condition is used only as a fast,
        SOUND necessary pruning. Acceptance of a matrix is decided by the exact
        rank check rank_K(M) <= 2, which here is sufficient for total positive
        semidefiniteness.
        """
        sq_classes = (self.square_class_reps() if expand_square_classes
                      else self.get_square_classes())
        if verbose:
            mode = "classical" if classical else "non-classical"
            cand = "complete" if expand_square_classes else "orbit-reps-only"
            print(f"Found {len(sq_classes)} distinct square classes "
                  f"({cand}). Mode: {mode}")

        cache = {}

        def get_valid_c(a, b):
            if (a, b) in cache:
                return cache[(a, b)]
            if (b, a) in cache:
                return cache[(b, a)]
            c_list = self.valid_off_diagonal(a, b, classical=classical)
            cache[(a, b)] = c_list
            return c_list

        max_size = 0
        largest_matrices = []

        def search(M, start_idx):
            nonlocal max_size, largest_matrices
            k = len(M)

            # --- acceptance / pruning -------------------------------------
            # k <= 3 is automatically a valid rank-<=2 totally PSD block given
            # diagonal totally positive, 2x2 minors totally >= 0 and the det = 0
            # pruning below. For k >= 4 the principal-minor conditions are no
            # longer sufficient, so verify EXACTLY. If invalid, the whole subtree
            # is invalid too (a principal submatrix of a PSD/rank-<=2 matrix is
            # PSD/rank-<=2), so we prune by returning.
            if k >= 4 and self.gram_rank(M) > 2:
                return

            if k > max_size:
                max_size = k
                largest_matrices = [M]
                if verbose:
                    print(f"Found valid Gram matrix of size {k}")
            elif k == max_size and k > 0:
                largest_matrices.append(M)

            for i in range(start_idx, len(sq_classes)):
                v_new = sq_classes[i]

                if k == 0:
                    search([[v_new]], i + 1)
                    continue

                valid_c_lists = [get_valid_c(M[j][j], v_new) for j in range(k)]

                def backtrack_col(idx, current_c):
                    if idx == k:
                        new_M = [row[:] + [current_c[j]] for j, row in enumerate(M)]
                        new_M.append(list(current_c) + [v_new])
                        search(new_M, i + 1)
                        return

                    for c in valid_c_lists[idx]:
                        valid = True
                        for j in range(idx):
                            v_j = M[j][j]
                            v_idx = M[idx][idx]
                            m_j_idx = M[j][idx]
                            c_j = current_c[j]
                            c_idx = c
                            # 3x3 principal minor {j, idx, new} must vanish
                            # (necessary for rank <= 2).
                            det = (v_j * v_idx * v_new
                                   + 2 * m_j_idx * c_j * c_idx
                                   - v_j * c_idx ** 2
                                   - v_idx * c_j ** 2
                                   - v_new * m_j_idx ** 2)
                            if det != 0:
                                valid = False
                                break
                        if valid:
                            current_c.append(c)
                            backtrack_col(idx + 1, current_c)
                            current_c.pop()

                backtrack_col(0, [])

        search([], 0)

        if verbose:
            print(f"Maximum Gram matrix size: {max_size}")
            print(f"Number of distinct matrices representing this limit: "
                  f"{len(largest_matrices)}")

        return largest_matrices


# ==========================================
# Testing Block (uncomment to execute)
# ==========================================
# if __name__ == "__main__":
#     for D in [2, 7, 41, 58, 214]:
#         print("=" * 42)
#         print(D)
#         print("=" * 42)
#         KK = RealQuadraticField(D)
#         mats = KK.largest_gram_matrices(verbose=True)
#         # spot-check the winners independently:
#         for M in mats:
#             chk = KK.verify_gram_matrix(M)
#             assert chk["is_tot_psd_rank_le_2"], (D, chk)
