import numpy as np
from polarcodes.Tables import Q_Nmax


class SCD:
    def __init__(self, N, K):
        # Params
        self.N = N
        self.n = int(np.log2(N))
        self.K = K
        # Get the Polar sequence for N
        self.Q_N = Q_Nmax[Q_Nmax < N]
        print("Q_N: {}".format(self.Q_N))
        # DUMMY: Get the information bit pattern.
        self.Q_I = np.ones(N, dtype=int)
        self.F = N - K
        self.Q_I[self.Q_N < self.F] = 0
        print("Q_I: {}".format(self.Q_I))
        # LLRs and partial sums
        self.L = np.zeros((self.N, self.n+1), dtype=float)
        self.s_hat = np.zeros((self.N, self.n), dtype=int)

    def encode(self, a):
        # Frozen bits insertion
        u = np.zeros(self.N, dtype=int)
        u[self.Q_I == 1] = a
        print("u:\n{}".format(u))
        # Polar encoding
        s = np.zeros((self.N, self.n+1), dtype=int)
        s[:, 0] = u
        for l in np.arange(self.n, dtype=int):
            s_block = 1 << l
            i_step = 1 << (l+1)
            for i_block in np.arange(0, self.N, i_step, dtype=int):
                for i in np.arange(s_block):
                    s[i_block + i, l+1] = \
                        s[i_block + i, l] ^ s[i_block + i + s_block, l]
                    s[i_block + i + s_block, l+1] = s[i_block + i + s_block, l]
        # print("s:\n{}".format(s))
        x = s[:, self.n]
        return x

    def decode(self, sbit):
        # Polar decoding (TODO: Partial sum!)
        self.L[:, self.n] = sbit
        for i in np.arange(self.N, dtype=int):
            # LLR and bit levels
            a_llr = active_llr_level(bit_reversed(i, self.n), self.n)
            a_bit = active_bit_level(bit_reversed(i, self.n), self.n)
            # LLR computation
            self.update_llrs(i, a_llr)
            # Decision
            if self.Q_I[i] == 0:
                # Frozen bit
                self.s_hat[i, 0] = 0
            else:
                # Hard decision
                self.s_hat[i, 0] = hard_decision(self.L[i, 0])
            # Partial sum computation
            self.update_bits(i, a_bit)
        a_hat = self.s_hat[self.Q_I == 1, 0]
        return a_hat

    def update_llrs(self, i, a_llr):
        for l in np.arange(a_llr, -1, -1):
            s_block = 1 << l
            for i_b in np.arange(s_block):
                # Select upper or lower computation
                B = (i >> l) & 1
                if B == 0:
                    # Upper (f)
                    self.L[i + i_b, l] = \
                        upper_llr_approx(
                            self.L[i + i_b, l+1],
                            self.L[i + i_b + s_block, l+1])
                else:
                    # Lower (g)
                    self.L[i + i_b, l] = \
                        lower_llr_fix(
                            self.s_hat[i + i_b - s_block, l],
                            self.L[i + i_b - s_block, l+1],
                            self.L[i + i_b, l+1])

    def update_bits(self, i, a_bit):
        for l in np.arange(a_bit):
            s_block = 1 << l
            for i_b in np.arange(s_block):
                self.s_hat[i - i_b - s_block, l+1] = \
                    self.s_hat[i - i_b - s_block, l] ^ self.s_hat[i - i_b, l]
                self.s_hat[i - i_b, l+1] = self.s_hat[i - i_b, l]


def bit_reversed(x, n):
    """
    Bit-reversal operation.

    Parameters
    ----------
    x: ndarray<int>, int
        a vector of indices
    n: int
        number of bits per index in ``x``

    Returns
    ----------
    ndarray<int>, int
        bit-reversed version of x

    """
    result = 0
    for i in range(n):  # for each bit number
        if (x & (1 << i)):  # if it matches that bit
            result |= (1 << (n - 1 - i))  # set the "opposite" bit in result
    return result


def active_llr_level(i, n):
    """
        Find the first 1 in the binary expansion of i.
    """
    mask = 2**(n-1)
    for k in range(n):
        if (mask & i) == 0:
            mask >>= 1
        else:
            break
    return k


def active_bit_level(i, n):
    """
        Find the first 0 in the binary expansion of i.
    """
    mask = 2**(n-1)
    for k in range(n):
        if (mask ^ (i & mask)) == 0:
            mask >>= 1
        else:
            break
    return k


def upper_llr_approx(l1, l2):
    """
    Update top branch LLR in the log-domain (approximation).
    This function supports shortening by checking for infinite LLR cases.

    Parameters
    ----------
    l1: float
        input LLR corresponding to the top branch
    l2: float
        input LLR corresponding to the bottom branch

    Returns
    ----------
    float
        the top branch LLR at the next stage of the decoding tree

    """
    # check for infinite LLR, used in shortening
    if l1 == np.inf and l2 != np.inf:
        return l2
    elif l1 != np.inf and l2 == np.inf:
        return l1
    elif l1 == np.inf and l2 == np.inf:
        return np.inf
    else:  # principal decoding equation
        return np.sign(l1) * np.sign(l2) * np.min((np.abs(l1), np.abs(l2)))


def lower_llr_fix(b, l1, l2):
    """
    Update bottom branch LLR in the log-domain.
    This function supports shortening by checking for infinite LLR cases.

    Parameters
    ----------
    b: int
        the decoded bit of the top branch
    l1: float
        input LLR corresponding to the top branch
    l2: float
        input LLR corresponding to the bottom branch

    Returns
    ----------
    float, np.nan
        the bottom branch LLR at the next stage of the decoding tree
    """
    if b == 0:  # check for infinite LLRs, used in shortening
        if l1 == np.inf or l2 == np.inf:
            return np.inf
        else:
            # g(b, l1, l2) = (-1)^b * l1 + l2
            return l1 + l2
    elif b == 1:
        # g(b, l1, l2) = (-1)^b * l1 + l2
        return -l1 + l2
    else:
        return np.nan


def hard_decision(y):
    """
        Hard decision of a log-likelihood.
    """
    if y >= 0:
        return 0
    else:
        return 1
