import numpy as np
import argparse

np.set_printoptions(precision=4, suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument("-N", "--length", type=int)


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
    # count = 1
    for k in range(n):
        if (mask & i) == 0:
            # count += 1
            mask >>= 1
        else:
            break
    # return min(count, n)
    return k


def active_bit_level(i, n):
    """
        Find the first 0 in the binary expansion of i.
    """

    mask = 2**(n-1)
    # count = 1
    for k in range(n):
        if (mask ^ (i & mask)) == 0:
            # count += 1
            mask >>= 1
        else:
            break
    # return min(count, n)
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


def lower_llr(b, l2, l1):
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
        else:  # principal decoding equation
            return l1 + l2
    elif b == 1:  # principal decoding equation
        return l1 - l2
    else:
        return np.nan


def main(N):
    # Arguments:
    n = int(np.log2(N))
    print("Polar: N:{}, n:{}".format(N, n))

    # Random seed
    np.random.seed(654321)

    # Polar sequence (TODO: Now fixed for N=8!)
    Q_N = np.array([0, 1, 2, 4, 3, 5, 6, 7], dtype=int)
    Q_I = np.array([0, 0, 0, 1, 0, 1, 1, 1], dtype=int)
    n_a = np.sum(Q_I)

    # Info bits
    # a = np.random.randint(2, size=n_a)
    a = np.ones(n_a, dtype=int)
    print("a:\n{}".format(a))

    # Frozen bits
    u = np.zeros(N, dtype=int)
    u[Q_I == 1] = a
    print("u:\n{}".format(u))

    # Polar encoding
    s = np.zeros((N, n+1), dtype=int)
    s[:, 0] = u
    for l in np.arange(n, dtype=int):
        n_block = N >> (l+1)
        s_block = 1 << l
        # print("l: {}, n_b: {}, s_b: {}".format(l, n_block, s_block))
        for i_block in np.arange(n_block, dtype=int):
            # print("  i_b: {}".format(i_block))
            for i_bly in np.arange(s_block, dtype=int):
                i_L = (i_block << s_block) + i_bly
                s[i_L, l+1] = s[i_L, l] ^ s[i_L + s_block, l]
                s[i_L + s_block, l+1] = s[i_L + s_block, l]
                # print("   i: {}".format(i))
    print("s:\n{}".format(s))
    x = s[:, n]
    print("x:\n{}".format(x))

    # AWGN
    EsN0dB = 3
    scale = 10**(-EsN0dB / 20)
    awgn = np.random.randn(N)
    y = (1.0 - 2.0 * x) + (scale * awgn)
    print("y:\n{}, scale: {:.2f}".format(y, scale))

    # Channel LLRs
    norm = 10**(-EsN0dB / 10)
    sbit = (2 / norm) * y
    print("sbit:\n{}, norm: {:.2f}".format(sbit, norm))

    # Polar decoding (TODO: Partial sum!)
    L = np.zeros((N, n+1), dtype=float)
    L[:, n] = sbit
    print("L:\n{}".format(L))
    s_hat = -np.ones((N, n+1), dtype=int)
    print("s_hat:\n{}".format(s_hat))
    for i in np.arange(N, dtype=int):
        a_llr = active_llr_level(bit_reversed(i, n), n)
        a_bit = active_bit_level(bit_reversed(i, n), n)
        # print("a_llr: {}, a_bit: {}".format(a_llr, a_bit))
        for l in np.arange(a_llr, -1, -1):
            s_block = 1 << l
            for i_b in np.arange(s_block):
                # B = (i // (1 << l)) % 2
                # B = (i >> l) % 2
                B = (i >> l) & 1
                # print("  i: {}, l: {}, s_block: {}, B: {}".format(
                #     i, l, s_block, B))
                if B == 0:
                    L[i + i_b, l] = \
                        upper_llr_approx(
                            L[i + i_b, l+1], L[i + i_b + s_block, l+1])
                else:
                    L[i + i_b, l] = \
                        lower_llr(
                            s[i + i_b - s_block, l], L[i + i_b - s_block, l+1], L[i + i_b, l+1])

    print("L:\n{}".format(L))
    u_hat = np.array((1-np.sign(L[:, 0]))//2, dtype=int)
    print("u:\n{}".format(u))
    print("u_hat:\n{}".format(u_hat))
    print("PASS: {}".format(np.all(u_hat == u)))
    print("SYMS: {}".format(np.all(np.sign(sbit) == (1.0 - 2.0 * x))))


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.length)
