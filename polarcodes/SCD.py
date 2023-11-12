import numpy as np


class SCD:
    def __init__(self):
        pass

    def decode(self):
        pass

    def update_llrs(self, l):
        pass

    def update_bits(self, l):
        pass


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


# TODO: The order of these inputs are changed!
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
        else:  # TODO: principal decoding equation - this is propably OK?!
            return l1 + l2
    elif b == 1:  # TODO: principal decoding equation - this is NOT OK?!
        # TODO: Should be "g(b, l1, l2) = (-1)^b * l1 + l2"!
        return l1 - l2
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
