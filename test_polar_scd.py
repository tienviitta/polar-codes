import numpy as np
import argparse

import polarcodes.SCD as PC

np.set_printoptions(precision=2, suppress=True)
# np.set_printoptions(precision=2)

parser = argparse.ArgumentParser()
parser.add_argument("-N", "--l_enc", type=int)
parser.add_argument("-K", "--l_info", type=int)
parser.add_argument("-E", "--esn0", type=float)
parser.add_argument("-S", "--seed", type=int)


def main(N, K, EsN0dB, seed=0):

    # Arguments:
    n = int(np.log2(N))
    print("Polar: N: {}, n: {}, K: {}, EsN0: {} dB, seed: {}".format(
        N, n, K, EsN0dB, seed))

    # Random seed
    np.random.seed(seed)

    # SCD for Polar codes PC(N, K)
    m_scd = PC.SCD(N, K)

    # Info bits
    a = np.random.randint(2, size=K)
    print("a:\n{}".format(a))

    # Encode
    x = m_scd.encode(a)
    print("x:\n{}".format(x))

    # AWGN
    scale = 10**(-EsN0dB / 20)
    awgn = scale * np.random.randn(N)
    y = (1.0 - 2.0 * x) + awgn
    print("y:\n{}, awgn: {:.2f}".format(y, scale))

    # Channel LLRs aka BPSK softbit demodulation and demapper
    norm = 10**(-EsN0dB / 10)
    sbit = (2 / norm) * y
    print("sbit:\n{}, norm: {:.2f}".format(sbit, norm))

    # Decode
    a_hat = m_scd.decode(sbit)

    # Results
    print("L:\n{}".format(m_scd.L))
    print("s_hat:\n{}".format(m_scd.s_hat))
    print("SYMS: {}".format(np.all(np.sign(sbit) == (1.0 - 2.0 * x))))
    print("a:\n{}".format(a))
    print("a_hat:\n{}".format(a_hat))
    print("PASS: {}".format(np.all(a_hat == a)))


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.l_enc, args.l_info, args.esn0, args.seed)
