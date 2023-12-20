import numpy as np
import argparse

import polarcodes.SCLD as PC

np.set_printoptions(precision=2, suppress=True, sign="+")
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

    # # Random seed
    # np.random.seed(seed)

    # SCD for Polar codes PC(N, K)
    m_scld = PC.SCLD(N, K, L=4)

    # # Info bits
    # a = np.random.randint(2, size=K)
    # print("a:\n{}".format(a))
    a = np.array([1, 0, 1, 0, ], dtype=int)

    # # Encode
    # x = m_scld.encode(a)
    # print("x:\n{}".format(x))
    x = np.array(
        [+1.00, -1.00, +1.00, -1.00, -1.00, +1.00, -1.00, +1.00], dtype=int)

    # # AWGN
    # scale = 10**(-EsN0dB / 20)
    # awgn = scale * np.random.randn(N)
    # y = (1.0 - 2.0 * x) + awgn
    # print("y:\n{}, awgn: {:.2f}".format(y, scale))

    # # Channel LLRs aka BPSK softbit demodulation and demapper
    # norm = 10**(-EsN0dB / 10)
    # sbit = (2 / norm) * y
    # print("sbit:\n{}, norm: {:.2f}".format(sbit, norm))

    # Fixed channel LLRs from C-model
    sbit = np.array([-0.26, -0.81, -0.33, -0.71, -0.30, +1.44, -2.13, +1.65])
    print("sbit:\n{}".format(sbit))

    # Decode
    a_hat = m_scld.decode(sbit)

    # Results
    print("L:\n{}".format(m_scld.LLR))
    print("s_hat:\n{}".format(m_scld.BIT))
    print("PM: {}".format(m_scld.PM))
    print("SYMS: {}".format(np.all(np.sign(sbit) == (1.0 - 2.0 * x))))
    print("a:\n{}".format(a))
    print("a_hat:\n{}".format(a_hat))
    print("PASS: {}".format(np.all(a_hat == a)))


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.l_enc, args.l_info, args.esn0, args.seed)
