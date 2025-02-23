import numpy as np
import argparse
from polarcodes import utils
from polarcodes import decoder_utils

np.set_printoptions(precision=4, suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument("-A", "--n_info", type=int)
parser.add_argument("-E", "--n_rm", type=int)
parser.add_argument("-N", "--esn0", type=float)


def main(A, E, EsN0):
    # Random seed
    np.random.seed(123456)

    # The CRC polynomial used in 3GPP PBCH and PDCCH channel is
    # D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1
    crc_polynomial_pattern = np.array(
        [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1])
    P = crc_polynomial_pattern.size - 1
    K = A + P
    print("A: {}, E: {}, P: {}, K: {}".format(A, E, P, K))

    # Determine the number of bits used at the input and output of the polar
    # encoder kernel.
    n_max = int(9)  # PBCH and PDCCH
    if E <= (9/8)*(2**(np.ceil(np.log2(E))-1)) and (K/E) < (9/16):
        n_1 = np.ceil(np.log2(E))-1
    else:
        n_1 = np.ceil(np.log2(E))
    R_min = 1/8
    n_min = 5
    n_2 = np.ceil(np.log2(K/R_min))
    n = int(np.max((n_min, np.min([n_1, n_2, n_max]))))
    N = 2 ** n
    print("n: {}, N: {}".format(n, N))

    # Get the 3GPP CRC interleaver pattern.
    Pi_IL_max = np.array(
        [0, 2, 3, 5, 6, 8, 11, 12, 13, 16, 19, 20, 22, 24, 28, 32, 33, 35, 37, 38,
         39, 40, 41, 42, 44, 46, 47, 49, 50, 54, 55, 57, 59, 60, 62, 64, 67, 69, 74,
         79, 80, 84, 85, 86, 88, 91, 94, 102, 105, 109, 110, 111, 113, 114, 116, 118,
         119, 121, 122, 125, 126, 127, 129, 130, 131, 132, 136, 137, 141, 142, 143,
         147, 148, 149, 151, 153, 155, 158, 161, 164, 166, 168, 170, 171, 173, 175,
         178, 179, 180, 182, 183, 186, 187, 189, 192, 194, 198, 199, 200, 1, 4, 7,
         9, 14, 17, 21, 23, 25, 29, 34, 36, 43, 45, 48, 51, 56, 58, 61, 63, 65, 68,
         70, 75, 81, 87, 89, 92, 95, 103, 106, 112, 115, 117, 120, 123, 128, 133,
         138, 144, 150, 152, 154, 156, 159, 162, 165, 167, 169, 172, 174, 176, 181,
         184, 188, 190, 193, 195, 201, 10, 15, 18, 26, 30, 52, 66, 71, 76, 82, 90,
         93, 96, 104, 107, 124, 134, 139, 145, 157, 160, 163, 177, 185, 191, 196,
         202, 27, 31, 53, 72, 77, 83, 97, 108, 135, 140, 146, 197, 203, 73, 78, 98,
         204, 99, 205, 100, 206, 101, 207, 208, 209, 210, 211, 212, 213, 214, 215,
         216, 217, 218, 219, 220, 221, 222, 223], dtype=int)
    crc_interleaver_pattern = np.zeros(K, dtype=int)
    k = int(0)
    for m in np.arange(Pi_IL_max.size):
        if Pi_IL_max[m] >= (Pi_IL_max.size - K):
            crc_interleaver_pattern[k] = Pi_IL_max[m] - (Pi_IL_max.size - K)
            k = k + 1

    # Get the 3GPP rate matching pattern.
    # Table 5.4.1.1-1: Sub-block interleaver pattern
    SBIP = np.array(
        [0,  1,  2, 4, 3, 5, 6, 7, 8, 16, 9, 17, 10, 18, 11, 19,
         12, 20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31], dtype=int)
    J = np.zeros(N, dtype=int)
    for i_n in np.arange(N, dtype=int):
        i = int(np.floor(32 * i_n // N))
        J[i_n] = SBIP[i] * (N // 32) + np.mod(i_n, N // 32)
    rate_matching_pattern = np.zeros(E, dtype=int)
    if (E >= N):
        for k in np.arange(E):
            rate_matching_pattern[k] = J[np.mod(k, N)]
        mode = "repetition"
    else:
        if (K / E) <= (7 / 16):
            for k in np.arange(E):
                rate_matching_pattern[k] = J[k + N - E]
            mode = "puncturing"
        else:
            for k in np.arange(E):
                rate_matching_pattern[k] = J[k]
            mode = "shortening"
    print("mode: {}".format(mode))
    print("rmp: {}".format(rate_matching_pattern))

    # Get the 3GPP sequence pattern.
    Q_Nmax = np.array(
        [0, 1, 2, 4, 8, 16, 32, 3, 5, 64, 9, 6, 17, 10, 18, 128, 12, 33, 65, 20,
         256, 34, 24, 36, 7, 129, 66, 512, 11, 40, 68, 130, 19, 13, 48, 14, 72,
         257, 21, 132, 35, 258, 26, 513, 80, 37, 25, 22, 136, 260, 264, 38, 514,
         96, 67, 41, 144, 28, 69, 42, 516, 49, 74, 272, 160, 520, 288, 528, 192,
         544, 70, 44, 131, 81, 50, 73, 15, 320, 133, 52, 23, 134, 384, 76, 137,
         82, 56, 27, 97, 39, 259, 84, 138, 145, 261, 29, 43, 98, 515, 88, 140, 30,
         146, 71, 262, 265, 161, 576, 45, 100, 640, 51, 148, 46, 75, 266, 273, 517,
         104, 162, 53, 193, 152, 77, 164, 768, 268, 274, 518, 54, 83, 57, 521, 112,
         135, 78, 289, 194, 85, 276, 522, 58, 168, 139, 99, 86, 60, 280, 89, 290,
         529, 524, 196, 141, 101, 147, 176, 142, 530, 321, 31, 200, 90, 545, 292,
         322, 532, 263, 149, 102, 105, 304, 296, 163, 92, 47, 267, 385, 546, 324,
         208, 386, 150, 153, 165, 106, 55, 328, 536, 577, 548, 113, 154, 79, 269,
         108, 578, 224, 166, 519, 552, 195, 270, 641, 523, 275, 580, 291, 59, 169,
         560, 114, 277, 156, 87, 197, 116, 170, 61, 531, 525, 642, 281, 278, 526,
         177, 293, 388, 91, 584, 769, 198, 172, 120, 201, 336, 62, 282, 143, 103,
         178, 294, 93, 644, 202, 592, 323, 392, 297, 770, 107, 180, 151, 209, 284,
         648, 94, 204, 298, 400, 608, 352, 325, 533, 155, 210, 305, 547, 300, 109,
         184, 534, 537, 115, 167, 225, 326, 306, 772, 157, 656, 329, 110, 117, 212,
         171, 776, 330, 226, 549, 538, 387, 308, 216, 416, 271, 279, 158, 337, 550,
         672, 118, 332, 579, 540, 389, 173, 121, 553, 199, 784, 179, 228, 338, 312,
         704, 390, 174, 554, 581, 393, 283, 122, 448, 353, 561, 203, 63, 340, 394,
         527, 582, 556, 181, 295, 285, 232, 124, 205, 182, 643, 562, 286, 585, 299,
         354, 211, 401, 185, 396, 344, 586, 645, 593, 535, 240, 206, 95, 327, 564,
         800, 402, 356, 307, 301, 417, 213, 568, 832, 588, 186, 646, 404, 227, 896,
         594, 418, 302, 649, 771, 360, 539, 111, 331, 214, 309, 188, 449, 217, 408,
         609, 596, 551, 650, 229, 159, 420, 310, 541, 773, 610, 657, 333, 119, 600,
         339, 218, 368, 652, 230, 391, 313, 450, 542, 334, 233, 555, 774, 175, 123,
         658, 612, 341, 777, 220, 314, 424, 395, 673, 583, 355, 287, 183, 234, 125,
         557, 660, 616, 342, 316, 241, 778, 563, 345, 452, 397, 403, 207, 674, 558,
         785, 432, 357, 187, 236, 664, 624, 587, 780, 705, 126, 242, 565, 398, 346,
         456, 358, 405, 303, 569, 244, 595, 189, 566, 676, 361, 706, 589, 215, 786,
         647, 348, 419, 406, 464, 680, 801, 362, 590, 409, 570, 788, 597, 572, 219,
         311, 708, 598, 601, 651, 421, 792, 802, 611, 602, 410, 231, 688, 653, 248,
         369, 190, 364, 654, 659, 335, 480, 315, 221, 370, 613, 422, 425, 451, 614,
         543, 235, 412, 343, 372, 775, 317, 222, 426, 453, 237, 559, 833, 804, 712,
         834, 661, 808, 779, 617, 604, 433, 720, 816, 836, 347, 897, 243, 662, 454,
         318, 675, 618, 898, 781, 376, 428, 665, 736, 567, 840, 625, 238, 359, 457,
         399, 787, 591, 678, 434, 677, 349, 245, 458, 666, 620, 363, 127, 191, 782,
         407, 436, 626, 571, 465, 681, 246, 707, 350, 599, 668, 790, 460, 249, 682,
         573, 411, 803, 789, 709, 365, 440, 628, 689, 374, 423, 466, 793, 250, 371,
         481, 574, 413, 603, 366, 468, 655, 900, 805, 615, 684, 710, 429, 794, 252,
         373, 605, 848, 690, 713, 632, 482, 806, 427, 904, 414, 223, 663, 692, 835,
         619, 472, 455, 796, 809, 714, 721, 837, 716, 864, 810, 606, 912, 722, 696,
         377, 435, 817, 319, 621, 812, 484, 430, 838, 667, 488, 239, 378, 459, 622,
         627, 437, 380, 818, 461, 496, 669, 679, 724, 841, 629, 351, 467, 438, 737,
         251, 462, 442, 441, 469, 247, 683, 842, 738, 899, 670, 783, 849, 820, 728,
         928, 791, 367, 901, 630, 685, 844, 633, 711, 253, 691, 824, 902, 686, 740,
         850, 375, 444, 470, 483, 415, 485, 905, 795, 473, 634, 744, 852, 960, 865,
         693, 797, 906, 715, 807, 474, 636, 694, 254, 717, 575, 913, 798, 811, 379,
         697, 431, 607, 489, 866, 723, 486, 908, 718, 813, 476, 856, 839, 725, 698,
         914, 752, 868, 819, 814, 439, 929, 490, 623, 671, 739, 916, 463, 843, 381,
         497, 930, 821, 726, 961, 872, 492, 631, 729, 700, 443, 741, 845, 920, 382,
         822, 851, 730, 498, 880, 742, 445, 471, 635, 932, 687, 903, 825, 500, 846,
         745, 826, 732, 446, 962, 936, 475, 853, 867, 637, 907, 487, 695, 746, 828,
         753, 854, 857, 504, 799, 255, 964, 909, 719, 477, 915, 638, 748, 944, 869,
         491, 699, 754, 858, 478, 968, 383, 910, 815, 976, 870, 917, 727, 493, 873,
         701, 931, 756, 860, 499, 731, 823, 922, 874, 918, 502, 933, 743, 760, 881,
         494, 702, 921, 501, 876, 847, 992, 447, 733, 827, 934, 882, 937, 963, 747,
         505, 855, 924, 734, 829, 965, 938, 884, 506, 749, 945, 966, 755, 859, 940,
         830, 911, 871, 639, 888, 479, 946, 750, 969, 508, 861, 757, 970, 919, 875,
         862, 758, 948, 977, 923, 972, 761, 877, 952, 495, 703, 935, 978, 883, 762,
         503, 925, 878, 735, 993, 885, 939, 994, 980, 926, 764, 941, 967, 886, 831,
         947, 507, 889, 984, 751, 942, 996, 971, 890, 509, 949, 973, 1000, 892, 950,
         863, 759, 1008, 510, 979, 953, 763, 974, 954, 879, 981, 982, 927, 995, 765,
         956, 887, 985, 997, 986, 943, 891, 998, 766, 511, 988, 1001, 951, 1002, 893,
         975, 894, 1009, 955, 1004, 1010, 957, 983, 958, 987, 1012, 999, 1016, 767,
         989, 1003, 990, 1005, 959, 1011, 1013, 895, 1006, 1014, 1017, 1018, 991,
         1020, 1007, 1015, 1019, 1021, 1022, 1023], dtype=int)
    Q_N = Q_Nmax[Q_Nmax < N]
    print("Q_N: {}".format(Q_N))

    # Get the 3GPP information bit pattern.
    info_bit_pattern = np.zeros(N, dtype=int)
    # Find the set difference of two arrays. Return the unique values in ar1 that are not in ar2.
    Q_Ftmp_N = np.setdiff1d(np.arange(N), rate_matching_pattern)
    if mode == "puncturing":
        if (E >= 3 * N // 4):
            Q_Ftmp_N = np.concatenate(
                (Q_Ftmp_N, np.arange(np.ceil(3 * N / 4 - E / 2), dtype=int)))
        else:
            Q_Ftmp_N = np.concatenate(
                (Q_Ftmp_N, np.arange(np.ceil(9 * N / 16 - E / 4), dtype=int)))
    Q_Itmp_N = np.setdiff1d(Q_N, Q_Ftmp_N, assume_unique=True)
    Q_I_N = Q_Itmp_N[-K:]
    info_bit_pattern[Q_I_N] = 1

    # Generate CRC matrix
    G_P = np.zeros((K, P), dtype=int)
    G_P[-1, :] = crc_polynomial_pattern[1:]
    for k in np.arange(K - 2, -1, -1, dtype=int):
        G_P[k, :] = np.bitwise_xor(
            np.concatenate((G_P[k+1, 1:], np.zeros(1, dtype=int))), G_P[k+1, 0] * crc_polynomial_pattern[1:])

    # Info bits
    # a = np.random.randint(2, size=A)
    a = np.array([1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1], dtype=int)
    print("a:", a)

    # RNTI = np.random.randint(2, size=16)
    RNTI = np.array(
        [0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1], dtype=int)
    print("RNTI:", RNTI)

    # Generate the scrambled CRC bits.
    crc_bits = np.mod(np.matmul(np.concatenate(
        (np.ones(P, dtype=int), a)), G_P), 2)
    scrambled_crc_bits = np.bitwise_xor(
        crc_bits, np.concatenate((np.zeros(P-16, dtype=int), RNTI)))

    # Append the scrambled CRC bits to the information bits.
    b = np.concatenate((a, scrambled_crc_bits))
    print("b:", b)

    # Interleave the information and CRC bits.
    c = b[crc_interleaver_pattern]
    print("c:", c)

    # Position the information and CRC bits within the input to the polar
    # encoder kernal.
    u = np.zeros(N, dtype=int)
    u[info_bit_pattern == 1] = c
    print("u:", u)

    # Perform the polar encoder kernal operation.
    G_2 = np.array([[1, 0], [1, 1]])
    G_N = G_2
    for i in range(n - 1):
        G_N = np.kron(G_2, G_N)
    d = np.mod(np.matmul(u, G_N), 2)
    print("d:", d)

    # Extract the encoded bits from the output of the polar encoder kernal.
    f = d[rate_matching_pattern]
    print("f:", f)

    # QPSK modulation
    tx = np.sqrt(1/2) * ((1-2*f[0::2]) + (1-2*f[1::2])*1j)

    # AWGN
    N0 = 1/(10 ** (EsN0/10))
    rx = tx + np.sqrt(N0/2) * (np.random.randn(tx.size) +
                               np.random.randn(tx.size)*1j)

    # QPSK demodulation
    f_tilde = np.zeros_like(f)
    f_tilde[0::2] = 4 * np.sqrt(1/2) * rx.real / N0
    f_tilde[1::2] = 4 * np.sqrt(1/2) * rx.imag / N0

    # Rate matching
    if mode == "repetition":
        d_tilde = np.zeros(N)
        for i in np.arange(E):
            d_tilde[rate_matching_pattern[i]] = \
                d_tilde[rate_matching_pattern[i]] + f_tilde[i]
    else:
        if mode == "puncturing":
            d_tilde = np.zeros(N)
        else:
            # shortening
            d_tilde = np.full(N, np.inf)
        d_tilde[rate_matching_pattern] = f_tilde
    print("d_tilde:", d_tilde)

    # Successive cancellation decoder
    LLRS = np.full((N, n + 1), np.nan)
    BITS = np.full((N, n + 1), np.nan)
    LLRS[:, 0] = d_tilde
    for i in np.arange(N):
        # Evaluate tree of LLRs for root index i
        act_llr_level = \
            decoder_utils.active_llr_level(utils.bit_reversed(i, n), n)
        for s in np.arange(n - act_llr_level, n):
            block_size = int(2 ** (n - s - 1))
            for j in np.arange(i, i + block_size):
                if j < block_size:  # upper branch
                    top_llr = LLRS[j, s]
                    btm_llr = LLRS[j + block_size, s]
                    LLRS[j, s + 1] = \
                        decoder_utils.upper_llr_approx(top_llr, btm_llr)
                else:  # lower branch
                    btm_llr = LLRS[j, s]
                    top_llr = LLRS[j - block_size, s]
                    top_bit = BITS[j - block_size, s + 1]
                    LLRS[j, s + 1] = \
                        decoder_utils.lower_llr(btm_llr, top_llr, top_bit)
        # Make hard decision at output
        if info_bit_pattern[i] == 0:
            BITS[i, n] = 0
        else:
            BITS[i, n] = \
                decoder_utils.hard_decision(LLRS[i, n])


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.n_info, args.n_rm, args.esn0)
