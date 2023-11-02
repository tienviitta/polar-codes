import numpy as np
from polarcodes.utils import *
from polarcodes.decoder_utils import *


class PSCD:
    def __init__(self, myPC):
        self.myPC = myPC
        self.L = np.full((self.myPC.N, self.myPC.n + 1),
                         np.nan, dtype=np.float64)
        self.B = np.full((self.myPC.N, self.myPC.n + 1), np.nan)
        self.L[:, 0] = self.myPC.likelihoods

    def decode(self):
        """
        Successive Cancellation Decoder. The decoded message is set to 
        ``message_received`` in ``myPC``. The decoder will use the frozen set as 
        defined by ``frozen`` in ``myPC``. 
        Depends on `update_llrs` and `update_bits`.

        Parameters
        ----------
        y: ndarray<float>
            a vector of likelihoods at the channel output

        -------------
        **References:**

        We first propose a new decoding algorithm on a new decoder graph, where 
        the various stages of the graph are permuted. We then observe that, even
        though the usage of the permuted graph doesn't affect the encoder, it 
        can significantly affect the decoding performance of a given polar code.
        The new permuted successive cancellation decoder (PSCD) typically 
        exhibits a performance degradation, since the polar code is optimized 
        for the standard SCD.

        *  Vangala, H., Viterbo, & Yi Hong. (2014). Permuted successive 
        cancellation decoder for polar codes. 2014 International Symposium on 
        Information Theory and Its Applications, 438-442. IEICE.

        """

        # decode bits in natural order
        for l in [bit_reversed(i, self.myPC.n) for i in range(self.myPC.N)]:
            # evaluate tree of LLRs for root index i
            self.update_llrs(l)

            # make hard decision at output
            if l in self.myPC.frozen:
                self.B[l, self.myPC.n] = 0
            else:
                self.B[l, self.myPC.n] = hard_decision(self.L[l, self.myPC.n])

            # propagate the hard decision just made
            self.update_bits(l)
        return self.B[:, self.myPC.n].astype(int)

    def update_llrs(self, l):
        for s in range(self.myPC.n - active_llr_level(l, self.myPC.n), self.myPC.n):
            block_size = int(2 ** (s + 1))
            branch_size = int(block_size / 2)
            for j in range(l, self.myPC.N, block_size):
                if j % block_size < branch_size:  # upper branch
                    top_llr = self.L[j, s]
                    btm_llr = self.L[j + branch_size, s]
                    # self.L[j, s + 1] = upper_llr(top_llr, btm_llr)
                    self.L[j, s + 1] = upper_llr_approx(top_llr, btm_llr)
                else:  # lower branch
                    btm_llr = self.L[j, s]
                    top_llr = self.L[j - branch_size, s]
                    top_bit = self.B[j - branch_size, s + 1]
                    self.L[j, s + 1] = lower_llr(btm_llr, top_llr, top_bit)

    def update_bits(self, l):
        if l < int(self.myPC.N / 2):
            return

        for s in range(self.myPC.n, self.myPC.n - active_bit_level(l, self.myPC.n), -1):
            block_size = int(2 ** s)
            branch_size = int(block_size / 2)
            for j in range(l, -1, -block_size):
                if j % block_size >= branch_size:  # lower branch
                    self.B[j - branch_size, s - 1] = \
                        int(self.B[j, s]) ^ int(self.B[j - branch_size, s])
                    self.B[j, s - 1] = self.B[j, s]
