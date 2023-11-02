#!/usr/bin/env python

"""
A polar decoder class. Currently only Successive Cancellation Decoder (SCD) is supported.
"""

import numpy as np
from polarcodes.utils import *
from polarcodes.PSCD import PSCD


class Decode:
    def __init__(self, myPC, decoder_name='pscd'):
        """
        Parameters
        ----------
        myPC: `PolarCode`
            a polar code object created using the :class:`PolarCode` class
        decoder_name: string
            name of decoder to use (default is 'scd')
        """

        self.myPC = myPC
        self.x_noisy = np.array([])

        # select decoding algorithm
        if decoder_name == 'pscd':
            scd = PSCD(myPC)
            self.x_noisy = scd.decode()
            self.myPC.message_received = self.noisy_message(
                self.x_noisy, False)
        elif decoder_name == 'systematic_pscd':
            scd = PSCD(myPC)
            self.x_noisy = scd.decode()
            self.myPC.message_received = self.noisy_message(self.x_noisy, True)

    def noisy_message(self, x_noisy, systematic_flag):
        if systematic_flag:
            x_noisy = self.systematic_decode(x_noisy)
        return x_noisy[self.myPC.frozen_lookup == 1]

    def systematic_decode(self, x_noisy):
        x = np.array([x_noisy], dtype=int)
        return np.transpose(np.mod(np.dot(self.myPC.T, x.T), 2))[0]
