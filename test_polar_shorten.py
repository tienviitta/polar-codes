import numpy as np
from polarcodes import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-M", "--block_length", type=int)
parser.add_argument("-K", "--code_dim", type=int)


def main(M, K):
    # initialise shortened polar code
    shorten_params = ('shorten', 'brs', None, None, False)
    myPC = PolarCode(M, K, shorten_params)

    # shortened construction
    design_SNR = 5.0
    Shorten(myPC, design_SNR)
    print(myPC, "\n\n")

    # set message
    my_message = np.random.randint(2, size=myPC.K)
    myPC.set_message(my_message)
    print("The message is:", my_message)

    # encode message
    Encode(myPC)
    print("The coded message is:", myPC.get_codeword())

    # transmit the codeword
    AWGN(myPC, design_SNR)
    print("The log-likelihoods are:", myPC.likelihoods)

    # decode the received codeword
    Decode(myPC)
    print("The decoded message is:", myPC.message_received)

    # check results
    print("PASS:", np.all(myPC.message == myPC.message_received))


if __name__ == "__main__":
    args = parser.parse_args()
    main(args.block_length, args.code_dim)
