import numpy as np
from polarcodes import *
import argparse

np.set_printoptions(precision=4, suppress=True)

parser = argparse.ArgumentParser()
parser.add_argument("-M", "--block_length", type=int)
parser.add_argument("-K", "--code_dim", type=int)


def main(M, K):
    # Random seed
    np.random.seed(123456)

    # initialise polar code
    myPC = PolarCode(M, K)
    myPC.construction_type = 'bb'

    # mothercode construction
    design_SNR = 0.0
    Construct(myPC, design_SNR)
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
