"""
The "chk" version: upper_llr(top_llr, btm_llr)

Design SNR: 0.0 dB

test_polar_nopunct.py -M 8 -K 5 

========== Polar Code ==========
N: 8
M: 8
K: 5
Mothercode Construction: bb
Ordered Bits (least reliable to most reliable): [0 4 2 1 6 5 3 7]
Frozen Bits: [2 4 0]
Puncturing Flag: False
Puncturing Parameters: {punct_type: 
                        punct_algorithm: 
                        punct_set: []
                        source_set: []
                        update_frozen_flag: None}
 


The message is: [1 0 0 1 0]
The coded message is: [0 1 1 0 1 0 1 0]
The log-likelihoods are: [ 0.16082923 -2.05397914 -1.45657378  0.55757278  0.02308327  1.96730347 -1.74552562  1.37325182]

Natural order?:
l = [0, 4, 2, 6, 1, 5, 3, 7]

l = 0:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.3414,     nan,     nan],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0174, -0.0073,     nan],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -0.8923,     nan,     nan],
       [ 1.3733,     nan,     nan,     nan]])
array([[nan, nan, nan,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan]])

l = 4:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.3414,     nan,     nan],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -0.8923,     nan,     nan],
       [ 1.3733,     nan,     nan,     nan]])
array([[nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan]])

l = 2:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.3414, -0.4655,  0.1887],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -0.8923, -0.8749,     nan],
       [ 1.3733,     nan,     nan,     nan]])
array([[nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan,  0.],
       [nan, nan, nan, nan],
       [nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan]])

l = 6:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.3414, -0.4655,  0.1887],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -0.8923, -0.8749, -1.3404],
       [ 1.3733,     nan,     nan,     nan]])
array([[nan,  1.,  0.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  0.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

l = 1:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 , -2.2148, -1.4308, -0.8986],
       [-1.4566, -0.3414, -0.4655,  0.1887],
       [ 0.5576,  2.0141,     nan,     nan],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,  1.9442,  1.6813,     nan],
       [-1.7455, -0.8923, -0.8749, -1.3404],
       [ 1.3733,  3.1188,     nan,     nan]])
array([[nan,  1.,  0.,  0.],
       [nan, nan, nan,  1.],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  0.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

l = 5:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 , -2.2148, -1.4308, -0.8986],
       [-1.4566, -0.3414, -0.4655,  0.1887],
       [ 0.5576,  2.0141,     nan,     nan],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,  1.9442,  1.6813,  3.1121],
       [-1.7455, -0.8923, -0.8749, -1.3404],
       [ 1.3733,  3.1188,     nan,     nan]])
array([[nan,  1.,  0.,  0.],
       [nan, nan,  1.,  1.],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  0.,  0.],
       [nan, nan,  0.,  0.],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

l = 3:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 , -2.2148, -1.4308, -0.8986],
       [-1.4566, -0.3414, -0.4655,  0.1887],
       [ 0.5576,  2.0141,  4.229 ,  3.8684],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,  1.9442,  1.6813,  3.1121],
       [-1.7455, -0.8923, -0.8749, -1.3404],
       [ 1.3733,  3.1188,  5.063 ,     nan]])
array([[nan,  1.,  0.,  0.],
       [nan, nan,  1.,  1.],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan,  0.],
       [nan,  1.,  0.,  0.],
       [nan, nan,  0.,  0.],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

l = 7:
array([[ 0.1608, -0.1242,  0.021 , -0.0001],
       [-2.054 , -2.2148, -1.4308, -0.8986],
       [-1.4566, -0.3414, -0.4655,  0.1887],
       [ 0.5576,  2.0141,  4.229 ,  3.8684],
       [ 0.0231,  0.0174, -0.0073,  0.0137],
       [ 1.9673,  1.9442,  1.6813,  3.1121],
       [-1.7455, -0.8923, -0.8749, -1.3404],
       [ 1.3733,  3.1188,  5.063 ,  9.292 ]])
array([[0., 1., 0., 0.],
       [1., 1., 1., 1.],
       [1., 1., 1., 0.],
       [0., 0., 0., 0.],
       [1., 1., 0., 0.],
       [0., 0., 0., 0.],
       [1., 1., 1., 1.],
       [0., 0., 0., 0.]])

self.x_noisy: array([0, 1, 0, 0, 0, 0, 1, 0])

myPC.message_received: array([1, 0, 0, 1, 0])

The decoded message is: [1 0 0 1 0]
PASS: True

"""

"""
The "chk" approximation: upper_llr_approx(top_llr, btm_llr)

Design SNR: 0.0 dB

test_polar_nopunct.py -M 8 -K 5 

========== Polar Code ==========
N: 8
M: 8
K: 5
Mothercode Construction: bb
Ordered Bits (least reliable to most reliable): [0 4 2 1 6 5 3 7]
Frozen Bits: [2 4 0]
Puncturing Flag: False
Puncturing Parameters: {punct_type: 
                        punct_algorithm: 
                        punct_set: []
                        source_set: []
                        update_frozen_flag: None}
 


The message is: [1 0 0 1 0]
The coded message is: [0 1 1 0 1 0 1 0]
The log-likelihoods are: [ 0.1608 -2.054  -1.4566  0.5576  0.0231  1.9673 -1.7455  1.3733]
Order: l = [0, 4, 2, 6, 1, 5, 3, 7]

j = 0, l = 0:
array([[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.5576,     nan,     nan],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0231, -0.0231,     nan],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -1.3733,     nan,     nan],
       [ 1.3733,     nan,     nan,     nan]])
BIT = [[nan, nan, nan,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan]])

j = 1, l = 4:
array([[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.5576,     nan,     nan],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -1.3733,     nan,     nan],
       [ 1.3733,     nan,     nan,     nan]])
BIT = [[nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan]])

j = 2, l = 2:
array([[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.5576, -0.7184,  0.7184],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -1.3733, -1.3502,     nan],
       [ 1.3733,     nan,     nan,     nan]])
BIT = [[nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan,  0.],
       [nan, nan, nan, nan],
       [nan, nan,  0.,  0.],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan],
       [nan, nan, nan, nan]])

j = 3, l = 6:
array([[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 ,     nan,     nan,     nan],
       [-1.4566, -0.5576, -0.7184,  0.7184],
       [ 0.5576,     nan,     nan,     nan],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,     nan,     nan,     nan],
       [-1.7455, -1.3733, -1.3502, -2.0686],
       [ 1.3733,     nan,     nan,     nan]])
BIT = [[nan,  1.,  0.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  0.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

j = 4, l = 1:
array([[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 , -2.2148, -2.0141, -1.9442],
       [-1.4566, -0.5576, -0.7184,  0.7184],
       [ 0.5576,  2.0141,     nan,     nan],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,  1.9442,  1.9442,     nan],
       [-1.7455, -1.3733, -1.3502, -2.0686],
       [ 1.3733,  3.1188,     nan,     nan]])
BIT = [[nan,  1.,  0.,  0.],
       [nan, nan, nan,  1.],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  0.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

j = 5, l = 5:
array([[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 , -2.2148, -2.0141, -1.9442],
       [-1.4566, -0.5576, -0.7184,  0.7184],
       [ 0.5576,  2.0141,     nan,     nan],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,  1.9442,  1.9442,  3.9584],
       [-1.7455, -1.3733, -1.3502, -2.0686],
       [ 1.3733,  3.1188,     nan,     nan]])
BIT = [[nan,  1.,  0.,  0.],
       [nan, nan,  1.,  1.],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan, nan],
       [nan,  1.,  0.,  0.],
       [nan, nan,  0.,  0.],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

l = 3:
LLR = [[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 , -2.2148, -2.0141, -1.9442],
       [-1.4566, -0.5576, -0.7184,  0.7184],
       [ 0.5576,  2.0141,  4.229 ,  4.229 ],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,  1.9442,  1.9442,  3.9584],
       [-1.7455, -1.3733, -1.3502, -2.0686],
       [ 1.3733,  3.1188,  5.063 ,     nan]])
BIT = [[nan,  1.,  0.,  0.],
       [nan, nan,  1.,  1.],
       [nan,  1.,  1.,  0.],
       [nan, nan, nan,  0.],
       [nan,  1.,  0.,  0.],
       [nan, nan,  0.,  0.],
       [nan,  1.,  1.,  1.],
       [nan, nan, nan, nan]])

l = 7:
LLR = [[ 0.1608, -0.1608,  0.1608, -0.0231],
       [-2.054 , -2.2148, -2.0141, -1.9442],
       [-1.4566, -0.5576, -0.7184,  0.7184],
       [ 0.5576,  2.0141,  4.229 ,  4.229 ],
       [ 0.0231,  0.0231, -0.0231,  0.1377],
       [ 1.9673,  1.9442,  1.9442,  3.9584],
       [-1.7455, -1.3733, -1.3502, -2.0686],
       [ 1.3733,  3.1188,  5.063 ,  9.292 ]])
BIT = [[0., 1., 0., 0.],
       [1., 1., 1., 1.],
       [1., 1., 1., 0.],
       [0., 0., 0., 0.],
       [1., 1., 0., 0.],
       [0., 0., 0., 0.],
       [1., 1., 1., 1.],
       [0., 0., 0., 0.]])

self.x_noisy: array([0, 1, 0, 0, 0, 0, 1, 0])

myPC.message_received: array([1, 0, 0, 1, 0])

The decoded message is: [1 0 0 1 0]
PASS: True

"""

...
