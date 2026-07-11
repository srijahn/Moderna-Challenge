import numpy as np

num_variables = 5

Q = np.zeros((num_variables, num_variables))

Q[0,0] = -1
Q[1,1] = -1
Q[2,2] = -1
Q[3,3] = -1
Q[4,4] = -1

Q[0,1] = 5
Q[1,0] = 5

Q[2,3] = 5
Q[3,2] = 5

print(Q)