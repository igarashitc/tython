import numpy as np

def grad(qq, x):
    nx = qq.shape[0]

    dqq = np.zeros(nx)

    for i in range(1, nx-1):
        dqq[i] = (qq[i+1] - qq[i-1])/(x[i+1] - x[i-1])

    return dqq

