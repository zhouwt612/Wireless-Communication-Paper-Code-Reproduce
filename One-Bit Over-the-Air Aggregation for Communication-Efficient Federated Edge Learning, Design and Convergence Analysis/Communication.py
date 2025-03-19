import numpy as np

import Functions

##### modeling communication envionrment #####
##### outputs are power allocated to optimal, CH and allocated sub-carriers #####
def power_allocation(h_k, params):
    Ei_g_th, _ = Functions.exponential_integratation(params.g_th)

    rho_0 = params.P0 / (params.M * Ei_g_th) # scaling factor
    rho_1 = 1

    p_k = np.zeros((params.K, params.M, 1)) # power allocation matrix
    for k in range(params.K):
        for m in range(params.M):
            h = abs(h_k[k][m])**2
            if h >= params.g_th:
                p_k[k][m] = np.sqrt(rho_0)
            else:
                p_k[k][m] = 0

    return p_k

def majority_vote_decoder(g_tilde):
    length = len(g_tilde)

    v = np.zeros((length,1))

    for l in range(length):
        if g_tilde[l] >= 0:
            v[l] = 1
        else:
            v[l] = -1

    return v