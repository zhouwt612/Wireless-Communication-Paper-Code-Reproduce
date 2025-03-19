import numpy as np
import Functions

##### modeling communication envionrment #####
##### outputs are power allocated to optimal, CH and allocated sub-carriers #####
def model(h_k, batch_set, params):
    K = h_k.shape[0]
    N = h_k.shape[1]

    SINR_kn = np.zeros((K, N))
    q_kn = np.zeros((K, N))
    # define error rate of each user in each sub-carrier
    for k in range(K):
        for n in range(N):
            SNR = params.P / params.sigma
            SINR_kn[k][n] = abs(h_k[k][n]) ** 2 * SNR
            q_kn[k][n] = 1 - np.exp(-params.m[0] / SINR_kn[k][n])

    # optimal power allocation
    P_i_opt = np.zeros((K, N))
    CPU_param = params.energy_consume * params.CPU_cycle * params.CPU_freq ** 2
    for k in range(K):
        for n in range(N):
            P_i_gammaE = Functions.fixed_point(CPU_param, params.Z, params.gamma_E, h_k[k][n], params.B_U, params.N0)
            P_i_opt[k][n] = min(params.P_max, P_i_gammaE)

    # define data rate, delay, energy of each user in each sub-carrier
    c_U = np.zeros((K, N))  # uplink data rate
    c_D = np.zeros((K, N))  # downlink data rate
    l_U = np.zeros((K, N))  # uplink delay
    l_D = np.zeros((K, N))  # downlink delay
    e = np.zeros((K, N))
    for k in range(K):
        for n in range(N):
            c_U[k][n] = params.B_U * np.log2(1 + np.abs(h_k[k][n]) ** 2 * P_i_opt[k][n] / (params.B_U * params.N0))
            c_D[k][n] = params.B_D * np.log2(1 + np.abs(h_k[k][n]) ** 2 * params.P_B / (params.B_D * params.N0))
            l_U[k][n] = params.Z / c_U[k][n]
            l_D[k][n] = params.Z / c_D[k][n]
            e[k][n] = CPU_param * params.Z + P_i_opt[k][n] * l_U[k][n]
    #print(1/c_U[0][0])
    # define Hungarian matrix
    psi = np.zeros((K, N))
    for k in range(K):
        for n in range(N):
            if l_U[k][n] + l_D[k][n] <= params.gamma_T and e[k][n] <= params.gamma_E:
                psi[k][n] = batch_set[k] * (q_kn[k][n] - 1)
            else:
                psi[k][n] = 0

    _, x = Functions.hungarian(psi)
    h_k_allocated = np.zeros(K, dtype=complex)
    for k in range(K):
        h_k_allocated[k] = h_k[k][x[k]]

    return P_i_opt, h_k_allocated, x