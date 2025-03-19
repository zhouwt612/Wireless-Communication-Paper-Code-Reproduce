import numpy as np
from scipy.optimize import linear_sum_assignment
import itertools
# ✅ 원래 비용 행렬 (음수 포함)

def hungarian(cost_matrix):
# ✅ 음수를 없애기 위해 최소값을 빼줌
    min_value = cost_matrix.min()
    cost_matrix += abs(min_value)  # 최소값의 절댓값을 더해줌

    # ✅ 헝가리안 알고리즘 적용
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # ✅ 최적의 작업 할당 결과 출력
    #print("Optimal Assignment:")
    #for i, j in zip(row_ind, col_ind):
    #    print(f"Worker {i} → Task {j}")

    # ✅ 최적 비용 출력
    optimal_cost = cost_matrix[row_ind, col_ind].sum() - (abs(min_value) * len(cost_matrix))
    #print(type(col_ind))
    #print(f"Total Minimum Cost: {optimal_cost}")  # 원래 비용 계산으로 복원

    return row_ind, col_ind

def exhaustive(cost_matrix):
    n = cost_matrix.shape[0]  # 행렬 크기 (n x n)

    # ✅ 모든 가능한 작업 할당 (순열)
    assignments = itertools.permutations(range(n))

    # ✅ 최소 비용 찾기
    min_cost = float('inf')
    best_assignment = None

    for assignment in assignments:
        cost = sum(cost_matrix[i, assignment[i]] for i in range(n))
        if cost < min_cost:
            min_cost = cost
            best_assignment = assignment

    return best_assignment

def fixed_point(CPU_param, Z, gamma_E, h_k, B_U, N0):
    iter_num = 100
    P = 1
    for i in range(iter_num):
        value = (gamma_E - CPU_param * Z) * B_U / Z * np.log2(1+P * abs(h_k)**2 / (B_U * N0))
        #print(value)
        P_new = value
        if abs(P_new - P) < 0.0001:
            P = P_new
            break
        P = P_new
    #print('fixed point :',P)
    return P