import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # TkAgg 백엔드 사용
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
import os
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import Communication
import FL
import Dataset


# ✅ 환경 변수 설정
@dataclass
class CommunicationParams:
    K: int = 5
    N: int = 10
    Nb: int = 1
    Nu: int = 1
    sigma: float = 1
    P: float = 10
    P_max: float = 0.01
    alpha: float = 2
    P_B: float = 1
    N0: float = field(default=10 ** (-50 / 10) * 10 ** (-3))
    B_D: float = field(default=20 * 10 ** 6)
    B_U: float = field(default=1 * 10 ** 6)
    gamma_T: float = field(default=500 * 10 ** (-3))
    gamma_E: float = field(default=3 * 10 ** (-3))
    energy_consume: float = field(default=10 ** (-27))
    CPU_freq: float = field(default=10 ** 9)
    CPU_cycle: int = 40
    Z: int = 5 * 10 ** 4

    # m_dB 값을 필드로 선언
    m_dB: np.ndarray = field(default_factory=lambda: np.array([5.782, 7.083, -0.983, 0.023, -4.401, -4.312]))

    # dB → 선형 변환 (default_factory 사용)
    m: np.ndarray = field(init=False)  # 자동 초기화 필드

    def __post_init__(self):
        """dB 값을 선형 값으로 변환 (m_dB → m)"""
        self.m = 10 ** (self.m_dB / 10)


params = CommunicationParams()
subset_size = 60000 // params.K
batch = 100
batch_set = [12, 10, 8, 4, 2]
num_epochs = 100
learning_rate = 0.01

device = "cuda" if torch.cuda.is_available() else "cpu"
train_loaders, test_loader = Dataset.MNIST_dataloader(batch, batch_set, subset_size, params)


def train_federated_learning(h_k_mode):
    models = [FL.FNN().to(device) for _ in range(params.K)]
    optimizers = [optim.Adam(models[k].parameters(), lr=learning_rate) for k in range(params.K)]
    criterions = [nn.CrossEntropyLoss() for _ in range(params.K)]
    central_model = FL.FNN().to(device)
    accuracies = []

    for epoch in range(num_epochs):
        print(f"\n=== Epoch {epoch + 1}/{num_epochs} ({h_k_mode}) ===")

        h_k = (np.random.normal(0, params.sigma, (params.K, params.N))
               + 1j * np.random.normal(0, params.sigma, (params.K, params.N))) / np.sqrt(2)



        if h_k_mode == "allocated":
            P_allocated, h_k_allocated, x = Communication.model(h_k, batch_set, params)
        else:
            h_k_random_allocation = (np.random.normal(0, params.sigma, params.K)
                                     + 1j * np.random.normal(0, params.sigma, params.K)) / np.sqrt(2)
            P_allocated = params.P_max * np.ones((params.K, params.N))
            h_k_allocated = h_k_random_allocation
            x = np.zeros((params.K), dtype=int)
        num_joining_data = 0
        total_loss = 0
        total_grad_sum = None

        num_joining_data = 0
        total_loss = 0
        total_grad_sum = None
        at_least_one_joined = False  # ✅ 최소 한 명이 참여했는지 체크

        for k in range(params.K):
            loss, join, total_grad_sum = FL.local_NN(k, h_k_allocated, P_allocated, x, batch_set, params, train_loaders,
                                                     models,
                                                     optimizers, criterions, device, total_grad_sum)
            total_loss += loss
            num_joining_data += join * batch_set[k]

            if join > 0:
                at_least_one_joined = True  # ✅ 최소 한 명이라도 참여

        # ✅ 모든 유저가 참여하지 못했을 경우, 마지막 유저라도 강제로 참여
        if not at_least_one_joined:
            print("⚠️ No users participated! Forcing last user to contribute...")
            loss, join, total_grad_sum = FL.local_NN(params.K - 1, h_k_allocated, P_allocated, x, batch_set, params,
                                                     train_loaders, models, optimizers, criterions, device,
                                                     total_grad_sum)
            total_loss += loss
            num_joining_data += batch_set[params.K - 1]



        FL.central_NN(num_joining_data, optimizers, central_model, models, learning_rate, total_grad_sum, total_loss,
                      params)

        # 중앙 모델 정확도 측정
        central_model.eval()
        correct = 0
        total = 0
        with torch.no_grad():
            for images, labels in test_loader:
                images, labels = images.to(device), labels.to(device)
                outputs = central_model(images)
                _, predicted = torch.max(outputs, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
        accuracy = 100 * correct / total
        accuracies.append(accuracy)
        print(f"{h_k_mode.capitalize()} Model Test Accuracy: {accuracy:.2f}%")

    return accuracies


# ✅ 두 가지 방법 비교
accuracies_allocated = train_federated_learning("allocated")
accuracies_random = train_federated_learning("random_allocation")

# ✅ 수렴 속도 비교 그래프 출력
plt.figure(figsize=(10, 5))
plt.plot(range(1, num_epochs + 1), accuracies_allocated, label="h_k_allocated", marker='o')
plt.plot(range(1, num_epochs + 1), accuracies_random, label="h_k_random_allocation", marker='s', linestyle='dashed')
plt.xlabel("Epochs")
plt.ylabel("Accuracy (%)")
plt.title("FL Convergence Speed: h_k_allocated vs. h_k_random_allocation")
plt.legend()
plt.grid(True)
plt.show()
