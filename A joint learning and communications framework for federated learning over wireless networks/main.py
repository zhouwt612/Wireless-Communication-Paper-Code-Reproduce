import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
import numpy as np
import random
from torch.utils.data import DataLoader, SubsetRandomSampler
from dataclasses import dataclass, field

import Functions
import Communication
import FL
import Dataset

# dB | Waterfall threshold ref. by On the Frame Error Rate of Transmission Schemes
# on Quasi-Static Fading Channels
# BPSK
# m[0] : uncoded block length = 256
# m[1] : uncoded block length = 1024
# m[2] : convolutioanl coding block length = 256
# m[3] : convolutioanl coding block length = 1024
# m[4] : turbo coding block length = 256
# m[5] : turbo coding block length = 1024

# ✅ 환경 변수 설정
@dataclass
# ✅ Communication parameters
class CommunicationParams:
    K: int = 5  # 사용자 수
    N: int = 10  # 서브 채널 수 (K와 동일)
    Nb: int = 1  # 기지국 안테나 수
    Nu: int = 1  # UE 안테나 수
    sigma: float = 1  # 노이즈 파워
    P: float = 10  # TX 파워
    P_max: float = 0.01  # 최대 전송 전력
    alpha: float = 2  # 페이딩 지수
    P_B: float = 1  # 기지국 전력
    N0: float = field(default=10 ** (-174 / 10) * 10 ** (-3))  # 노이즈 밀도
    B_D: float = field(default=20 * 10 ** 6)  # 다운링크 대역폭 (20MHz)
    B_U: float = field(default=1 * 10 ** 6)  # 업링크 대역폭 (1MHz)
    gamma_T: float = field(default=500 * 10 ** (-3))  # 시간 제약 (0.5s)
    gamma_E: float = field(default=3 * 10 ** (-3))  # 에너지 제약 (0.003J)
    energy_consume: float = field(default=10 ** (-27))  # 에너지 소비율
    CPU_freq: float = field(default=10 ** 9)  # CPU 클럭 속도 (1GHz)
    CPU_cycle: int = 40  # CPU 사이클 수
    Z: int = 5*10**4  # 데이터 크기 (예제 값)

    # m_dB 값을 필드로 선언
    m_dB: np.ndarray = field(default_factory=lambda: np.array([5.782, 7.083, -0.983, 0.023, -4.401, -4.312]))

    # dB → 선형 변환 (default_factory 사용)
    m: np.ndarray = field(init=False)  # 자동 초기화 필드

    def __post_init__(self):
        """dB 값을 선형 값으로 변환 (m_dB → m)"""
        self.m = 10 ** (self.m_dB / 10)

# ✅ Learning parameters
params = CommunicationParams()
subset_size = 60000 // params.K  # 각 사용자당 할당할 데이터 개수
batch = 100
batch_set = [12, 10, 8, 4, 2]
num_epochs = 100
learning_rate = 0.01

device = "cuda" if torch.cuda.is_available() else "cpu"
# ✅ 1. 데이터셋 로딩 (MNIST)
train_loaders, test_loader = Dataset.MNIST_dataloader(batch, batch_set, subset_size, params)

# ✅ K명의 사용자 모델 생성 (각 사용자가 자기 모델 학습)
models = [FL.FNN().to(device) for _ in range(params.K)]
optimizers = [optim.Adam(models[k].parameters(), lr=learning_rate) for k in range(params.K)]
criterions = [nn.CrossEntropyLoss() for _ in range(params.K)]

# ✅ 중앙 서버(기지국) 모델 초기화
central_model = FL.FNN().to(device)

# ✅ Federated Learning 학습
for epoch in range(num_epochs):
    print(f"\n=== Epoch {epoch + 1}/{num_epochs} ===")

    h_k = (np.random.normal(0, params.sigma, (params.K, params.N))  # downlink CH; channel of kth users in nth BW
         + 1j * np.random.normal(0, params.sigma, (params.K, params.N))) / np.sqrt(2)

    h_k_random_allocation = (np.random.normal(0, params.sigma, (params.K))  # downlink CH; channel of kth users in nth BW
         + 1j * np.random.normal(0, params.sigma, (params.K))) / np.sqrt(2)

    # allocating power and CH through communication model
    P_i_opt, h_k_allocated, x = Communication.model(h_k, batch_set, params)

    ##### FL process #####
    num_joining_data = 0
    total_loss = 0
    total_grad_sum = None
    # local training and communication
    for k in range(params.K):
        # ✅ local_NN 호출 시 total_grad_sum을 인자로 넘김
        loss, join, total_grad_sum = FL.local_NN(k, h_k_allocated, P_i_opt, x, batch_set, params, train_loaders, models, optimizers,
                                              criterions, device, total_grad_sum)
        total_loss += loss
        num_joining_data += join * batch_set[k]

    # aggregating local models into central model
    FL.central_NN(num_joining_data, optimizers, central_model, models, learning_rate, total_grad_sum, total_loss, params)


# ✅ 중앙 모델(BS) 테스트 정확도 측정
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
print(f"\n=== Central model Test Accuracy: {accuracy:.2f}% ===")

models[0].eval()
correct = 0
total = 0
with torch.no_grad():
    for images, labels in test_loader:
        images, labels = images.to(device), labels.to(device)
        outputs = models[0](images)
        _, predicted = torch.max(outputs, 1)
        # print(labels, predicted.data)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

accuracy = 100 * correct / total
print(f"\n=== Individual model Test Accuracy: {accuracy:.2f}% ===")
