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

# ✅ 환경 변수 설정
@dataclass
# ✅ Communication parameters
class CommunicationParams:
    K: int = 10  # 사용자 수
    M: int = 10 # orthogonal sub-CHs
    q: int = 10  # size of gradient vector
    Nb: int = 1  # 기지국 안테나 수
    Nu: int = 1  # UE 안테나 수
    g_th: float = 0.3 # truncation threshold
    sigma: float = 1  # 노이즈 파워
    P0: int = 10  # TX 파워
    N0: float = field(default=10 ** (-174 / 10) * 10 ** (-3))  # 노이즈 밀도

# ✅ Learning parameters
params = CommunicationParams()
subset_size = int(60000 // params.K)  # 각 사용자당 할당할 데이터 개수
batch = 100
batch_set = []
for k in range(params.K): # set mini batch size to each user
    batch_set.append(batch)

num_epochs = 100
learning_rate = 1/np.sqrt(num_epochs)

device = "cuda" if torch.cuda.is_available() else "cpu"
# ✅ 1. 데이터셋 로딩 (MNIST)
train_loaders, test_loader = Dataset.MNIST_dataloader(batch, batch_set, subset_size, params)

# ✅ K명의 사용자 모델 생성 (각 사용자가 자기 모델 학습)
models = [FL.FNN().to(device) for _ in range(params.K)]
optimizers = [optim.Adam(models[k].parameters(), lr=learning_rate) for k in range(params.K)]
criterions = [nn.CrossEntropyLoss() for _ in range(params.K)]

# ✅ 중앙 서버(기지국) 모델 초기화
central_model = FL.FNN().to(device)
num_weights = central_model.count_parameters()
network_structure = central_model.get_gradient_structure
print(network_structure)
print(num_weights)
# ✅ Federated Learning 학습
for epoch in range(num_epochs):
    print(f"\n=== Epoch {epoch + 1}/{num_epochs} ===")

    ##### FL process #####
    num_joining_data = 0
    total_loss = 0
    g_k_tilde = np.zeros((params.K, num_weights, 1))
    shapes = None
    total_grad_sum = None
    for k in range(params.K):
        # ✅ local_NN 호출
        loss, total_grad, total_grad_sum = FL.local_NN(k, train_loaders, models, optimizers, criterions, device, total_grad_sum)
        total_loss += loss
        num_joining_data += batch_set[k]

        g_k_tilde[k], shapes = Functions.extract_and_sign(total_grad)
    num_iter = num_weights // params.M
    g_tilde = np.zeros((num_weights, 1))
    for iter in range(num_iter):
    ##### communication process #####
        h_k = (np.random.normal(0, 1, (params.K, params.M))  # downlink CH; channel of kth users in nth BW
             + 1j * np.random.normal(0, 1, (params.K, params.M))) / np.sqrt(2)

        # power allocation
        p_k = Communication.power_allocation(h_k, params)

        # pre-coding & transmitting
        g_tilde_t = np.zeros((params.M, 1), dtype=float)
        for k in range(params.K):
            g_k_tilde_Tr_iter = g_k_tilde[k][iter*params.M:(iter+1)*params.M]
            g_k_tilde_Tr_k = (g_k_tilde_Tr_iter * p_k[k]).astype(np.float64)
            g_tilde_t += g_k_tilde_Tr_k

    # aggregation
        g_tilde[iter*params.M:(iter+1)*params.M] = g_tilde_t

    # noise
    z = (np.random.normal(0, params.sigma, (num_weights, 1))  # downlink CH; channel of kth users in nth BW
         + 1j * np.random.normal(0, params.sigma, (num_weights, 1))) / np.sqrt(2)

    # majority-vote decoding
    v = Communication.majority_vote_decoder(g_tilde + z)
    #print(f"Shapes dictionary: {shapes}")  # shapes 값 출력

    ##### Central model update
    v_structure = Functions.restore_original_shape(v, shapes)
    #print(v_structure)
    # aggregating local models into central model
    FL.central_NN(num_joining_data, optimizers, central_model, models, learning_rate, v_structure, total_loss, params)


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
