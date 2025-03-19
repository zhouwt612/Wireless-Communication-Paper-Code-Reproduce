import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
import numpy as np
import random
import matplotlib
matplotlib.use("TkAgg")  # TkAgg 백엔드 사용
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader, SubsetRandomSampler
from dataclasses import dataclass, field

import Functions
import Communication
import FL
import Dataset
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


# ✅ 환경 변수 설정
class CommunicationParams:
    K = 10 # 사용자 수
    M = 10  # orthogonal sub-CHs
    q = 10  # size of gradient vector
    Nb = 1  # 기지국 안테나 수
    Nu = 1  # UE 안테나 수
    g_th = 0.3  # truncation threshold
    sigma = 1  # 노이즈 파워
    P0 = 10  # TX 파워
    N0 = 10 ** (-174 / 10) * 10 ** (-3)  # 노이즈 밀도


# ✅ Learning parameters
params = CommunicationParams()
subset_size = int(60000 // params.K)
batch = 100
#batch_set = [12, 10, 8, 4, 2]
batch_set = []
for k in range(params.K):  # set mini batch size to each user
    batch_set.append(batch)
print(batch_set)
num_epochs = 100
learning_rate = 0.01

device = "cuda" if torch.cuda.is_available() else "cpu"

# ✅ 데이터셋 로딩 (MNIST)
train_loaders, test_loader = Dataset.MNIST_dataloader(100, batch_set, subset_size, params)

# ✅ 중앙 모델 및 사용자 모델 초기화
models = [FL.FNN().to(device) for _ in range(params.K)]
optimizers = [optim.Adam(models[k].parameters(), lr=learning_rate) for k in range(params.K)]
criterions = [nn.CrossEntropyLoss() for _ in range(params.K)]
central_model = FL.FNN().to(device)

num_weights = central_model.count_parameters()

# ✅ 두 가지 실험을 위한 정확도 기록 리스트
central_model_accuracies_v_structure = []
central_model_accuracies_total_grad_sum = []

for mode in ["v_structure"]:
    print(f"\n=== Running Experiment with X = {mode} ===")

    central_model.load_state_dict(FL.FNN().state_dict())  # 중앙 모델 초기화
    accuracies = []

    for epoch in range(num_epochs):
        print(f"\n=== Epoch {epoch + 1}/{num_epochs} ===")

        total_grad_sum = None
        num_joining_data = 0
        total_loss = 0
        g_k_tilde = np.zeros((params.K, num_weights, 1))
        shapes = None

        for k in range(params.K):
            loss, total_grad, total_grad_sum = FL.local_NN(
                k, train_loaders, models, optimizers, criterions, device, total_grad_sum
            )
            total_loss += loss
            num_joining_data += batch_set[k]
            g_k_tilde[k], shapes = Functions.extract_and_sign(total_grad)

        num_iter = num_weights // params.M
        g_tilde = np.zeros((num_weights, 1))

        for iter in range(num_iter):
            h_k = (np.random.normal(0, 1, (params.K, params.M))
                   + 1j * np.random.normal(0, 1, (params.K, params.M))) / np.sqrt(2)

            p_k = params.P0 * np.ones((params.K, params.M, 1))
            g_tilde_t = np.zeros((params.M, 1), dtype=float)

            for k in range(params.K):
                g_k_tilde_Tr_iter = g_k_tilde[k][iter * params.M:(iter + 1) * params.M]
                g_k_tilde_Tr_k = (g_k_tilde_Tr_iter * p_k[k]).astype(np.float64)
                g_tilde_t += g_k_tilde_Tr_k

            g_tilde[iter * params.M:(iter + 1) * params.M] = g_tilde_t

        z = (np.random.normal(0, params.sigma, (num_weights, 1))
             + 1j * np.random.normal(0, params.sigma, (num_weights, 1))) / np.sqrt(2)

        v = Communication.majority_vote_decoder(g_tilde + z)
        v_structure = Functions.restore_original_shape(v, shapes)

        # ✅ X 값 변경 (v_structure 또는 total_grad_sum)
        X = v_structure if mode == "v_structure" else total_grad_sum
        FL.central_NN(num_joining_data, optimizers, central_model, models, learning_rate, X, total_loss, params)

        # ✅ 중앙 모델 정확도 측정
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
        print(correct)
        print(total)
        accuracy = 100 * correct / total
        accuracies.append(accuracy)
        print(f"\n=== Central model Test Accuracy ({mode}): {accuracy:.2f}% ===")

    # ✅ 결과 저장
    if mode == "v_structure":
        central_model_accuracies_v_structure = accuracies
    else:
        central_model_accuracies_total_grad_sum = accuracies

# ✅ Accuracy 비교 플롯
plt.plot(range(1, num_epochs + 1), central_model_accuracies_v_structure, marker='o', linestyle='-',
         label="Digital OTA")
# plt.plot(range(1, num_epochs + 1), central_model_accuracies_total_grad_sum, marker='s', linestyle='--',
#          label="OTA")
plt.xlabel("Epoch")
plt.ylabel("Accuracy (%)")
plt.ylim([0, 100])
plt.title(f"g_th={params.g_th}, SNR={int(10*np.log10(params.P0/params.N0))}, batch_size={batch_set}")
plt.legend()
plt.grid()
plt.show()
