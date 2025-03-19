import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
import numpy as np
import random
from torch.utils.data import DataLoader, SubsetRandomSampler


class FNN(nn.Module):
    def __init__(self):
        super(FNN, self).__init__()
        self.fc1 = nn.Linear(28 * 28, 50)  # 입력: 784(28x28), 은닉층: 50개 뉴런
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(50, 10)  # 출력: 10개 클래스 (0~9)

    def forward(self, x):
        x = x.view(-1, 28 * 28)  # 이미지 데이터를 1D 벡터로 변환
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x


def local_NN(k, h_k, P_allocated, x, batch_set, params, train_loaders, models, optimizers, criterions, device,
             total_grad_sum):
    SNR = P_allocated[k][x[k]] / (params.N0 * params.B_U)
    SINR = abs(h_k[k]) ** 2 * SNR
    p_e = 1 - np.exp(-params.m[3] / SINR)

    models[k].train()
    batch_iterator = iter(train_loaders[k])

    # ✅ 예외 처리: 데이터가 부족할 경우 방지
    try:
        images, labels = next(batch_iterator)
    except StopIteration:
        return 0, 0, total_grad_sum

    images, labels = images.to(device), labels.to(device)
    optimizers[k].zero_grad()
    outputs = models[k](images)
    loss = criterions[k](outputs, labels)
    loss.backward()
    #print(loss.item())
    optimizers[k].step()

    # ✅ Gradient 값 저장
    total_grad = {name: param.grad.clone() for name, param in models[k].named_parameters() if param.grad is not None}

    # ✅ main에서 받은 total_grad_sum이 None이면 초기화
    if total_grad_sum is None:
        total_grad_sum = {name: torch.zeros_like(grad) for name, grad in total_grad.items()}

    # ✅ CH modeling : p_e를 넘기면 데이터 전송 성공
    join = 0
    if random.uniform(0, 1) > p_e:
        join = 1
        for name in total_grad.keys():
            total_grad_sum[name] += batch_set[k] * total_grad[name]  # 사용자별 평균 Gradient를 누적, weigthed by batch size

    return loss.item() , join, total_grad_sum

def central_NN(num_joining_data, optimizers, central_model, models, learning_rate, total_grad_sum, total_loss, params):
    print(f'joining users {num_joining_data}, average loss: {total_loss / params.K}')
    # ✅ 중앙 모델(Central Model) 업데이트 (FedAvg 방식)
    central_optimizer = optim.Adam(central_model.parameters(), lr=learning_rate)  # 중앙 모델용 Optimizer 추가
    with torch.no_grad():
        for name, param in central_model.named_parameters():
            param.grad = total_grad_sum[name] / num_joining_data  # 모든 사용자 평균 Gradient 적용

    # ✅ Optimizer를 사용하여 중앙 모델 업데이트
    central_optimizer.step()
    central_optimizer.zero_grad()

    # ✅ 모든 사용자 모델을 중앙 모델과 동일한 가중치로 동기화
    with torch.no_grad():
        for k in range(params.K):
            for name, param in models[k].named_parameters():
                param.copy_(central_model.state_dict()[name])  # 중앙 모델의 가중치로 동기화

            optimizers[k].step()