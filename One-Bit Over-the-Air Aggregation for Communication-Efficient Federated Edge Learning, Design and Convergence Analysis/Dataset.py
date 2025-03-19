import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
import numpy as np
import random
from torch.utils.data import DataLoader, SubsetRandomSampler
from dataclasses import dataclass, field

def MNIST_dataloader(batch, batch_set, subset_size, params):
    transform = transforms.Compose([transforms.ToTensor()])
    full_train_dataset = torchvision.datasets.MNIST(root='./data', train=True, transform=transform, download=True)
    test_dataset = torchvision.datasets.MNIST(root='./data', train=False, transform=transform, download=True)
    test_loader = DataLoader(test_dataset, batch_size=batch, shuffle=False)

    # ✅ K명의 사용자 데이터 분배 (랜덤 샘플링 적용)
    indices = np.arange(len(full_train_dataset))
    np.random.shuffle(indices)  # 전체 데이터를 랜덤하게 섞음

    # ✅ 각 사용자에게 `SubsetRandomSampler`로 랜덤한 데이터 분배
    # ✅ 사용자별 동일한 데이터셋을 사용하기 위해 미리 샘플링된 데이터셋 저장
    train_datasets = []
    for k in range(params.K):  # 전체 MNIST 데이터셋 불러오기
        indices = list(range(k * subset_size, (k + 1) * subset_size))  # 고정된 인덱스 설정
        train_datasets.append(torch.utils.data.Subset(full_train_dataset, indices))

    # ✅ 고정된 데이터셋을 가진 DataLoader 생성 (shuffle=False)
    train_loaders = [torch.utils.data.DataLoader(dataset, batch_size=batch, shuffle=False) for dataset in
                     train_datasets]

    for k in range(params.K):
        print(f"User {k + 1} has {len(train_loaders[k])} batches.")

    return train_loaders, test_loader