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
    train_loaders = []
    for i in range(params.K):
        subset_indices = indices[i * subset_size: (i + 1) * subset_size]  # 랜덤 인덱스 분배
        sampler = SubsetRandomSampler(subset_indices)
        train_loader = DataLoader(full_train_dataset, batch_size=batch_set[i], sampler=sampler)
        train_loaders.append(train_loader)

    for k in range(params.K):
        print(f"User {k + 1} has {len(train_loaders[k])} batches.")

    return train_loaders, test_loader