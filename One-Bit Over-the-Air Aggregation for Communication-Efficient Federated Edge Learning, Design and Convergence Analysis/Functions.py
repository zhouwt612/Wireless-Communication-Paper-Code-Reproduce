import numpy as np
import scipy.integrate as integrate
import torch
# Define the function to integrate with lower limit k
def exponential_integratation(k):
    result, error = integrate.quad(lambda x: (1/x) * np.exp(-x), k, np.inf)
    return result, error


# total_grad_sum을 sign 벡터로 변환하는 함수
def extract_and_sign(total_grad_sum):
    shapes = {name: grad.shape for name, grad in total_grad_sum.items()}  # 각 gradient의 shape 저장
    grad_values = [grad.flatten().cpu() for grad in total_grad_sum.values()]  # GPU에서 CPU로 이동
    grad_vector = torch.cat(grad_values).view(-1, 1)  # 열 벡터 변환
    signed_vector = torch.sign(grad_vector)  # sign 적용
    vector_length = signed_vector.numel()  # 벡터의 길이 계산

    #print(f"Signed Vector Length: {vector_length}")  # 벡터의 길이 출력
    return signed_vector, shapes
 # sign 벡터와 원래 shape 반환

# sign 벡터를 원래 구조로 복원하는 함수
def restore_original_shape(signed_vector, shapes):
    torch_vector = torch.from_numpy(signed_vector)
    restored_dict = {}
    idx = 0
    for name, shape in shapes.items():
        num_elements = torch.prod(torch.tensor(shape))  # 요소 개수
        restored_dict[name] = torch_vector[idx:idx + num_elements].view(shape)  # 원래 shape로 변환
        idx += num_elements

    return restored_dict

# 예제 total_grad_sum (가상의 gradient 값)
# total_grad_sum = {
#     "layer1.weight": torch.tensor([[0.5, -0.3], [0.2, -0.7]]),
#     "layer1.bias": torch.tensor([0.1, -0.4]),
# }
#
# # Step 1: Sign 벡터 변환
# signed_grad_vector, shapes = extract_and_sign(total_grad_sum)
# print("Sign 벡터:\n", shapes)
#
# # Step 2: 원래 shape로 복원
# restored_total_grad_sum = restore_original_shape(signed_grad_vector, shapes)
# print("\n복원된 total_grad_sum:")
# for name, tensor in restored_total_grad_sum.items():
#     print(f"{name}: \n{tensor}")

