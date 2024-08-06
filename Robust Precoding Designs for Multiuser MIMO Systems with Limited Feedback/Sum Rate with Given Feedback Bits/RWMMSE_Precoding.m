function [Pk,Pk_itea,Rate_Itea] = RWMMSE_Precoding(Hhatk,H,Pinit,gamma,mu,rho,sigma2,itea_num)
% Robust WMMSE precoding
%This code refers to the following scientific article:
%
% Wentao Zhou, Di Zhang, Mérouane Debbah, and Inkyu Lee,
% "Robust Precoding Designs for Multiuser MIMO Systems with Limited Feedback,
%" IEEE Transactions on Wireless Communications, To appear.
% 
% This is version 1.0 (last edited: 2024-04-08)
% 
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
% Input:
% Hhatk: the quantized channels (size: M,N,K)
% H: real channels (for the calculation of achievable rate)
% Pinit: initial precoding matrices
% gamma: quantization distortion
% mu: weight matrix (size: 1,K)
% rho: transmit power
% sigma2: noise power
% itea_num: # of iterations
% Output:
% Pk: precoding matrices
% Pk_itea: precoding matrices concatenations
% Rate_Itea: concatenation of each iteration
[M,N,K] = size(Hhatk);
[delta,beta] = Amplitudes(M,N,gamma);
P = zeros(M,M);
for idx1 = 1:1:K
    P(:,(idx1-1)*N+1:idx1*N) = Pinit(:,:,idx1);
end
Pk_itea = zeros(M,N,K,itea_num+1);
Pk_itea(:,:,:,1) = Pinit;
for idx2 = 1:1:itea_num
    D = ReceiveFilterD(Hhatk,P,delta,beta,sigma2);
    W = WeightMatrix(Hhatk,P,delta,beta,sigma2,mu);
    P = TransmitFilterP(Hhatk,D,W,delta,beta,sigma2,rho);
    for idx4 = 1:1:K
        Pk_itea(:,:,idx4,idx2+1) = P(:,(idx4-1)*N+1:idx4*N);
    end
end
Pk = zeros(M,N,K);
for idx3 = 1:1:K
    Pk(:,:,idx3) = P(:,(idx3-1)*N+1:idx3*N);
end
Rate_Itea = zeros(1,itea_num+1);
for idx5 = 1:1:itea_num+1
    Rate_Itea(idx5) = SumRateMIMOforK(H,Pk_itea(:,:,:,idx5));
end





function [delta,beta] = Amplitudes(M,N,gamma)
% Calculate the amplitudes of the quantized CSI and its quantization
% distortion.
% delta: the amplitude for Hhatk
% beta: the amplitude for Sk
delta = sqrt(M-M^2*gamma/(M-N));
beta = sqrt(M^2*gamma/(M-N));

function D = ReceiveFilterD(Hhatk,P,delta,beta,sigma2)
% Hhatk: the quantized CSI (size: M,N,K)
% Pk: the concatenation of precoding matrices (size: M,NK)
% delta: amplitude for quantized CSI
% beta: amplitude for quantization distortion
% sigma2: noise power
[M,N,K] = size(Hhatk);
% Calculate Phi_k for all K UEs
% Every UE has the same Phi_k, we calculate it once
Phi_k = diag(((1/M)*ones(N,M))*diag(P*P'));
D = zeros(N*K,N*K);
for idx = 1:1:K
    D((idx-1)*N+1:idx*N,(idx-1)*N+1:idx*N) = delta*P(:,(idx-1)*N+1:idx*N)'*Hhatk(:,:,idx)*...
        inv(delta^2*Hhatk(:,:,idx)'*P*P'*Hhatk(:,:,idx)+beta^2*Phi_k+sigma2*eye(N));
end

function W = WeightMatrix(Hhatk,P,delta,beta,sigma2,mu)
% Hhatk: the quantized CSI (size: M,N,K)
% Pk: the concatenation of precoding matrices (size: M,NK)
% delta: amplitude for quantized CSI
% beta: amplitude for quantization distortion
% sigma2: noise power
[M,N,K] = size(Hhatk);
% Calculate Phi_k for all K UEs
% Every UE has the same Phi_k, we calculate it once
Phi_k = diag(((1/M)*ones(N,M))*diag(P*P'));
Ebark = zeros(N,N,K);
for idx1 = 1:1:K
    Ebark(:,:,idx1) = eye(N) - delta^2*P(:,(idx1-1)*N+1:idx1*N)'*Hhatk(:,:,idx1)*...
        inv(delta^2*Hhatk(:,:,idx1)'*P*P'*Hhatk(:,:,idx1)+beta^2*Phi_k+sigma2*eye(N))*...
        Hhatk(:,:,idx1)'*P(:,(idx1-1)*N+1:idx1*N);
end
W = zeros(N*K,N*K);
for idx2 = 1:1:K
    W((idx2-1)*N+1:idx2*N,(idx2-1)*N+1:idx2*N) = mu(idx2)*inv(Ebark(:,:,idx2));
end

function P = TransmitFilterP(Hhatk,D,W,delta,beta,sigma2,rho)
% Hhatk: the quantized CSI (size: M,N,K)
% P: the concatenation of precoding matrices (size: M,NK)
% delta: amplitude for quantized CSI
% beta: amplitude for quantization distortion
% sigma2: noise power
[M,N,K] = size(Hhatk);
G = (1/M)*ones(M,N*K)*diag(D'*W*D);
dgG = diag(G);
Hhat = zeros(M,N*K);
for idx1 = 1:1:K
    Hhat(:,(idx1-1)*N+1:idx1*N) = Hhatk(:,:,idx1);
end
Pbar = inv(delta^2*Hhat*D'*W*D*Hhat'+beta^2*dgG+(sigma2*trace(W*D*D'))/rho*eye(M))*Hhat*D'*W;
P = sqrt(rho/trace(Pbar*Pbar'))*Pbar;

function R = SumRateMIMOforK(H,P)
% H: channel realization (size: (M,N,K))
% P: precoder matrix (size: (M,N,K))
% R: sum rate of MIMO with K users
[M,N,K] = size(H);
interference = zeros(N,N,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        if idx2 ~= idx1
            interference(:,:,idx1) = interference(:,:,idx1) + ...
                H(:,:,idx1)'*P(:,:,idx2)*P(:,:,idx2)'*H(:,:,idx1);
        else
            interference(:,:,idx1) = interference(:,:,idx1) + zeros(N);
        end
    end
end
interference = interference + eye(N);
T = zeros(1,K);
for idx3 = 1:1:K
    T(idx3) = log2(det(eye(N) + ...
        P(:,:,idx3)'*H(:,:,idx3)/interference(:,:,idx3)*H(:,:,idx3)'*P(:,:,idx3)));
end
R = sum(real(T));