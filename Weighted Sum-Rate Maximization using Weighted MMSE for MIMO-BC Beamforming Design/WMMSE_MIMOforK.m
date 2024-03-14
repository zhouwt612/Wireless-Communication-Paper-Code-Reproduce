function [P,Ritea] = WMMSE_MIMOforK(H,P,pow,itea_num)
Ritea = zeros(1,itea_num+1);
Ritea(1) = SumRateMIMO(H,P);
for idx = 1:1:itea_num
    A = receive_filter(H,P);
    W = weight(H,P,A);
    P = precoding(H,A,W,pow);
    Ritea(idx+1) = SumRateMIMO(H,P);
end


function Ak = receive_filter(H,P)
[M,N,K] = size(H);
Ak = zeros(N,N,K);
Pcon = zeros(M,N*K);
for idx1 = 1:1:K
    Pcon(:,(idx1-1)*N+1:idx1*N) = P(:,:,idx1);
end
for idx2 = 1:1:K
    Ak(:,:,idx2) = P(:,:,idx2)'*H(:,:,idx2)*inv(H(:,:,idx2)'*Pcon*Pcon'*H(:,:,idx2)+eye(N));
end


function Wk = weight(H,P,A)
[M,N,K] = size(H);
Mk = zeros(N,N,K);
for idx1 = 1:1:K
    Mk(:,:,idx1) = eye(N) - A(:,:,idx1)*H(:,:,idx1)'*P(:,:,idx1);
end
Wk = zeros(N,N,K);
for idx2 = 1:1:K
    Wk(:,:,idx2) = inv(Mk(:,:,idx2));
end


function P = precoding(H,A,W,pow)
[M,N,K] = size(H);
Hcon = zeros(M,N*K);
Wcon = zeros(N*K,N*K);
Acon = zeros(N*K,N*K);
for idx1 = 1:1:K
    Hcon(:,(idx1-1)*N+1:idx1*N) = H(:,:,idx1);
    Wcon((idx1-1)*N+1:idx1*N,(idx1-1)*N+1:idx1*N) = W(:,:,idx1);
    Acon((idx1-1)*N+1:idx1*N,(idx1-1)*N+1:idx1*N) = A(:,:,idx1);
end
Pun = inv(Hcon*Acon'*Wcon*Acon*Hcon'+trace(Wcon*Acon*Acon')/pow*eye(M))*Hcon*Acon'*Wcon;
Pn = sqrt(pow/trace(Pun*Pun'))*Pun;
P = zeros(M,N,K);
for idx2 = 1:1:K
    P(:,:,idx2) = Pn(:,(idx2-1)*N+1:idx2*N);
end