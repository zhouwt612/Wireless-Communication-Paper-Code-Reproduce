function [Pc,Pp] = SNS_Precoding_nu(H,pow,sigma2,nu)
itea_num = 1000;
[M,N,K] = size(H);
Xk_la = Xk_Initial(M,N,K,pow);
Phik = SNS_BD(H);
Pc_itea = zeros(M,N,itea_num+1);
Pp_itea = zeros(M,N,K,itea_num+1);
for idx1 = 1:1:itea_num
    [Pc,Pp,Xk_la] = QX_Calculator(H,Phik,Xk_la,pow,sigma2);
    Pc_itea(:,:,idx1+1) = Pc;
    Pp_itea(:,:,:,idx1+1) = Pp;
    Diff = 0;
    DiffC = Pc_itea(:,:,idx1) - Pc;
    Diff = Diff + trace(DiffC'*DiffC);
    for idx2 = 1:1:K
        DiffP = Pp_itea(:,:,idx2,idx1) - Pp(:,:,idx2);
        Diff = Diff + trace(DiffP'*DiffP);
    end
    if Diff <= nu
        break
    end
end


function Phik = SNS_BD(H)
[M,N,K] = size(H);
Phik = cell(1,K);
Hcon = [];
for idx = 1:1:K
    if idx == 1
        Phik{idx} = eye(M);
    else
        Hcon = [Hcon H(:,:,idx-1)];
        Phik{idx} = null(Hcon');
    end
end


function Xk = Xk_Initial(M,N,K,pow)
Xk = cell(1,4);
for idx = 1:1:K
    Xk{idx} = pow/K/(M-(idx-1)*N)*eye(M-(idx-1)*N);
end


function [Pc,Pp,Xk,cvx_optval] = QX_Calculator(H,Phik,Xk_la,pow,sigma2)
[M,N,K] = size(H);
sig_inv = 1/sigma2;
cvx_begin quiet
    variable Qc(M,M) complex semidefinite
    variable X1(M,M) complex semidefinite
    variable X2(M-N,M-N) complex semidefinite
    variable X3(M-2*N,M-2*N) complex semidefinite
    variable X4(M-3*N,M-3*N) complex semidefinite
    
    Xkcon = cvx(zeros(M,M,K));
    Xkcon(:,:,1) = X1;
    Xkcon(1:M-N,1:M-N,2) = X2;
    Xkcon(1:M-2*N,1:M-2*N,3) = X3;
    Xkcon(1:M-3*N,1:M-3*N,4) = X4;
    
    Qk = cvx(zeros(M,M,K));
    for idx_pre_1 = 1:1:K
        Qk(:,:,idx_pre_1) = Phik{idx_pre_1}*Xkcon(1:M-N*(idx_pre_1-1),1:M-N*(idx_pre_1-1),idx_pre_1)*Phik{idx_pre_1}';
    end
    
    SumQk = cvx(zeros(M,M,K+1));
    Qk_breve = cvx(zeros(M,M,K+1));
    for idx_pre_2 = 1:1:K
        SumQk(:,:,idx_pre_2+1) = SumQk(:,:,idx_pre_2) + Qk(:,:,idx_pre_2);
        Qk_breve(:,:,idx_pre_2+1) = Qk_breve(:,:,idx_pre_2) + Phik{idx_pre_2}*Xk_la{idx_pre_2}*Phik{idx_pre_2}';
    end
    
    Rck = cvx(zeros(1,K));
    Rpk = cvx(zeros(1,K));
    for idx1 = 1:1:K
        Rck_minus = 0;
        for idx2 = 1:1:idx1-1
            Rck_minus = Rck_minus + sig_inv/log(2)*trace(Phik{idx2}'*H(:,:,idx1)...
                /(eye(N)+sig_inv*H(:,:,idx1)'*Qk_breve(:,:,idx1+1)*H(:,:,idx1))*H(:,:,idx1)'*Phik{idx2}*(Xkcon(1:M-N*(idx2-1),1:M-N*(idx2-1),idx2)-Xk_la{idx2}));
        end
        Rk_minus = 0;
        for idx3 = 1:1:idx1-1
            Rk_minus = Rk_minus + sig_inv/log(2)*trace(Phik{idx3}'*H(:,:,idx1)...
                /(eye(N)+sig_inv*H(:,:,idx1)'*Qk_breve(:,:,idx1)*H(:,:,idx1))*H(:,:,idx1)'*Phik{idx3}*(Xkcon(1:M-N*(idx3-1),1:M-N*(idx3-1),idx3)-Xk_la{idx3}));
        end
        A = eye(N)+sig_inv*H(:,:,idx1)'*(Qc+SumQk(:,:,idx1+1))*H(:,:,idx1);
        Rck(idx1) = 1/log(2)*log_det(0.5*(A+A'))...
            - 1/log(2)*log_det((eye(N)+sig_inv*H(:,:,idx1)'*(Qk_breve(idx1+1))*H(:,:,idx1))) - real(Rck_minus);
        B = eye(N)+sig_inv*H(:,:,idx1)'*SumQk(:,:,idx1+1)*H(:,:,idx1);
        Rpk(idx1) = 1/log(2)*log_det(0.5*(B+B'))...
            - 1/log(2)*log_det(real(eye(N)+sig_inv*H(:,:,idx1)'*(Qk_breve(idx1))*H(:,:,idx1))) - real(Rk_minus);
    end
    maximize( min(real(Rck)) + sum(real(Rpk)))
    subject to
        real(trace(Qc) + trace(X1) + trace(X2) + trace(X3) + trace(X4)) <= pow;
cvx_end
[Uc, Sc, ~] = svd(Qc);
Pc = Uc(:,1:N)*sqrt(Sc(1:N,1:N));
Pp = zeros(M,N,K);
[U1, S1, ~] = svd(X1);
Pp(:,:,1) = Phik{1}*U1(:,1:N)*sqrt(S1(1:N,1:N));
[U2, S2, ~] = svd(X2);
Pp(:,:,2) = Phik{2}*U2(:,1:N)*sqrt(S2(1:N,1:N));
[U3, S3, ~] = svd(X3);
Pp(:,:,3) = Phik{3}*U3(:,1:N)*sqrt(S3(1:N,1:N));
[U4, S4, ~] = svd(X4);
Pp(:,:,4) = Phik{4}*U4(:,1:N)*sqrt(S4(1:N,1:N));
Xk = cell(1,K);
Xk{1} = X1;
Xk{2} = X2;
Xk{3} = X3;
Xk{4} = X4;