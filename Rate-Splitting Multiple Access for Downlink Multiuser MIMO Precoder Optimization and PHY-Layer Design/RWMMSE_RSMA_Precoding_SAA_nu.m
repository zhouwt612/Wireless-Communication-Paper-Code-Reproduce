function [Pc,Pp,Pc_itea,Pp_itea] = RWMMSE_RSMA_Precoding_SAA_nu(Hhatk,Pcinit,Ppinit,rho,sigma2,sigmae,samp_num,nu)
itea_num = 1000;
[M,N,K] = size(Hhatk);
Pc_itea = zeros(M,N,itea_num+1);
Pp_itea = zeros(M,N,K,itea_num+1);
Pc = Pcinit;
Pp = Ppinit;
Pc_itea(:,:,1) = Pc;
Pp_itea(:,:,:,1) = Pp;
Hhat_Gen = zeros(M,N,K,samp_num);
for idx1 = 1:1:samp_num
    Hhat_Gen(:,:,:,idx1) = Hhatk + sqrt(sigmae)*(randn(M,N,K) + 1i*randn(M,N,K))/sqrt(2);
end
for idx3 = 1:1:itea_num
    [Pc,Pp] = TransmitFilterCVX_SAA(Hhat_Gen,Pc,Pp,rho,sigma2);
    Pc_itea(:,:,idx3+1) = Pc;
    Pp_itea(:,:,:,idx3+1) = Pp;
    Diff = 0;
    DiffC = Pc_itea(:,:,idx3) - Pc;
    Diff = Diff + trace(DiffC'*DiffC);
    for idx4 = 1:1:K
        DiffP = Pp_itea(:,:,idx4,idx3) - Pp(:,:,idx4);
        Diff = Diff + trace(DiffP'*DiffP);
    end
    if Diff <= nu
        break
    end
end




function [Pc,Pp] = TransmitFilterCVX_SAA(Hhat_Gen,Pc,Pp,rho,sigma2)
[M,N,K,samp_num] = size(Hhat_Gen);
Pp_con = zeros(M,N*K);
for IDX0 = 1:1:K
    Pp_con(:,(IDX0-1)*N+1:IDX0*N) = Pp(:,:,IDX0);
end
P = [Pc,Pp_con];
% Generation of Dck
Dck = zeros(N,N,K,samp_num);
for IDX1 = 1:1:samp_num
    for IDX2 = 1:1:K
        Dck(:,:,IDX2,IDX1) = Pc'*Hhat_Gen(:,:,IDX2,IDX1)*...
            inv(Hhat_Gen(:,:,IDX2,IDX1)'*P*P'*Hhat_Gen(:,:,IDX2,IDX1)+sigma2*eye(N));
    end
end
% Generation of Dpk
for IDX3 = 1:1:samp_num
    for IDX4 = 1:1:K
        Dpk(:,:,IDX4,IDX3) = Pp(:,:,IDX4)'*Hhat_Gen(:,:,IDX4,IDX3)*...
            inv(Hhat_Gen(:,:,IDX4,IDX3)'*Pp_con*Pp_con'*Hhat_Gen(:,:,IDX4,IDX3)+sigma2*eye(N));
    end
end
% Generation of Wck
Wck = zeros(N,N,K,samp_num);
for IDX5 = 1:1:samp_num
    for IDX6 = 1:1:K
        Wck(:,:,IDX6,IDX5) = inv(eye(N) - Dck(:,:,IDX6,IDX5)*Hhat_Gen(:,:,IDX6,IDX5)'*Pc);
    end
end
% Generation of Wpk
Wpk = zeros(N,N,K,samp_num);
for IDX7 = 1:1:samp_num
    for IDX8 = 1:1:K
        Wpk(:,:,IDX8,IDX7) = inv(eye(N) - Dpk(:,:,IDX8,IDX7)*Hhat_Gen(:,:,IDX8,IDX7)'*Pp(:,:,IDX8));
    end
end
% Calculation of P=[Pc,Pp]
cvx_begin quiet
    % variables
    variable Pc(M,N) complex
    variable Pp(M,N,K) complex
    variable xi(K)
    
    % Preparation
    A_ck = cvx(zeros(M*N,M*N,K));
    a_ck = cvx(zeros(M*N,1,K));
    A_pk = cvx(zeros(M*N,M*N,K));
    a_pk = cvx(zeros(M*N,1,K));
    Phi_ck = zeros(1,1,K);
    Phi_pk = zeros(1,1,K);
    for idx1 = 1:1:samp_num
        for idx2 = 1:1:K
            A_ck_block = Hhat_Gen(:,:,idx2,idx1)*Dck(:,:,idx2,idx1)'*Wck(:,:,idx2,idx1)*Dck(:,:,idx2,idx1)*Hhat_Gen(:,:,idx2,idx1)';
            A_ck_block = 0.5*(A_ck_block + A_ck_block');
            A_ck(:,:,idx2) = A_ck(:,:,idx2) + kron(eye(N),A_ck_block);
            a_ck(:,:,idx2) = a_ck(:,:,idx2) + vec(Hhat_Gen(:,:,idx2,idx1)*Dck(:,:,idx2,idx1)'*Wck(:,:,idx2,idx1));
            Apk_block = Hhat_Gen(:,:,idx2,idx1)*Dpk(:,:,idx2,idx1)'*Wpk(:,:,idx2,idx1)*Dpk(:,:,idx2,idx1)*Hhat_Gen(:,:,idx2,idx1)';
            Apk_block = 0.5*(Apk_block + Apk_block');
            A_pk(:,:,idx2) = A_pk(:,:,idx2) + kron(eye(N),Apk_block);
            a_pk(:,:,idx2) = a_pk(:,:,idx2) + vec(Hhat_Gen(:,:,idx2,idx1)*Dpk(:,:,idx2,idx1)'*Wpk(:,:,idx2,idx1));
            Phi_ck(:,:,idx2) = Phi_ck(:,:,idx2) + sigma2*trace(Wck(:,:,idx2,idx1)*Dck(:,:,idx2,idx1)*Dck(:,:,idx2,idx1)')+...
                trace(Wck(:,:,idx2,idx1)) - log2(det(Wck(:,:,idx2,idx1)));
            Phi_pk(:,:,idx2) = Phi_pk(:,:,idx2) + sigma2*trace(Wpk(:,:,idx2,idx1)*Dpk(:,:,idx2,idx1)*Dpk(:,:,idx2,idx1)')+...
                trace(Wpk(:,:,idx2,idx1)) - log2(det(Wpk(:,:,idx2,idx1)));
        end
    end
    A_ck = A_ck/samp_num;
    a_ck = a_ck/samp_num;
    A_pk = A_pk/samp_num;
    a_pk = a_pk/samp_num;
    Phi_ck = Phi_ck/samp_num;
    Phi_pk = Phi_pk/samp_num;
    
    obj_sum1 = 0;
    obj_sum2 = 0;
    obj_sum3 = 0;
    for idx3 = 1:1:K
        for idx4 = 1:1:K
            obj_sum1 = obj_sum1 + vec(Pp(:,:,idx4))'*A_pk(:,:,idx3)*vec(Pp(:,:,idx4));
        end
        obj_sum2 = obj_sum2 + (vec(Pp(:,:,idx3))'*a_pk(:,:,idx3) + a_pk(:,:,idx3)'*vec(Pp(:,:,idx3)));
        obj_sum3 = obj_sum3 + real(Phi_pk(:,:,idx3));
    end
    
    cons_sum1 = cvx(zeros(1,K));
    cons_sum2 = cvx(zeros(1,K));
    cons_sum3 = cvx(zeros(1,K));
    cons_sum4 = cvx(zeros(1,K));
    for idx5 = 1:1:K
        cons_sum1(idx5) = vec(Pc)'*A_ck(:,:,idx5)*vec(Pc);
        for idx6 = 1:1:K
            cons_sum2(idx5) = cons_sum2(idx5) + vec(Pp(:,:,idx6))'*A_ck(:,:,idx5)*vec(Pp(:,:,idx6));
        end
        cons_sum3(idx5) = vec(Pc)'*a_ck(:,:,idx5) + a_ck(:,:,idx5)'*vec(Pc);
        cons_sum4(idx5) = real(Phi_ck(:,:,idx5));
    end
    Pp_pow = 0;
    for idx7 = 1:1:K
        Pp_pow = Pp_pow + vec(Pp(:,:,idx7))'*vec(Pp(:,:,idx7));
    end
    minimize( sum(xi) + obj_sum1 - obj_sum2 + obj_sum3 )
    subject to
        for idx8 = 1:1:K
            cons_sum1(idx8) + cons_sum2(idx8) - cons_sum3(idx8) + cons_sum4(idx8) <= sum(xi) + N;
            xi(idx8) <= 0;
        end
        vec(Pc)'*vec(Pc) + Pp_pow <= rho;
%         xi <= 0;
cvx_end