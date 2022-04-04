function [u, p] = SINR_new(N,Mi,S,K,H,sigma2)

channel_correlation=zeros(N,N,K);
HHH=zeros(Mi,N);
for i1=1:K,
    HHH=H((i1-1)*Mi+1:i1*Mi,:);
    channel_correlation(:,:,i1)=(HHH'*HHH)/(Mi*sigma2);
end

% pre-coding
u=zeros(N,K);
p=zeros(K,1);
q=zeros(K,1);
lamda_old=0;
lamda_new=0;
sigma2_new=ones(K,1);   % normalized noise power
D=zeros(K,K);
fy=zeros(K,K);
while lamda_old-lamda_new>=10^(-10) | lamda_old==0,     % convergence or initialization
    lamda_old=lamda_new;
    for i1=1:K,
        Q=eye(N);
        for i2=1:K
            if i2~=i1,
                Q=Q+q(i2,1)*channel_correlation(:,:,i2);
            end
        end
        
        eig_vec = zeros(N,N);
        eig_val = zeros(N,N);
        [eig_vec,eig_val] = eig(channel_correlation(:,:,i1),Q);
        
        temp1=zeros(1,N);
        for i2=1:N,
            temp1(1,i2)=eig_val(i2,i2);
        end
        
        [max_val,max_pos] = max(temp1);
        u(:,i1)=eig_vec(:,max_pos);
        u(:,i1)=u(:,i1)/sqrt(trace(u(:,i1)'*u(:,i1)));
    end
    
    for i1=1:K,
        D(i1,i1)=1/(u(:,i1)'*channel_correlation(:,:,i1)*u(:,i1));
        for i2=1:K,
            if i2~=i1,
                fy(i1,i2)=u(:,i2)'*channel_correlation(:,:,i1)*u(:,i2);
            end
        end
    end
    
    A=zeros(K+1,K+1);
    A=[D*fy.' D*sigma2_new;(1/K)*ones(1,K)*D*fy.' (1/K)*ones(1,K)*D*sigma2_new];
    eig_vec=zeros(K+1,K+1);
    eig_val=zeros(K+1,K+1);
    [eig_vec,eig_val]=eig(A);
    
    temp1=zeros(1,K+1);
    for i1=1:(K+1),
        temp1(1,i1)=eig_val(i1,i1);
    end
    
    [max_val,max_pos]=max(temp1);
    q=eig_vec(1:K,max_pos)/eig_vec(K+1,max_pos);
    lamda_new=max_val;
end

B=[D*fy D*sigma2_new;(1/K)*ones(1,K)*D*fy (1/K)*ones(1,K)*D*sigma2_new];
eig_vec=zeros(K+1,K+1);
eig_val=zeros(K+1,K+1);
[eig_vec,eig_val]=eig(B);

temp1=zeros(1,K+1);
for i1=1:(K+1),
    temp1(1,i1)=eig_val(i1,i1);
end

[max_val,max_pos]=max(temp1);
p=eig_vec(1:K,max_pos)/eig_vec(K+1,max_pos);