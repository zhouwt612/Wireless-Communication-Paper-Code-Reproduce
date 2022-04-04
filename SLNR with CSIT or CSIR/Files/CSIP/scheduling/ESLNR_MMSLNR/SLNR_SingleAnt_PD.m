function u = SLNR_SingleAnt_PD(MT,R,K,KK,H,sigma2);
% only the first receiving antenna is used

% 根据调度结果，读入被调度者的信道系数
H_now=zeros(K,MT);
for i1=1:K,
    H_now(i1,:)=H((KK(1,i1)-1)*R+1,:);
end

% pre-coding
% The transmitter only knows the quantized channel correlation matrix.
u=zeros(MT,K);      % single stream
Hi=zeros(1,MT);
Hk=zeros(1,MT);
temp1=zeros(1,MT);
SLNR_reciprocal=zeros(1,K);
for i1=1:K,
    Hi=H_now(i1,:);
    C=Hi'*Hi;

    D=sigma2*eye(MT);
    for i2=1:K,
        if i2~=i1,
            Hk=H_now(i2,:);
            D=D+Hk'*Hk;
        end
    end
    
    [eig_vec,eig_value]=eig(C,D);
    for n=1:MT,
        temp1(1,n)=eig_value(n,n);
    end
    [max_val,max_pos]=max(temp1);
    SLNR_reciprocal(1,i1)=1/max_val;
    v=eig_vec(:,max_pos);
    u(:,i1)=v/(sqrt((v')*v));
end

p=zeros(K,K);
for i1=1:K,
    p(i1,i1)=K*SLNR_reciprocal(1,i1)/sum(SLNR_reciprocal);
end

u=u*sqrt(p);