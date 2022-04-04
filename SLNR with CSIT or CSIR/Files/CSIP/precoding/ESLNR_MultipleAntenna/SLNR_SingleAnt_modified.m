function u = SLNR_SingleAnt_modified(MT,R,K,KK,H,sigma2,FB_bits)

err_UpperBound=2^(-FB_bits/(MT-1));

% 根据调度结果，读入被调度者的信道系数
H_now=zeros(K,MT);
for i1=1:K,
    H_now(i1,:)=H((KK(1,i1)-1)*R+1,:);
end

% pre-coding
u=zeros(MT,K);      % single stream
Hi=zeros(1,MT);
Hk=zeros(1,MT);
temp1=zeros(1,MT);
for i1=1:K,
    Hi=H_now(i1,:);
    C=(MT-MT*err_UpperBound)*Hi'*Hi+err_UpperBound*eye(MT);

    D=sigma2*eye(MT);
    for i2=1:K,
        if i2~=i1,
            Hk=H_now(i2,:);
            D=D+(MT-MT*err_UpperBound)*Hk'*Hk+err_UpperBound*eye(MT);
        end
    end
    
    [eig_vec,eig_value]=eig(C,D);
    for n=1:MT,
        temp1(1,n)=eig_value(n,n);
    end
    [max_value,max_pos]=max(temp1);
    v=eig_vec(:,max_pos);
    u(:,i1)=v/(sqrt((v')*v));
end