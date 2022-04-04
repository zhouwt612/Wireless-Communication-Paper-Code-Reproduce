function [count,capacity] = SLNR(MT,K_all,K,R,H,sigma2,KK,x,w,Tc)

count = 0;
capacity = 0;

% 根据调度结果，读入被调度者的待传输数据以及信道系数
x_now=zeros(K,Tc);
w_now=zeros(R*K,Tc);
H_now=zeros(R*K,MT);
for i1=1:K,
    x_now(i1,:)=x(KK(1,i1),:);      % 各用户仍是只发送单路数据
    w_now((i1-1)*R+1:i1*R,:)=w((KK(1,i1)-1)*R+1:KK(1,i1)*R,:);      % 但噪声维数决定于各移动台端的天线数
    H_now((i1-1)*R+1:i1*R,:)=H((KK(1,i1)-1)*R+1:KK(1,i1)*R,:);
end

% pre-coding
u=zeros(MT,K);
temp1=zeros(1,MT);
Hi=zeros(R,MT);
Hk=zeros(R,MT);
for i1=1:K,
    Hi=H_now((i1-1)*R+1:i1*R,:);
    C=Hi'*Hi;

    D=R*sigma2*eye(MT);
    for i2=1:K,
        if i2~=i1,
            Hk=H_now((i2-1)*R+1:i2*R,:);
            D=D+Hk'*Hk;
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

y=H_now*u*x_now+w_now;
% matched filter
r=zeros(1,K*R);
temp2=zeros(1,R);
for i1=1:K,
    HHH=H_now((i1-1)*R+1:i1*R,:);
    temp2=(HHH*u(:,i1))';
    temp2=temp2/(sqrt(trace(temp2*temp2')));
    r(1,(i1-1)*R+1:i1*R)=temp2;    
end

yy=zeros(K,Tc);
for i1=1:K,
    temp3=zeros(1,K*R);
    temp3(1,(i1-1)*R+1:i1*R)=r(1,(i1-1)*R+1:i1*R);
    yi=temp3*y;
    yy(i1,:)=yi;
end

s_now=zeros(K,2*Tc);
for i1=1:K,
    for i2=1:Tc,
        if real(x_now(i1,i2))>=0,
            s_now(i1,i2)=0;
        else
            s_now(i1,i2)=1;
        end
        if imag(x_now(i1,i2))>=0,
            s_now(i1,i2+Tc)=0;
        else
            s_now(i1,i2+Tc)=1;
        end
    end
end
        
ss=zeros(K,2*Tc);
for i1=1:K,
    for i2=1:Tc,
        if real(yy(i1,i2))>=0,
            ss(i1,i2)=0;
        else
            ss(i1,i2)=1;
        end
        if imag(yy(i1,i2))>=0,
            ss(i1,i2+Tc)=0;
        else
            ss(i1,i2+Tc)=1;
        end
    end
end

for i1=1:K,
    for i2=1:2*Tc,
        if s_now(i1,i2)~=ss(i1,i2),
            count=count+1;
        end
    end
end

tempC=0;
for i1=1:K,
    temp4=R*sigma2;
    HHH=H_now((i1-1)*R+1:i1*R,:);
    for i2=1:K,
        if i2~=i1,
            temp4=temp4+u(:,i2)'*HHH'*HHH*u(:,i2);
        end
    end
    tempC=tempC+log2(1+u(:,i1)'*HHH'*HHH*u(:,i1)/temp4);
end
capacity=tempC;