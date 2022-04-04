function [count,capacity] = coherence_receive_SingleStream(N,K_all,K,Mi,KK,Tc,s,x,H,sigma2,w,u);

count=0;
capacity=0;

% single stream
s_now=zeros(K,2*Tc);    % Attention: the number of columns is influenced by the modulation method
x_now=zeros(K,Tc);
H_now=zeros(K*Mi,N);
w_now=zeros(K*Mi,Tc);

for i1=1:K,
    s_now(i1,:)=s(KK(1,i1),:);
    x_now(i1,:)=x(KK(1,i1),:);
    H_now((i1-1)*Mi+1:i1*Mi,:)=H((KK(1,i1)-1)*Mi+1:KK(1,i1)*Mi,:);
    w_now((i1-1)*Mi+1:i1*Mi,:)=w((KK(1,i1)-1)*Mi+1:KK(1,i1)*Mi,:);
end

y=zeros(K*Mi,Tc);
y=H_now*u*x_now+w_now;
% matched filter
r=zeros(K,K*Mi);
temp2=zeros(1,Mi);
for i1=1:K,
    HHH=H_now((i1-1)*Mi+1:i1*Mi,:);
    temp2=(HHH*u(:,i1))';
    temp2=temp2/(sqrt(trace(temp2*temp2')));
    r(i1,(i1-1)*Mi+1:i1*Mi)=temp2;    
end

yy=zeros(K,Tc);
yy=r*y;

% If the modulation method is different, the detection strategy should also be changed.
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
    temp4=Mi*sigma2;
    HHH=H_now((i1-1)*Mi+1:i1*Mi,:);
    for i2=1:K,
        if i2~=i1,
            temp4=temp4+u(:,i2)'*HHH'*HHH*u(:,i2);
        end
    end
    tempC=tempC+log2(1+u(:,i1)'*HHH'*HHH*u(:,i1)/temp4);
end
capacity=tempC;