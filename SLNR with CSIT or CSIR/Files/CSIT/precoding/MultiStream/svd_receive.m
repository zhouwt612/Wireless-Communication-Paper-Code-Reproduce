function [count,capacity] = svd_receive(N,Mi,S,K,KK,Tc,s,x,H,sigma2,w,u,p,r_temp)

count=0;
capacity=0;

% single stream
s_now=zeros(K*S,2*Tc);      % If the modulation method is different, the amount of columns may be changed.
x_now=zeros(K*S,Tc);
H_now=zeros(K*Mi,N);
w_now=zeros(K*Mi,Tc);

for i1=1:K,
    s_now((i1-1)*S+1:i1*S,:)=s((KK(1,i1)-1)*S+1:KK(1,i1)*S,:);
    x_now((i1-1)*S+1:i1*S,:)=x((KK(1,i1)-1)*S+1:KK(1,i1)*S,:);
    H_now((i1-1)*Mi+1:i1*Mi,:)=H((KK(1,i1)-1)*Mi+1:KK(1,i1)*Mi,:);
    w_now((i1-1)*Mi+1:i1*Mi,:)=w((KK(1,i1)-1)*Mi+1:KK(1,i1)*Mi,:);
end

y=H_now*u*sqrt(diag(p))*x_now+w_now;
% matched filter
r=zeros(K*S,K*Mi);
for i1=1:K,
    r_temp((i1-1)*S+1:i1*S,:)=r_temp((i1-1)*S+1:i1*S,:)/(sqrt(trace(r_temp((i1-1)*S+1:i1*S,:)*r_temp((i1-1)*S+1:i1*S,:)')));
    r((i1-1)*S+1:i1*S,(i1-1)*Mi+1:i1*Mi)=r_temp((i1-1)*S+1:i1*S,:);
end

yy=zeros(K*S,Tc);
yy=r*y;
        
ss=zeros(K*S,2*Tc);
for i1=1:K*S,
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

for i1=1:K*S,
    for i2=1:2*Tc,
        if s_now(i1,i2)~=ss(i1,i2),
            count=count+1;
        end
    end
end

tempC=0;
for i1=1:K,
    temp4=Mi*sigma2*eye(S);
    HHH=H_now((i1-1)*Mi+1:i1*Mi,:);
    for i2=1:K,
        if i2~=i1,
            temp4=temp4+diag(p(1,(i2-1)*S+1:i2*S))*u(:,(i2-1)*S+1:i2*S)'*HHH'*HHH*u(:,(i2-1)*S+1:i2*S);
        end
    end
    tempC=tempC+log2(det(eye(S)+diag(p(1,(i1-1)*S+1:i1*S))*u(:,(i1-1)*S+1:i1*S)'*HHH'*HHH*u(:,(i1-1)*S+1:i1*S)*inv(temp4)));
end
capacity=tempC;