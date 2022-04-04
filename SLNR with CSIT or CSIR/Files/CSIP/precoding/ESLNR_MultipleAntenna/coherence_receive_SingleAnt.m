function [count,capacity] = coherence_receive_SingleAnt(MT,R,K_all,K,KK,Tc,s,x,H,sigma2,w,u,index);

count=0;
capacity=0;

% single stream
s_now=zeros(K,Tc);    % Attention: the number of columns is influenced by the modulation method
x_now=zeros(K,Tc);
H_now=zeros(K,MT);
w_now=zeros(K,Tc);

for i1=1:K,
    s_now(i1,:)=s(KK(1,i1),:);
    x_now(i1,:)=x(KK(1,i1),:);
    H_now(i1,:)=H((KK(1,i1)-1)*R+index(1,(i1-1)*R+1),:);
    w_now(i1,:)=w((KK(1,i1)-1)*R+index(1,(i1-1)*R+1),:);
end

y=zeros(K,Tc);
y=H_now*u*x_now+w_now;
% matched filter
r=zeros(K,K);
temp2=zeros(1,1);
for i1=1:K,
    HHH=H_now(i1,:);
    temp2=(HHH*u(:,i1))';
    temp2=temp2/(sqrt(trace(temp2*temp2')));
    r(i1,i1)=temp2;    
end

yy=zeros(K,Tc);
yy=r*y;

% If the modulation method is different, the detection strategy should also be changed.
ss=zeros(K,Tc);
for i1=1:K,
    for i2=1:Tc,
        if real(yy(i1,i2))>=0,
            ss(i1,i2)=0;
        else
            ss(i1,i2)=1;
        end
        
%         if imag(yy(i1,i2))>=0,
%             ss(i1,i2+Tc)=0;
%         else
%             ss(i1,i2+Tc)=1;
%         end

    end
end

for i1=1:K,
    for i2=1:Tc,
        if s_now(i1,i2)~=ss(i1,i2),
            count=count+1;
        end
    end
end

tempC=0;
for i1=1:K,
    temp4=sigma2;
    HHH=H_now(i1,:);
    for i2=1:K,
        if i2~=i1,
            temp4=temp4+u(:,i2)'*HHH'*HHH*u(:,i2);
        end
    end
    tempC=tempC+log2(1+u(:,i1)'*HHH'*HHH*u(:,i1)/temp4);
end
capacity=tempC;