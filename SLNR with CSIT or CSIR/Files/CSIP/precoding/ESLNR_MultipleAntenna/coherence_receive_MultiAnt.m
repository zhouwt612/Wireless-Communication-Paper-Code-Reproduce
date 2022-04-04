function [count,capacity] = coherence_receive_MultiAnt(MT,R,K_all,K,KK,Tc,s,x,H,sigma2,w,u,index,ant_num);

count=0;
capacity=0;

% single stream
s_now=zeros(K,Tc);    % Attention: the number of columns is influenced by the modulation method
x_now=zeros(K,Tc);
H_now=zeros(sum(ant_num),MT);
w_now=zeros(sum(ant_num),Tc);

for i1=1:K,
    s_now(i1,:)=s(KK(1,i1),:);
    x_now(i1,:)=x(KK(1,i1),:);
    temp1=0;
    if i1~=1,
        temp1=sum(ant_num(1,1:(i1-1)));
    end
    for i2=1:ant_num(1,i1),
        H_now(temp1+i2,:)=H((KK(1,i1)-1)*R+index(1,(i1-1)*R+i2),:);
        w_now(temp1+i2,:)=w((KK(1,i1)-1)*R+index(1,(i1-1)*R+i2),:);
    end
end

y=zeros(sum(ant_num),Tc);
y=H_now*u*x_now+w_now;
% matched filter
r=zeros(K,sum(ant_num));
for i1=1:K,
    temp1=0;
    if i1~=1,
        temp1=sum(ant_num(1,1:(i1-1)));
    end
    HHH=zeros(ant_num(1,i1),MT);
    HHH=H_now(temp1+1:temp1+ant_num(1,i1),:);
    temp2=zeros(1,ant_num(1,i1));
    temp2=(HHH*u(:,i1))';
    temp2=temp2/(sqrt(trace(temp2*temp2')));
    r(i1,temp1+1:temp1+ant_num(1,i1))=temp2;
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
    temp4=ant_num(1,i1)*sigma2;
    temp1=0;
    if i1~=1,
        temp1=sum(ant_num(1,1:(i1-1)));
    end
    HHH=zeros(ant_num(1,i1),MT);
    HHH=H_now(temp1+1:temp1+ant_num(1,i1),:);
    for i2=1:K,
        if i2~=i1,
            temp4=temp4+u(:,i2)'*HHH'*HHH*u(:,i2);
        end
    end
    tempC=tempC+log2(1+u(:,i1)'*HHH'*HHH*u(:,i1)/temp4);
end
capacity=tempC;