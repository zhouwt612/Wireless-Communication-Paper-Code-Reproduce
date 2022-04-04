function [u, r_temp, gama] = BD_PD(N,Mi,S,K,H,sigma2)

u=zeros(N,K*S);
r_temp=zeros(K*S,Mi);
gama=zeros(1,K*S);
for i1=1:K,
    H_others=zeros((K-1)*Mi,N);
    i3=1;
    for i2=1:K,
        if i2~=i1,
            H_others((i3-1)*Mi+1:i3*Mi,:)=H((i2-1)*Mi+1:i2*Mi,:);
            i3=i3+1;
        end
    end
    
    U=zeros((K-1)*Mi,(K-1)*Mi);
    D=zeros((K-1)*Mi,N);
    V=zeros(N,N);
    [U,D,V]=svd(H_others);
    num=(K-1)*Mi+1;
    for i2=(K-1)*Mi:(-1):1,
        if D(i2,i2)==0,
            num=i2;
        end
    end
    T=V(:,num:N);
    
    HHH=H((i1-1)*Mi+1:i1*Mi,:);
    U=zeros(Mi,Mi);
    D=zeros(Mi,N-num+1);
    V=zeros(N-num+1,N-num+1);
    [U,D,V]=svd(HHH*T);
    
    u(:,(i1-1)*S+1:i1*S)=T*V(:,1:S);
    r_temp((i1-1)*S+1:i1*S,:)=U(:,1:S)';
    for i2=1:S,
        gama(1,(i1-1)*S+i2)=D(i2,i2);
    end
end

for i1=1:K,
    u(:,i1)=u(:,i1)/(sqrt((u(:,i1)')*u(:,i1)));
end