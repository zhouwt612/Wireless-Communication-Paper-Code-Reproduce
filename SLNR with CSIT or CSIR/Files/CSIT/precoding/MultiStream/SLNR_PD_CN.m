function [u, channel_norm] = SLNR_PD_CN(N,Mi,S,K,H,sigma2)

u=zeros(N,K*S);
channel_norm=zeros(1,K);

temp1=zeros(1,N);
Hi=zeros(Mi,N);
Hk=zeros(Mi,N);

for i1=1:K,
    Hi=H((i1-1)*Mi+1:i1*Mi,:);
    C=Hi'*Hi;
    channel_norm(1,i1)=real(trace(C));

    D=Mi*sigma2*eye(N);
    for i2=1:K,
        if i2~=i1,
            Hk=H((i2-1)*Mi+1:i2*Mi,:);
            D=D+Hk'*Hk;
        end
    end
    
    [eig_vec,eig_val]=eig(C,D);
    for i2=1:N,
        temp1(1,i2)=eig_val(i2,i2);
    end
    
    for i2=1:S,
        [max_val,max_pos]=max(temp1);
        v=eig_vec(:,max_pos);
        u(:,(i1-1)*S+i2)=v/(sqrt((v')*v));
        temp1(1,max_pos)=-1;
    end
end