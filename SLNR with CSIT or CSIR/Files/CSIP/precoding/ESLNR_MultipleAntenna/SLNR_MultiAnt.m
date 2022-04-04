function [u, ant_num] = SLNR_MultiAnt(N,Mi,K,KK,H,sigma2)

ant_num=ones(1,K);

% 根据调度结果，读入被调度者的信道系数
H_now=zeros(K,N);
for i1=1:K,
    H_now(i1,:)=H((KK(1,i1)-1)*Mi+1,:);
end

% pre-coding
u=zeros(N,K);      % single stream
Hi=zeros(1,N);
Hk=zeros(1,N);
temp1=zeros(1,N);
SLNR_temp=zeros(1,K);
for i1=1:K,
    Hi=H_now(i1,:);
    C=Hi'*Hi;

    D=sigma2*eye(N);
    for i2=1:K,
        if i2~=i1,
            Hk=H_now(i2,:);
            D=D+Hk'*Hk;
        end
    end
    
    [eig_vec,eig_value]=eig(C,D);
    for n=1:N,
        temp1(1,n)=eig_value(n,n);
    end
    [max_val,max_pos]=max(temp1);
    SLNR_temp(1,i1)=max_val;
    v=eig_vec(:,max_pos);
    u(:,i1)=v/(sqrt((v')*v));
end

while 0==0,
    [min_val, min_pos]=min(SLNR_temp);

    ant_num_temp=zeros(1,K);
    ant_num_temp=ant_num;
    
    if ant_num(1,min_pos)<Mi,
        ant_num(1,min_pos)=ant_num(1,min_pos)+1;
    else
        break;
    end
    
    H_now=zeros(sum(ant_num),N);
    for i1=1:K,
        for i2=1:ant_num(1,i1),
            temp1=0;
            if i1~=1,
                temp1=sum(ant_num(1,1:(i1-1)));
            end
            H_now(temp1+1:temp1+ant_num(1,i1),:)=H((KK(1,i1)-1)*Mi+1:(KK(1,i1)-1)*Mi+ant_num(1,i1),:);
        end
    end

    % pre-coding
    u_temp=zeros(N,K);      % single stream
    temp1=zeros(1,N);
    SLNR_new=zeros(1,K);
    for i1=1:K,
        Hi=zeros(ant_num(1,i1),N);
        temp2=0;
        if i1~=1,
            temp2=sum(ant_num(1,1:(i1-1)));
        end
        Hi=H_now(temp2+1:temp2+ant_num(1,i1),:);

        C=Hi'*Hi;

        D=ant_num(1,i1)*sigma2*eye(N);
        for i2=1:K,
            if i2~=i1,
                Hk=zeros(ant_num(1,i2),N);
                temp2=0;
                if i2~=1,
                    temp2=sum(ant_num(1,1:(i2-1)));
                end
                Hk=H_now(temp2+1:temp2+ant_num(1,i2),:);
                D=D+Hk'*Hk;
            end
        end

        [eig_vec,eig_val]=eig(C,D);
        for n=1:N,
            temp1(1,n)=eig_val(n,n);
        end
        [max_val,max_pos]=max(temp1);
        SLNR_new(1,i1)=max_val;
        v=eig_vec(:,max_pos);
        u_temp(:,i1)=v/(sqrt((v')*v));
    end
    
    if min(real(SLNR_new))<real(min_val),
        ant_num=ant_num_temp;
        break;
    else
        SLNR_temp=SLNR_new;
        u=u_temp;
    end
    
end