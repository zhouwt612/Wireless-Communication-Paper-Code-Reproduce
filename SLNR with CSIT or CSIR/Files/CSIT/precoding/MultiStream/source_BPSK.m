function [s,x] = source_BPSK(K,S,Tc)

s=zeros(K*S,Tc);
x=zeros(K*S,Tc);

for i1=1:K*S,
    for i2=1:Tc,
        temp=randn+0.5;
        if temp<=0.5
            s(i1,i2)=0;
        else
            s(i1,i2)=1;
        end
    end
end

for i1=1:K*S,
    for i2=1:Tc,
        switch s(i1,i2),
            case 0,
                x(i1,i2)=1*(1/sqrt(S));
            case 1,
                x(i1,i2)=-1*(1/sqrt(S));
        end
    end
end