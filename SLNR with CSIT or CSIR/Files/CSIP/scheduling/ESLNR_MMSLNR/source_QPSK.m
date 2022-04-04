function [s,x] = source_QPSK(K,S,Tc)

s=zeros(K*S,2*Tc);
x=zeros(K*S,Tc);

for i1=1:K*S,
    for i2=1:2*Tc,
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
                x_real=(sqrt(2))/2;
            case 1,
                x_real=-(sqrt(2))/2;
        end
        switch s(i1,i2+Tc),
            case 0,
                x_imag=(sqrt(2))/2;
            case 1,
                x_imag=-(sqrt(2))/2;
        end
        x(i1,i2)=complex(x_real,x_imag);
    end
end