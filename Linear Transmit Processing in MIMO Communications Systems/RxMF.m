function G = RxMF(P,H,R_s,R_n,alpha)
% Receive Matched Filter
% Decoder
% G_mf = alpha*R_s*P'*H'*R_eta^(-1)

G = alpha*R_s*P'*H'/R_n;
