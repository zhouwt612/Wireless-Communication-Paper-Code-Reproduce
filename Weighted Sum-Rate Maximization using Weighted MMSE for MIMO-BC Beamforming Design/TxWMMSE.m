function Bwmmse = TxWMMSE(H,A,W,Etx)
% H (QKxP): Channel realization
% A (QKxQK): MMSE receive filter
% W (QKxQK): MSE weights
% Etx: Power constraint
% B (BxQK): Receive filter
B_tilde = inv(H'*A'*W*A*H + (trace(W*A*A')/Etx)*eye(size(H'*H)))*H'*A'*W;
b = sqrt(Etx/trace(B_tilde*B_tilde'));
Bwmmse = b*B_tilde;