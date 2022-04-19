function [F,s_proj,gamma] = ACalgorithm(H,C)
% H: channel realization (size: M,N,K)
% C: codebook (size: M, L=2^B, K)
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of users
% L: the number of codewords
% F: quantized channel
% s_proj: the projection of quantization vector onto span H
% gamma: weighting vector for all users

spanH = Space_SpanH(H);
F = quantizedchannel_AC(spanH,C);
s_proj = quantizedvector_proj(spanH,F);
gamma = weighting_vector(H,s_proj);