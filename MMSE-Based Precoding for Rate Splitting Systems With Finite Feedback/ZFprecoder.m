function P = ZFprecoder(H,pow)
Pb = H';
P = sqrt(pow/trace(Pb*Pb'))*Pb;