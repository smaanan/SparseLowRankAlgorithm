function [sbc, fpe, rnml, aicc] = arord_mod(data, R, m, mcor, ne, pmin, pmax)
imax     = pmax-pmin+1;
sbc      = zeros(1, imax);
fpe      = zeros(1, imax);
rnml     = zeros(1, imax);
aicc     = zeros(1, imax);
logdp    = zeros(1, imax);
np       = zeros(1, imax);
np(imax) = m*pmax;

R22      = R(np(imax)+1 : np(imax)+m, np(imax)+1 : np(imax)+m);
Delta    = R22'*R22;
invR22   = inv(R22);
Mp       = invR22*invR22';

logdp(imax) = 2.*log(abs(prod(diag(R22))));
i = imax;
Y = data(pmax+1:end,:);
trYY = trace(Y'*Y);

for p = pmax:-1:pmin,
    np(i) = m*p;
    if p < pmax
        Rp = R(np(i)+1:np(i)+m, np(imax)+1:np(imax)+m);
        Delta_prime = Delta + Rp'*Rp;
        L  = chol(eye(m) + Rp*Mp*Rp')';
        N  = L \ Rp*Mp;
        Mp = Mp - N'*N;
        logdp(i) = logdp(i+1) + 2.* log(abs(prod(diag(L))));
    end
    sbc(i) = logdp(i)/m - log(ne) * (ne-np(i))/ne;
    fpe(i) = logdp(i)/m - log(ne*(ne-np(i))/(ne+np(i)));
    T      = ne;
    ell    = np(i);
    term   = zeros(4,1);
    term(1)= ((T-ell-m+1)/2)*(logdp(i)-m*log(T));
    term(2)= -sum(gammaln((T-ell:-1:T-ell+1-m)/2));
    term(3)= -gammaln(ell*m/2);
    F      = (trYY - trace(Delta))/T;
    if p < pmax
        Delta = Delta_prime;
    end
    term(4)= ((ell*m)/2)*log(F);
    rnml(i)= sum(term);
    b      = T/(T-(p*m+m+1));
    aicc(i)= T*logdp(i) + 2*b*(p*m^2+m*(m+1)/2);
    i      = i-1;
end
end