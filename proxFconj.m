function [out1, out2, out3] = proxFconj(z1, z2, z3, b, sigma, eta)

%We know that the form of prox of F2* is the same as prox of F3* so we
%reuse the function
    out1 = proxSigF1conj(z1, b, sigma);
    out2 = proxSigfF2conj(z2, eta);
    out3 = proxSigfF2conj(z3, eta);
    
end

function out = proxSigF1conj(z1, b, sigma)
    out = (z1 - sigma*b)/(1+sigma);
end

function out = proxSigfF2conj(z2, eta)
    out = min(z2, eta);
    out = max(out, -eta);
end

