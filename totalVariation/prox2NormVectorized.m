function out = prox2NormVectorized(xHat,t)

[numRows,numCols,check] = size(xHat);
coeffs = zeros(numRows,numCols);
nrms = sqrt(xHat(:,:,1).^2 + xHat(:,:,2).^2);
fnd = find(nrms >= t);
coeffs(fnd) = 1 - t(fnd)./nrms(fnd);
coeffs = cat(3,coeffs,coeffs);
out = coeffs.*xHat;

% nrm = norm(xHat(:));
% if t <= nrm
%     out = (1-t/nrm)*xHat;
% else
%     out = zeros(size(xHat));
% end

end