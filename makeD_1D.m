
function D = makeD_1D(M, mask, dz)

  tmp = circshift(mask, -1);
  tmp(end) = 1;
  masked = tmp .* mask;

  %Construct D
  D = diag(-1*masked);
  for i = 1:M-1
    D(i,i+1) = masked(i);
  end   
  D(end) = 0;
  D = D ./ dz;


%   D = -eye(M);
%   D(M,M) = 0;
%   for i=1:M-1
%     D(i,i+1) = 1;
%   end
% 
%   D = D ./ dz;

end
