clear; close all;


n = 20;
m = 10;
%Ablock = sparse(eye(n));
%Ablock(:,end) = ones(n,1);
%Ablock(end,:) = 0;

Ablock = triu(ones(n,n));
Ablock(end,:) = 0;

A = zeros(n*m);

s = 1:n:n*m;
e = n:n:n*m;

for i = 1:m
    A(s(i):e(i), s(i):e(i)) = Ablock;
end

figure, spy(A)
set(gca, 'XtickLabel', [], 'YtickLabel', [], ...
  'xtick', [], 'ytick', [] )
xlabel([])
