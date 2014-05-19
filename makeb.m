function b = makeb(I)

  nDimI = ndims(I);

  switch nDimI
    case 1
      m = length(I);
      ut = triu(ones(m-1, m-1));
      b = makeb1D(I, ut);
    case 2
      [m, ~] = size(I);
      ut = triu(ones(m-1, m-1));
      b = makeb2D(I, ut);
    case 3
      [m, ~, ~] = size(I);
      ut = triu(ones(m-1, m-1));
      b = makeb3D(I, ut);
    otherwise
        error('Improper size of I');
  end

    
end

function b = makeb1D(I, ut)
  m = numel(I);
  b = zeros(m,1);
  b(1:m-1) = ut * I(1:m-1);
end

function b = makeb2D(I, ut)
    [m, n] = size(I);
    b = zeros(m, n);
    for i = 1:n
        b(:, i) = makeb1D(I(:,i), ut);
    end
end

function b = makeb3D(I, ut)
  [m, n, p] = size(I);
  b = zeros(m, n, p);
  for i = 1:n
    for j = 1:p
      b(:,i,j) = makeb1D(I(:,i,j), ut);
    end
  end
end