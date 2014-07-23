function div = computeDiv_neumannBCs(u)

u1 = u(:,:,1);
u2 = u(:,:,2);

[numRows,numCols] = size(u1);

du1dx = u1(2:end-1,:) - u1(1:end-2,:);
du1dx = [u1(1,:) ; du1dx; -u1(end-1,:)];

du2dy = u2(:,2:end-1) - u2(:,1:end-2);
du2dy = [u2(:,1), du2dy,-u2(:,end-1)];

div = du1dx + du2dy;

end