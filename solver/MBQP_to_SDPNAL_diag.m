function [blk,At,C,bb]=MBQP_to_SDPNAL_diag(Q,c,A,b,bidx)
%% pars
[m,~] = size(A);
n = size(Q,1);
nc = 2*m+length(bidx)+1;
%% vars
blk{1,1} = 's'; blk{1,2} = n+1;
Acell = cell(1,nc);
for i = 1:m
    Acell{i} = [0,A(i,:)/2;A(i,:)'/2,sparse(n,n)];
    Acell{m+i} = [0,sparse(1,n);sparse(n,1),A(i,:)'*A(i,:)];
end
for i = 1:length(bidx)
    ei = zeros(n,1); ei(bidx(i)) = 1;
    eii = zeros(n,n); eii(bidx(i),bidx(i)) = 1;
    Acell{2*m+i} = [0,ei'/2;ei/2,-eii];
end
Acell{nc} = [1,sparse(1,n);sparse(n,1),sparse(n,n)];
At = svec(blk,Acell,1);
bb = zeros(nc,1);
bb(1:m,1)   = b;
bb(m+1:2*m,1) = b.^2;
bb(2*m+1:end-1) = 0;
bb(end) = 1;
C = {[0,c';c,Q]};
end