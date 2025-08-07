function [blk,At,C,bb]=MBQP_to_SDPNAL(Q,c,A,b,bidx,zid)
%% pars
[m,~] = size(A);
n = size(Q,1);
nc = m+m*n+length(bidx)+1+length(zid)/2;
%% vars
blk{1,1} = 's'; blk{1,2} = n+1;
Acell = cell(1,nc);
%% Ax=b
for i = 1:m
    Acell{i} = [0,A(i,:)/2;A(i,:)'/2,sparse(n,n)];
end
%% AX=bx'
for i = 1:m
    for j = 1:n
        bi = sparse(j,1,-b(i)/2,n,1);
        AX = sparse(j,1:n,A(i,:),n,n);
        AX = (AX+AX')/2;
        Acell{m+(i-1)*n+j} = [0,bi';bi,AX];
    end
end
%% diag_B(X)=x_B
for i = 1:length(bidx)
    ei = sparse(bidx(i),1,1,n,1); 
    eii = sparse(bidx(i),bidx(i),1,n,n);
    Acell{m+m*n+i} = sparse([0,ei'/2;ei/2,-eii]);
end
%% Y_{11} = 1
Acell{m+m*n+length(bidx)+1} = [1,sparse(1,n);sparse(n,1),sparse(n,n)];
%% X_{ij}=0, ij in zid
id0 = m+m*n+length(bidx)+1;
M = zeros(n+1,n+1);
M(zid) = 1;
M = triu(M);
[rowid,colid] = find(M);
for i = 1:length(zid)/2
    Acell{id0+i} = sparse([rowid(i),colid(i)],[colid(i),rowid(i)],[1/2,1/2],n+1,n+1);
end
At = svec(blk,Acell,1);
bb = zeros(nc,1);
bb(1:m,1) = b;
bb(m+m*n+length(bidx)+1) = 1;
C = {[0,c';c,Q]};
end