function parPR = BQP_parameter_diag(A,b,bidx,par)
%%
[m,n] = size(A);
AAt = A*A';
At = A';
Atid = At(bidx,:);
AAtinv = inv(A*At);
diagAAtinv = 1./diag(AAt);
Ainv = A'/(A*A');%pinv(A);
Ainvid = Ainv(bidx,:);
Aid = A(:,bidx);
projA = eye(n)-Ainv*A; projA = (projA+projA')/2;
projAid = projA(bidx,bidx);
%% manifold type
if isempty(A)&&isempty(bidx)
    type = 3; % no constraints
elseif isempty(bidx)
    type = 2; % only affine
elseif isempty(A)
    type = 1; % only diag, but B=[n]
else
    type = 0; % both constraints 
end
%% binary type 
if length(bidx) == n
    id_full = 1;
else
    id_full = 0;
end
%% proj type
if isfield(par,'projtype')
    projtype = par.projtype; 
    parPR.projtype = projtype;
end
%% Amap
if isfield(par,'Amap1')
    parPR.Amaptype = 1;
    parPR.Amap1 = par.Amap1;
else
    parPR.Amaptype = 0;
end
%% Atmap
if isfield(par,'Atmap1')
    parPR.Atmaptype = 1;
    parPR.Atmap1 = par.Atmap1;
else
    parPR.Atmaptype = 0;
end
%% dimension
parPR.m = m;
parPR.n = n;
parPR.l = length(bidx);
%% assign
parPR.b  = b;
parPR.A  = A;
parPR.At = A';
parPR.Atid = Atid;
parPR.Ae = sum(A,2);
parPR.AAt = AAt;
parPR.Ainv = Ainv;%pinv(A)
parPR.Ainvid = Ainvid;
parPR.Aid = Aid;
parPR.Ainvt = Ainv';
parPR.AAtinv = AAtinv;
parPR.projA = projA;
parPR.projAid = projAid;
parPR.diagAAtinv = diagAAtinv;
parPR.bidx = bidx;
parPR.id_full = id_full;
parPR.bidxdiff = setdiff(1:n,bidx);
parPR.type = type;
end