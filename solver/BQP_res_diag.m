function [pfeas,dfeas,comp,pdgap,pfval,dfval,V2,d2] = BQP_res_diag(X,Y,lambda1,lambda2,lambda3,lam,Q,c,par)
l = par.l;
m = par.m;
A = par.A;
b = par.b;
id = par.bidx;
Ainv = par.Ainv;
R0 = X.R; 
[n,r] = size(R0);
R = R0;
R(:,1) = R(:,1)+1;
R = R/2;
RR = R*R';
RR1 = [1,R(:,1)';R(:,1),RR];
RR1n = norm(RR1,"fro");
dRR = diag(RR)-R(:,1);
ARR = (A*R)*R';
QR = Q*R;
%% recover dual variables v1 (in the form of C-At(y)-W)
mut = zeros(n,1); mut(id) = lam;
if size(lambda1,1) ~= 0
    S21 = c-Atmap(lambda1,par,0)/2+lambda2'*b/2+mut/2;
else
    S21 = c+lambda2'*b/2+mut/2;
end
S22 = Q-(Atmap(lambda2,par,0)+Atmap(lambda2,par,0)')/2-diag(mut);
S = [0,S21';S21,S22]-lambda3;
alpha = Produc(S,RR1);
S(1,1) = -alpha;
%%
[V,d] = eig(S,'vector');
[d,ind] = sort(d);
d = d.*(d<0);
%% escaping direction 
d2 = d;
V2 = V(ind,ind);
V2 = [-R(:,1),eye(n)]*V2;
%%
nS = norm(S,'fro');
nY = norm(Y,'fro');
%% compute KKT
pfval = Produc(QR,R)+2*c'*R(:,1);
dfval = lambda1'*b + alpha;
pfeas = max(sqrt(norm(A*R(:,1)-b,'fro')^2 + norm(ARR-b*R(:,1)','fro')^2 ...
        + norm(dRR(id))^2) / (1+norm(b)), norm(RR1-Y,'fro') / (1+RR1n+nY));
dfeas = norm(d,'fro')/(1+nS); 
comp = max(abs(Produc(RR1,S))/(1+norm(RR1,'fro')+nS),...
    abs(Produc(Y,lambda3))/(1+nY+norm(lambda3,'fro')));
pdgap = abs(pfval-dfval)/(1+abs(pfval)+abs(dfval));
end
