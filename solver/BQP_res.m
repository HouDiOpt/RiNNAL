function [pfeas,dfeas,comp,pdgap,pfval,dfval,V,d] = BQP_res(X,Y,lambda,lam,Q,c,par)
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
e1t = zeros(1,r); e1t(1) = 1;
R1 = [e1t;R];
dRR = diag(RR)-R(:,1);
ARR = (A*R)*R';
QR = Q*R;
%% recover subproblem dual variable: lambda_hat
muhat  = lam(1:l);
lamhat = lam(l+1:end);
lammat = reshape(lamhat,m,r); 
%% recover dual variables
JA = eye(n)-Ainv*A; JA = (JA+JA')/2;
muB = zeros(n,1); muB(id) = muhat;
L = Q-diag(muB)-lambda(2:end,2:end);
h = 2*c+muB-2*lambda(2:end,1);
S = [0,-lambda(1,2:end);-lambda(2:end,1),L];
S21 = c + muB/2 + L*(Ainv*b);
alphaA = b'*lammat(:,1)/2 + h'*R(:,1)/2 - b'*Ainv'*h/2 - (h'/2 + b'*Ainv'*L)*(Ainv*b);
S = S + [-alphaA,S21';S21,zeros(n,n)];
S = [1,zeros(1,n);zeros(n,1),JA]*S*[1,zeros(1,n);zeros(n,1),JA];
S = (S+S')/2;
nS = norm(S,'fro');
nY = norm(Y,'fro');
%% escaping S direction
[V,d] = eig(S,'vector');
[d,ind] = sort(d);
V = V(ind,ind);
V = JA*[-R(:,1),eye(n)]*V;
d = d.*(d<0);
alpha = -S(1,1);
lam1 = Ainv'*(h+L*(Ainv*b));
%% compute KKT
pfval = Produc(QR,R)+2*c'*R(:,1);
dfval = lam1'*b + alpha;
pfeas = max(sqrt(norm(A*R(:,1)-b,'fro')^2 + norm(ARR-b*R(:,1)','fro')^2 ...
        + norm(dRR(id))^2) / (1+norm(b)), norm(RR1-Y,'fro') / (1+RR1n+nY));
dfeas = norm(d,'fro')/(1+nS); 
comp = max(abs(Produc(RR1,S))/(1+norm(RR1,'fro')+nS),...
    abs(Produc(Y,lambda))/(1+nY+norm(lambda,'fro')));
pdgap = abs(pfval-dfval)/(1+abs(pfval)+abs(dfval));
end