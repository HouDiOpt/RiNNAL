function [fval,X,info] = RiNNAL_diag(Q,c,A,b,bidx,r,par)
% Problem:
%            min  <Q,X>+2<c,x>  s.t.
%            Ax = b            (recovered)
%            AX = bx'          (recovered)
%            diag(X)_B = x_B   (lam)
%            X  >= 0, X PSD    (lambda)
%            X_E = 0
%
% Algorithm: ALM + BM + Rie GD with BB stepsize

%% User defiend options
tstart = clock;
rng('default');
verbose = 1;
vergap = 1;
tol = 1e-6;
maxiter = 1e5;
beta = 1;
ra = 100;
scale = 1;
ramin = 1e-2;
betamin = 1e-6;
betamax = 1e2;
tolrrate = 1;
subscale = 1;
uptol = 1e-8;
timelimit = 60*60;
if isfield(par,'verbose'); verbose = par.verbose; end
if isfield(par,'vergap'); vergap = par.vergap; end
if isfield(par,'tol'); tol = par.tol; end
if isfield(par,'beta'); beta = par.beta; end
if isfield(par,'betamin'); betamin = par.betamin; end
if isfield(par,'betamax'); betamax = par.betamax; end
if isfield(par,'ra'); ra = par.ra; end
if isfield(par,'subscale'); subscale = par.subscale; end
if isfield(par,'ramin'); ramin = par.ramin; end
if isfield(par,'scale'); scale = par.scale; end
if isfield(par,'tolrrate'); tolrrate = par.tolrrate; end
if isfield(par,'uptol'); uptol = par.uptol; end
if isfield(par,'timelimit'); timelimit = par.timelimit; end
%% Scale the problem (Q,c)
if (scale==1)
    nCorg = sqrt(norm(Q,'fro')^2+norm(c,'fro')^2);
    nC = max(1,nCorg);
elseif (scale==0)
    nC = 1;
end
Qs = Q/nC;
cs = c/nC;
%% Construct projection and retraction operator
parPR = BQP_parameter_diag(A,b,bidx,par);
parPR.tolrrate = tolrrate;
if isfield(par,'incmax'); parPR.incmax = par.incmax; end
Rie_proj = @(C,X,parPR)BQP_proj_diag(C,X,parPR);
Rie_retrac = @(X,H,s,tol,parPR)BQP_retrac_diag(X,H,s,tol,parPR);
%% Initialization
n = size(Q,1);
m = size(A,1);
parPR.checkres = 1;
if isfield(par,'x0')
    R0 = 1e-5*randn(n,r);
    % R0 = zeros(n,r);
    R0(:,1) = 2*par.x0-1;
    X.R = R0;
    X = Rie_retrac(X,0,0,1e-5,parPR);
elseif isfield(par,'R')
    X.R = par.R;
else
    X.R = randn(n,r);
    X = Rie_retrac(X,0,0,1e-5,parPR);
end
parPR.checkres = 0;
%% Parameter for BB
parBB.tolg = 1e16;
parBB.verbose = 0;
parBB.safegard = 1;
parBB.eigratio = 1e2;
parBB.maxiter = 1e4;
parBB.alp = 1;
parBB.adstop = 1;
parBB.checkrank = 0;
parBB.connec = 1;
parBB.plotyes = 0;
parBB.scaler = nC;
parBB.uptol  = uptol;
if isfield(par,'BBmaxiter'); parBB.maxiter = par.BBmaxiter; end
%% Print header
if verbose
    fprintf('\n ite |   fval     |  beta  |  pfeas    dfeas    comp     pdgap |   ra   |gradnorm|Bite|  time  |rank|');
end
%% Initialization
lambda1 = zeros(m,1);     % dual var for Ax=b
lambda2 = zeros(m,n);     % dual var for AX=bx'
lambda3 = zeros(n+1,n+1); % dual var for Y>=0
if isfield(par,'lambda1'); lambda1 = par.lambda1; end
if isfield(par,'lambda2'); lambda1 = par.lambda2; end
if isfield(par,'lambda3'); lambda3 = par.lambda3; end
lambda1 = lambda1/nC;
lambda2 = lambda2/nC;
lambda3 = lambda3/nC;
pfeas  = inf;
BBiterTotal = 0;
%% Main loop
for i = 1:maxiter
    %% pre
    pfeas0  = pfeas;
    f_g = @(X)f_and_g(X,Qs,cs,lambda1,lambda2,lambda3,beta,parPR);
    %% inner loop
    %% update primal variable
    [X,lam,iter,~,~,gradnorm] = Rie_BB_Sub(f_g,Rie_proj,Rie_retrac,X,ra,parBB,parPR);
    BBiterTotal = BBiterTotal + iter;
    %% update  dual  variable
    R = X.R;
    R(:,1) = R(:,1)+1;
    R = R/2;
    RR = R*R';
    RR1 = [1,R(:,1)';R(:,1),RR];
    Z = RR1-lambda3/beta;
    be1t = zeros(m,size(R,2)); be1t(:,1) = b;
    lambda1 = lambda1-beta*(Amap(R(:,1),parPR,0)-b);
    lambda2 = lambda2-beta*(Amap(R,parPR,0) - be1t)*R';
    lambda3 = -beta*Z.*(Z<=0);
    %% compute residue
    [pfeas,dfeas,comp,pdgap,pfval,dfval,V,d] = BQP_res_diag(X,Z.*(Z>=0),lambda1*nC,lambda2*nC,lambda3*nC,2*lam*nC,Q,c,parPR);
    subfeas  = subscale*max(dfeas,comp);
    res = max([pfeas,dfeas,comp]);
    ttime = etime(clock,tstart);
    if verbose && mod(i,vergap) == 0
        fprintf('\n %3d |%+6.5e|%3.2e|%3.2e %3.2e %3.2e %3.2e|%3.2e|%3.2e|%4d|%3.2e|',i,pfval,beta,pfeas,dfeas,comp,pdgap,ra,gradnorm,iter,ttime);
        fprintf('%3d |', size(X.R,2));
    end
    %% check stopping criterion
    if res < tol
        fprintf('\n residue = %3.2e < tol = %3.2e, Stop! \n',res,tol);
        break;
    end
    if ttime > timelimit
        fprintf('\n Reached time limit! residue = %3.2e > tol = %3.2e, Stop! \n',res,tol);
        break;
    end
    %% update ra
    if subfeas < 0.01*pfeas
        ra = ra*2;
    elseif subfeas > 0.8*pfeas
        ra = max(ramin,ra/2);
        [X.R,~] = rankinc_diag(X.R,V,d/nC,f_g,Rie_proj,Rie_retrac,parPR);
    end
    %% update beta
    if pfeas > 0.9*pfeas0 || pfeas > 10*subfeas
        beta = min(betamax,beta*1.25);
    elseif pfeas<subfeas*0.1 %pfeas < 0.2*pfeas0 && ra <= ramin && iter>1000
        beta = max(betamin,beta/1.25);
    end
    %% reduce rank
    if size(X.R,2)>1
        [X,~] = rankred(X,Rie_retrac,parPR);
    end
end
%% Record results
fval(1) = pfval;
fval(2) = dfval;
info.pfeas = pfeas;
info.dfeas = dfeas;
info.comp  = comp;
info.pdgap = pdgap;
info.ttime = ttime;
info.ALMite = i;
info.BBite  = BBiterTotal;
info.lambda1 = lambda1*nC;
info.lambda2 = lambda2*nC;
info.lambda3 = lambda3*nC;
info.lam = 2*lam*nC;
%% Print
fprintf('\n pobj  = %6.7e',fval(1));
fprintf('\n dobj  = %6.7e',fval(2));
fprintf('\n rank  = %3d',size(X.R,2));
fprintf('\n pfeas = %3.2e',pfeas);
fprintf('\n dfeas = %3.2e',dfeas);
fprintf('\n comp  = %3.2e',comp);
fprintf('\n pdgap = %3.2e',pdgap);
fprintf('\n time  = %3.2e\n\n',ttime);
end

%% Functions
function [f,eugf,pfeas,RR1] = f_and_g(X,Q,c,lambda1,lambda2,lambda3,beta,par)
b = par.b;
m = size(b,1);
R = X.R;
R(:,1) = R(:,1)+1;
R = R/2;
r = size(R,2);
e1t = zeros(1,r); e1t(1) = 1;
be1t = zeros(m,r); be1t(:,1) = b;
R1 = [e1t;R];
RR1 = R1*R1';
Z = RR1-lambda3/beta;
Y = Z.*(Z>=0);
QR = Q*R;
pres1 = Amap(R(:,1),par,0)-b;
ARmbe1t = Amap(R,par,0) - be1t;
pres2 = ARmbe1t*R';
tmp1 = pres1-lambda1/beta;
tmp2 = pres2-lambda2/beta;
tmp3 = Z.*(Z<=0);
%%
f = Produc(QR,R) + 2*c'*R(:,1) + (beta/2)*((norm(tmp1,'fro')^2+norm(tmp2,'fro')^2+norm(tmp3,'fro')^2));
eugf = QR + beta*tmp3(2:end,:)*R1 + (beta/2)*Atmap(tmp1*e1t,par,0) + (beta/2)*(Atmap(tmp2*R,par,0)+tmp2'*ARmbe1t);
eugf(:,1) = eugf(:,1)+c; % no need eugf*2 because here R=(R+...)/2
pfeas = max(norm(RR1-Y,'fro')/(1+norm(RR1,'fro')), sqrt(norm(pres1,'fro')^2+norm(pres2,'fro')^2)/(1+norm(b)));
end