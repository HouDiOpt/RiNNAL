function [fval,X,info] = RiNNAL(Q,c,A,b,bidx,r,zid,par)
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
ramax = 1e8;
betamin = 1e-6;
betamax = 1e2;
tolrrate = 1;
betascale = 1;
betascale2 = 1;
subscale = 1;
uptol = 1e-8;
maxtime = 3600;
redfrequence = 1;
if isfield(par,'verbose'); verbose = par.verbose; end
if isfield(par,'vergap'); vergap = par.vergap; end
if isfield(par,'tol'); tol = par.tol; end
if isfield(par,'beta'); beta = par.beta; end
if isfield(par,'betamin'); betamin = par.betamin; end
if isfield(par,'betamax'); betamax = par.betamax; end
if isfield(par,'ra'); ra = par.ra; end
if isfield(par,'betascale'); betascale = par.betascale; end
if isfield(par,'betascale2'); betascale2 = par.betascale2; end
if isfield(par,'subscale'); subscale = par.subscale; end
if isfield(par,'ramin'); ramin = par.ramin; end
if isfield(par,'ramax'); ramax = par.ramax; end
if isfield(par,'scale'); scale = par.scale; end
if isfield(par,'tolrrate'); tolrrate = par.tolrrate; end
if isfield(par,'uptol'); uptol = par.uptol; end
if isfield(par,'maxtime'); maxtime = par.maxtime; end
if isfield(par,'redfrequence'); redfrequence = par.redfrequence; end
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
parPR = BQP_parameter(A,b,bidx,par);
parPR.tolrrate = tolrrate;
if isfield(par,'incmax'); parPR.incmax = par.incmax; end
Rie_proj = @(C,X,parPR)BQP_proj(C,X,parPR);
Rie_retrac = @(X,H,s,tol,parPR)BQP_retrac(X,H,s,tol,parPR);
%% Initialization
n = size(Q,1);
parPR.checkres = 1;
if isfield(par,'x0')
    R0 = 1e-5*randn(n,r);
    % R0 = zeros(n,r);
    R0(:,1) = 2*par.x0-1;
    X.R = R0;
    X = Rie_retrac(X,0,0,1e-5,parPR);
else
    X.R = randn(n,r);
    X = Rie_retrac(X,0,0,1e-5,parPR);
end
parPR.checkres = 0;
%% Parameter for BB
parBB.tolg = 1e2;
parBB.verbose = 0;
parBB.safegard = 1;
parBB.eigratio = 1e2;
parBB.maxiter = 5e4;
parBB.alp = 1;
parBB.adstop = 1;
parBB.checkrank = 0;
parBB.connec = 1;
parBB.plotyes = 0;
parBB.scaler = nC;
parBB.uptol  = uptol;
if isfield(par,'safegard'); parBB.safegard = par.safegard; end
%% Print header
if verbose
    fprintf('\n ite |   fval     |  beta  |  pfeas    dfeas    comp     pdgap |   ra   |gradnorm|Bite|  time  |rank|');
end
%% Initialization
lambda = zeros(n+1,n+1);
pfeas  = inf;
BBiterTotal = 0;
retrac_flag = 0;
l = parPR.l;
%% Main loop
for i = 1:maxiter
    %% pre
    lambda0 = lambda;
    pfeas0  = pfeas;
    f_g = @(X)f_and_g(X,Qs,cs,lambda0,beta,zid);
    %% update primal variable
    [X,lam,iter,~,~,gradnorm] = Rie_BB_Sub(f_g,Rie_proj,Rie_retrac,X,ra,parBB,parPR);
    BBiterTotal = BBiterTotal + iter;
    %% update dual   variable
    R = X.R;
    R(:,1) = R(:,1)+1;
    R = R/2;
    RR = R*R';
    RR1 = [1,R(:,1)';R(:,1),RR];
    Z = RR1-lambda0/beta;
    Zp = Projpoly(Z,zid);
    lambda = beta*(Zp-Z);%-beta*Z.*(Z<=0);
    %% compute residue
    [pfeas,dfeas,comp,pdgap,pfval,dfval,V,d] = BQP_res(X,Zp,lambda*nC,2*lam*nC,Q,c,parPR);
    subfeas  = subscale*max(dfeas,comp);
    res = max([pfeas,dfeas,comp]);
    ttime = etime(clock,tstart);
    if verbose && mod(i,vergap) == 0
        fprintf('\n %3d |%+6.5e|%3.2e|%3.2e %3.2e %3.2e %3.2e|%3.2e|%3.2e|%4d|%3.2e|',i,pfval,beta,pfeas,dfeas,comp,pdgap,ra,gradnorm,iter,ttime);
        fprintf('%3d |', size(X.R,2));
    end
    %% check dual feas
    if (0.9*pfeas < subfeas) && (res > tol)
        ra = max(ramin,ra/2);
        if retrac_flag == 0
            [X.R,~,retrac_flag] = rankinc(X.R,V,d/nC,f_g,Rie_proj,Rie_retrac,parPR);
        end
    end
    %% check stopping criterion
    if res < tol
        fprintf('\n Residue = %3.2e < tol = %3.2e, Stop! \n',res,tol);
        info.dualAR = 2*reshape(lam(l+1:end),[],size(X.R,2))*nC;
        info.retrac_flag = 0;
        break;
    end
    if retrac_flag == 1
        fprintf('\n Retrac singular, Stop! \n');
        info.retrac_flag = 1;
        break
    end
    if ttime > maxtime
        fprintf('\n Reach time limit, Stop! \n');
        break
    end
    %% update beta
    if (iter > 1*betascale    && (pfeas <0.01*subfeas*betascale2 || pfeas < 0.1*tol)) ||...
            (iter > 10*betascale   && (pfeas <0.05*subfeas*betascale2 || pfeas < 0.5*tol)) ||...
            (iter > 100*betascale  && (pfeas <0.1 *subfeas*betascale2 || pfeas < 0.8*tol)) ||...
            (iter > 1000*betascale && (pfeas <0.7 *subfeas*betascale2 || pfeas < 0.9*tol))
        beta = max(betamin,beta/2);
        ra = min(ramax,ra*2);
    elseif pfeas > 0.7*pfeas0
        beta = min(betamax,beta*1.25);
    end
    %% update ra
    if 0.1*pfeas > subfeas
        ra = min(ramax,ra*2);
    end
    %% reduce rank
    if size(X.R,2)>1 && mod(iter,redfrequence) == 0
        [X,~] = rankred(X,Rie_retrac,parPR);
    end
end
%% Record results
fval(1) = pfval;
fval(2) = dfval;
X.RR1   = RR1;
info.pfeas = pfeas;
info.dfeas = dfeas;
info.comp  = comp;
info.pdgap = pdgap;
info.ttime = ttime;
info.ALMite = i;
info.BBite  = BBiterTotal;
info.lambda = lambda*nC;
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
function [f,eugf,pfeas,RR1] = f_and_g(X,Q,c,lambda,beta,zid)
R = X.R;
R(:,1) = R(:,1)+1;
R = R/2;
r = size(R,2);
e1t = zeros(1,r); e1t(1) = 1;
R1 = [e1t;R];
RR1 = R1*R1';
Z = RR1-lambda/beta;
Zp = Projpoly(Z,zid);
Zm = Z-Zp;
QR = Q*R;
f = Produc(QR,R)+2*c'*R(:,1)+beta*norm(Zm,'fro')^2/2;
eugf = QR+beta*Zm(2:end,:)*R1;
eugf(:,1) = eugf(:,1)+c;
pfeas = norm(RR1-Zp,'fro')/(1+norm(RR1,'fro'));
end