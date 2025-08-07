%%*********************************************************
%%This is a matlab implemention of Riemannian Barzilai and Borwein algorithm
%%on fixed rank manifold
%%  min_R {f(X):  X\in M}. M is a manifold
%%  
%% Input: f = objective function
%%        eu_gradf = euclidean gradient of f
%%        X0 = initial point of X
%%        par some parameters
%% Output:  solution X
function [X,lam,iter,ttime,f1,gradnorm] = Rie_BB_Sub(f_g,Rie_proj,Rie_retrac,X0,ra,par,parPR)
tstart = clock;  
if ~exist('par'); par = []; end  
tolg = 1e-3;
uptol = 1e-8;%
maxiter = 100;  
verbose = 1;  
M = 10; %% linesearch parameter  
gamma = 1e-4;                                                     
epsilon = 1e-10; 
sigma = 0.5;
alp = 1; 
checkrank = 10; 
eigratio = 10;    
safegard = 20;
plotyes = 0;
printfre = 1;
adstop = 1;
connec = 0;
scaler = 1;
if isfield(par,'tolg') tolg = par.tolg; end
if isfield(par,'maxiter') maxiter = par.maxiter; end
if isfield(par,'verbose') verbose = par.verbose; end
if isfield(par,'M') M = par.M; end
if isfield(par,'gamma') gamma = par.gamma;end
if isfield(par,'epsilon') epsilon = par.epsilon; end
if isfield(par,'alp') alp = par.alp; end
if isfield(par,'checkrank') checkrank = par.checkrank; end
if isfield(par,'rankratio') eigratio = par.eigratio; end
if isfield(par,'safegard') safegard = par.safegard; end
if isfield(par,'plotyes') plotyes = par.plotyes; end
if isfield(par,'printfre') printfre = par.printfre; end
if isfield(par,'adstop') adstop = par.adstop; end
if isfield(par,'connec') connec = par.connec; end
if isfield(par,'uptol') uptol = par.uptol; end
if isfield(par,'scaler') scaler = par.scaler; end
X = X0;
%%
%% initial values
%%
[f0,eugf,~,X.RR1] = f_g(X);
f1 = f0; 
histf = [f1];  
[grad0,lam] = Rie_proj(eugf,X,parPR); 
grad1 = grad0;
gradnorm = norm(grad1,'fro');
infeas = inf;
newtol = 0;
Nback = 0;
Slowp = 0;
lastwarn('');
iter = 0;
%%
%% main loop
%%
if verbose
    fprintf('\niter |fval     |reldiff  |gradnorm |stepsize  |  time');
end
for i = 1:maxiter
    iter = i;
    %%
    %% when alp is too large or too small
    %%
    if (alp<=epsilon)||(alp>=1/epsilon)
        if gradnorm>1
            alp =1;
        elseif (gradnorm<=1)&&(gradnorm>=1e-5)
            alp = 1/gradnorm;
        else
            alp = 1e5;
        end
    end
    %%
    %% none monotone line search
    %%
    step = 1/alp;
    tolr = max(min(gradnorm/1e2,1/i^3),1e-9);
    Xp = Rie_retrac(X,grad1,-step,tolr,parPR);
    f0 = f1;
    [f1,eugf,pfeas,Xp.RR1] = f_g(Xp);
    j = length(histf);
    ag = 1;
    while f1> max(histf(max(1,j-M):j))-gamma*step*gradnorm^2
        step = step*sigma;
        Xp = Rie_retrac(X,grad1,-step,tolr,parPR);
        [f1,eugf,pfeas,Xp.RR1] = f_g(Xp);
        ag = ag+1;
        if ag>10
            break;
        end
    end
    %%
    %% update iterates
    %%
    X = Xp;
    histf = [histf,f1];
    grad0 = grad1;
    [grad1,lam] = Rie_proj(eugf,X,parPR);
    gradnorm = norm(grad1,'fro');
    infeas0 = infeas;
    infeas = gradnorm;
    reldiff = (f1-f0)/(eps()+abs(f0));
    %%
    %% plot iterates
    %%
    ttime  = etime(clock,tstart);
    if (verbose)
        if mod(i,printfre)==0 % reduce plot times
            fprintf('\n %2d  | %6.10e|%3.2e|%3.2e| %3.2e| %3.2e',i,f1,reldiff,gradnorm,step,ttime);
        end
    end
    [warnMsg, ~] = lastwarn;
    %%
    %% check stopping criterion
    %%
    tol = min([tolg,ra*max(1e-6,pfeas)]);%
    tol = max(tol,1e-10);%!!
    % tol = 1e-7;%!!
    if pfeas <= 1e-6
        tol = max(uptol,tol);
    end
    if (gradnorm <= tol)||(gradnorm <= newtol)||~isempty(warnMsg)||(i == maxiter)
        if (verbose)
            if gradnorm <= tol
                fprintf('\n gradnorm=%3.2e<tolg=%3.2e, gradnorm = %3.2e ,Convergent!',gradnorm,tol,gradnorm);
            elseif (gradnorm <= newtol)
                fprintf('\n gradnorm=%3.2e<newtol=%3.2e, with little progress',gradnorm,newtol);
            elseif i == maxiter
                fprintf('\n maxiter reached!');
            else
                fprintf('\n warning detected stop!');
            end
        end
        break;
    end
    
    %%
    %% increase tolerance while little progress
    %%
    if adstop
        if infeas < 1e3*tol
            if infeas0 > infeas % grad increase
                Nback = Nback+1;
            end
            if infeas < 10*tol % grad close to tol
                Slowp = Slowp+1;
            end
        end
        newtol = floor((Slowp+Nback)/50)*2*tol;%
    end
    %%
    %% riemannian connection
    %%
    %clear S;%
    S = -grad0;
    if connec == 1
        [S,lam] = Rie_proj(S,X,parPR);
    end
    K = S+grad1;
    if mod(i,2)==1
        alp = abs(Produc(S,K))/(step*Produc(S,S));
    else
        alp = Produc(K,K)/(step*abs(Produc(S,K)));
    end
    %%
    %% update rank
    %%
    [~,r] = size(X.R);
    if (mod(i,checkrank)==0)&&(r>safegard)
        [X,reduced] = rankred(X,Rie_retrac,parPR);
        if size(X.R,2) == 1 % rank is reduced to 1, end in advance
            fprintf('\n rank is reduced to 1, stop in advance!');
            X = Rie_retrac(X,0,0,1e-8,parPR);
            break;
        end
        if reduced == 1
            X = Rie_retrac(X,0,0,1e-8,parPR);
            [f1,eugf,~,X.RR1] = f_g(X);
            histf = [f1];
            [grad0,lam] = Rie_proj(eugf,X,parPR);
            grad1 = grad0;
            alp = par.alp;
%             Nback = 0;
%             Slowp = 0;
%             newtol = 0;
        end           
    end
end
% fprintf('\n gradnorm = %3.2e, pfeas = %3.2e',gradnorm,pfeas);
end



