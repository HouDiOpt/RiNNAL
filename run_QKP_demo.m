%% preparation
close all;
clear
lastwarn('');
warning off;
rng('default');
setup_path;

%% solver
useALM = 1;
useALMdiag = 0;
useSDPNAL = 0;

%% pars
tol = 1e-6;
record = 0;
probtype = 'QKP';
fname = 'QKP-n500-p01-beta09';

%% record results
rrALM      = {"data","bscale","Qdensity","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time"};
rrALM_diag = {"data","bscale","Qdensity","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time"};
rrNAL      = {"data","bscale","Qdensity","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time"};

%% data
load(strcat(fname,'.mat'),'data');
Q = data.Q;
c = data.c;
A = data.A;
b = data.b;
bidx = data.bidx;
zid  = [];
m = data.m;
n = data.n;
bs = data.bscale;
Qd = data.Qdensity;
l = length(bidx);
nc = m+m*n+l+1;
fprintf('\n Start testing %s problems: n=%2d, m=%2d, binary=%2d\n',probtype,n,m,l);

%% ALM
if useALM
    par.tol = tol;
    par.projtype = 3;  % 0.SchurChol 3.PCG
    par.ra  = 1e2;     % bscale small 1e1  other 1e2
    par.ramin  = 1e-2; % bscale small 1e-4 other 1e-2
    par.tolrrate = 1e-2;
    r = max(5,round(sqrt(2*nc)+2));
    fprintf("\n ----- START RADNN: initial rank = %4i",r);
    [obj_ALM,X_ALM,info_ALM] = RiNNAL(Q,c,A,b,bidx,r,zid,par);
    if record
        nlALM = {fname,bs,Qd,n,m,length(bidx),size(X_ALM.R,2),obj_ALM(1),obj_ALM(2),...
            info_ALM.pfeas,info_ALM.dfeas,info_ALM.comp,info_ALM.pdgap,info_ALM.ALMite,info_ALM.BBite,info_ALM.ttime};
        rrALM = [rrALM;nlALM];
        cname = strcat(pwd,"/results/QKP/",probtype,"-ALM-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
        save(cname,"rrALM");
    end
end

%% ALM (diag)
if useALMdiag && n<2000
    opts.tol = tol;
    opts.beta = 1e-3;%1
    opts.betamax = 1e5;
    opts.ra = 1e2;
    opts.ramin = 1e0;
    opts.timelimit = 60*60;
    r = max(5,round(sqrt(2*nc)+2));
    fprintf("\n ----- START ALM-diag: initial rank = %4i",r);
    [obj_ALM_diag,X_ALM_diag,info_ALM_diag] = RiNNAL_diag(Q,c,A,b,bidx,r,opts);
    if record
        nlALM_diag = {fname,bs,Qd,n,m,length(bidx),size(X_ALM_diag.R,2),...
            obj_ALM_diag(1),obj_ALM_diag(2),info_ALM_diag.pfeas,info_ALM_diag.dfeas,...
            info_ALM_diag.comp,info_ALM_diag.pdgap,info_ALM_diag.ALMite,info_ALM_diag.BBite,info_ALM_diag.ttime};
        rrALM_diag = [rrALM_diag;nlALM_diag];
        cname = strcat(pwd,"/results/QKP/",probtype,"-ALM-diag-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
        save(cname,"rrALM_diag");
    end
end

%% SDPNAL
if useSDPNAL 
    % convert data
    [blk,At,C,bb] = MBQP_to_SDPNAL(Q,c,A,b,bidx,zid);
    [blk,At,C,bb] = SDPR_to_SDPNAL(blk,At,C,bb);
    % solver
    LL = 0;
    OPTIONS.tol = tol;
    OPTIONS.maxtime = 3600;
    if nc < 1000
        OPTIONS.AATsolve.method = 'direct';
    else
        OPTIONS.AATsolve.method = 'iterative';
    end
    OPTIONS.stopoption = 0;
    tic;
    [obj_NAL,X_NAL,~,~,~,~,~,~,info_NAL,~] = ...
        sdpnalplus(blk,At,C,bb,LL,[],[],[],[],OPTIONS);
    ttime_NAL = toc;
    fprintf('\n SDPNAL+ rank  = %3.0d\n',rank(X_NAL{1},1e-6));
    if record
        nlNAL = {fname,bs,Qd,n,m,length(bidx),rank(X_NAL{1},1e-6),obj_NAL(1),obj_NAL(2),...
            max([info_NAL.etaRp,info_NAL.etaK1,info_NAL.etaK2]),...
            info_NAL.etaRd,max([info_NAL.etaC1,info_NAL.etaC1]),...
            info_NAL.relgap,info_NAL.iterSSN,info_NAL.iterSSNsub,info_NAL.iterADM,ttime_NAL};
        rrNAL = [rrNAL;nlNAL];
        cname = strcat(pwd,"/results/QKP/",probtype,"-NAL-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
        save(cname,"rrNAL");
    end
end
