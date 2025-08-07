%% preparation
close all;
clear
lastwarn('');
warning off;
rng('default');
setup_path;

%% solver
useALM = 1;
useSDPNAL = 0;

%% pars
tol = 1e-6;
record = 0;
probtype = 'BIQ';
fname = 'bqp1000.1';

%% record results
rrALM   = {"data","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","time"};
rrNAL   = {"data","n","m","binary","rank","pobj","dobj","pfeas","dfeas","comp","pdgap","iter","itersub","iterA","time"};

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
l = length(bidx);
nc = m+m*n+l+1;
fprintf('\n Start testing %s problems: n=%2d, m=%2d, binary=%2d\n',probtype,n,m,l);

%% RNNAL
if useALM
    par.tol = tol;
    par.beta = 1e0;
    par.ra  = 1e2;
    r = max(5,round(1.5*sqrt(2*nc)+2));
    fprintf("\n ----- START RADNN: initial rank = %4i",r);
    [obj_ALM,X_ALM,info_ALM] = RiNNAL(Q,c,A,b,bidx,r,zid,par);
    if record
        nlALM = {fname,n,m,length(bidx),size(X_ALM.R,2),obj_ALM(1),obj_ALM(2),...
            info_ALM.pfeas,info_ALM.dfeas,info_ALM.comp,info_ALM.pdgap,info_ALM.ALMite,info_ALM.BBite,info_ALM.ttime};
        rrALM = [rrALM;nlALM];
        cname = strcat(pwd,"/results/BIQ/",probtype,"-ALM-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
        save(cname,"rrALM");
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
        nlNAL = {fname,n,m,length(bidx),rank(X_NAL{1},1e-6),obj_NAL(1),obj_NAL(2),...
            max([info_NAL.etaRp,info_NAL.etaK1,info_NAL.etaK2]),...
            info_NAL.etaRd,max([info_NAL.etaC1,info_NAL.etaC1]),...
            info_NAL.relgap,info_NAL.iterSSN,info_NAL.iterSSNsub,info_NAL.iterADM,ttime_NAL};
        rrNAL = [rrNAL;nlNAL];
        cname = strcat(pwd,"/results/BIQ/",probtype,"-NAL-",string(datetime('now','Format','yyyy-MM-dd--HH-mm-ss')),".mat");
        save(cname,"rrNAL");
    end
end