function [mpfn,flag] = basicidentify(q,pect,thro)
%% SETTINGS
% throshold
%pect = 0.03;
%thro = 10^-8;
% maximum iteration
MAX_ITER_BASIC = 100;
%% import data
M0 = length(q);
[n,~] = size(q(1).bus);
[m,~] = size(q(1).branch);
% correct topology and parameter
correct_topo = full(makeYbus(q(1))~=0);
tic;

%% build pv pq
PV = zeros(M0,n);
QV = zeros(M0,n);
V  = zeros(M0,n);
pdivv = @(s,x)(real(s)./x.bus(:,8))';
qdivv = @(s,x)(imag(s)./x.bus(:,8))';
vl = @(x)(x.bus(:,8))';
nodeSmatrix = zeros(n,M0);  % matrix records 
for i = 1:M0
    s = makeSm(q(i));
    nodeSmatrix(:,i) = s;
    PV(i,:) = pdivv(s,q(i));
    QV(i,:) = qdivv(s,q(i));
    V(i,:)  = vl(q(i));
end

%% build Gij# and Bij#
% Gijsp: Gij sharp
Gijsp = (V'*V)\V'*PV;
% Bijsp: Bij sharp
Bijsp = -(V'*V)\V'*QV;

%% reduce noise, based on Gij#
% Gij# Bij# temporary
% s.t. i==j: positive, i~=j: negative
Gijtp = Gijsp;
Bijtp = -Bijsp;

% record last Gij#, used for judgement of convergence
Gijlast = Gijtp;
dev_topo = [];
for i = 1:MAX_ITER_BASIC                   
    %fprintf('[Basic Identification] No.%d starts\n',i);
    Gijlast = Gijtp;
    % regression part   
    for j = 1:n
        bit = find(Gijtp(j,:)~=0);          % find non-zero items in i row
        greg = regress(PV(:,j),V(:,bit));   % regress PV,V => G
        breg = regress(QV(:,j),V(:,bit));   % regress QV,V => -B
        % renew the regression to Gijtp/Bijtp
        Gijtp(j,bit) = greg;
        Bijtp(j,bit) = breg;
    end
    
    % de-noise part
    for k = 1:n      
        % backup diag items
        gdiag = Gijtp(k,k);
        bdiag = Bijtp(k,k);
        
        % rows
        % Gij#: find item less than pect
        grec = -Gijtp(k,:)/abs(gdiag);
        ricc = find(grec<pect);  
        % columns
        % Gij#: find item less than pect
        grec = -Gijtp(:,k)/abs(gdiag);
        cicc = find(grec<pect);
        % both less than pect
        icc = intersect(ricc,cicc);
        
        % remove these items for both Gij# and Bij#
        Gijtp(k,icc) = 0;                   % set as zero              
        Bijtp(k,icc) = 0;
        % remove these items for both Gij# and Bij#
        Gijtp(icc,k) = 0;
        Bijtp(icc,k) = 0;
        
        % recover i == j items
        Gijtp(k,k) = gdiag;       
        Bijtp(k,k) = bdiag;
    end
    % symmetrization
    Gijsym = Gijtp;
    Bijsym = Bijtp;
    for ii = 1:n
        for jj = 1:n
            Gijsym(ii,jj) = Gijtp(jj,ii);
            Bijsym(ii,jj) = Bijtp(jj,ii);
        end
    end
    Gijtp = (Gijsym + Gijtp)/2;
    Bijtp = (Bijsym + Bijtp)/2;
    
    % check throshold
    score = sum(sum(abs(Gijlast-Gijtp)));           % deviation
    current_topo = Gijtp~=0;                        % topology now
    dev_topo = current_topo-correct_topo;           % contradictory matrix
    wrong_add = length(find(dev_topo == 1))/2;      % rebundant branches
    wrong_loss = length(find(dev_topo == -1))/2;    % missing branches
    if score<thro                                   % converged
        fprintf('[Basic Identification] Converged after %d iterations with error %e\n Wrong branch +/-: %d,%d\n',i,score,wrong_add,wrong_loss);
        break
    else
        fprintf('[Basic Identification] No.%d ends with error %f\n Wrong branch +/-: %d,%d\n',i,score,wrong_add,wrong_loss);
    end
end
toc;

%% build new matpower model for new branches
% rebundant branches
[rp,cp] = find(triu(dev_topo) == 1);
pitem = find(triu(dev_topo) == 1);
% missing branches
[rn,cn] = find(triu(dev_topo) == -1);
nitem = find(triu(dev_topo) == -1);

% wrong identification indication
flag = 1;
if ~isempty(nitem)
     fprintf('[Basic Identification] Reduce pect\n');
     flag = 0;
end

% rebundant branches in list
fromtolist = q(1).branch(:,1)+1j*q(1).branch(:,2);
tofromlist = q(1).branch(:,2)+1j*q(1).branch(:,1);
addbran = rp + 1j*cp;
[~,iaddbran1,ifromtolist] = intersect(addbran,fromtolist);   % find already in list
[~,iaddbran2,itofromlist] = intersect(addbran,tofromlist);
addbranalreadylist =[iaddbran1; iaddbran2];
fromtolist = [ifromtolist;itofromlist];
rp(addbranalreadylist) = [];
cp(addbranalreadylist) = [];    % remove item already in list

% branch model
mpfbranch = q(1).branch;

% fake branch: newbranitem
egbranch = q(1).branch(1,:);
egbranch(end-3:end) = 0;
fakebranchlist = zeros(length(rp),length(egbranch));
for i = 1:length(rp)
    fakebranchlist(i,:) = egbranch;
end
fakebranchlist(:,1:2) = [rp cp];

% new branch model
mpfbranch = [mpfbranch; fakebranchlist];
[mnew,~] = size(mpfbranch);
yfake = zeros(mnew,1);
for i = 1:mnew
    % obtain y on branches from Gij and Bij
    yfake(i) = -Gijtp(mpfbranch(i,1),mpfbranch(i,2))+1j*Bijtp(mpfbranch(i,1),mpfbranch(i,2));
    if yfake(i) == 0
        mpfbranch(i,11) = 0;
    else
        mpfbranch(i,11) = 1;
    end
end
zfake = 1./yfake;
mpfbranch(:,3) = real(zfake);
mpfbranch(:,4) = imag(zfake);

% renew to all dataset
for i = 1:M0
    q(i).branch = mpfbranch;
end

%% output
% nodeSmatrix: P+iQ
%nodeSmatrix;
% node V
nodeVmatrix = V';
% branch Mat
%mpfbranch;
%% mpfn 
% PQVtheta
mpfn.smat = nodeSmatrix;
mpfn.vmat = nodeVmatrix;
mpfn.thetamat = zeros(n,M0);
% y
mpfn.gb = [real(yfake) imag(yfake)];
% branch
mpfn.branch = mpfbranch;
% bus
mpfn.bus = q(1).bus;
% q
mpfn.q = q;
end