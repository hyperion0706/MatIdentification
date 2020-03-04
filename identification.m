%% basic and fine identifications
clear;
fprintf('[Identification] Programme starts .........\n')
%% import data
fprintf('[Identification] Data processing .........\n')
pf = loadcase('case33bw');                % test case
load = readLD;                                  % power load origin file
freq = 1/10;                                    % times per hour
err  = 0.001;                                   % additional error
tol = 10^-8;                                    % tol for ac flow calculation
qlist = data_processing(pf,load,freq,err,tol);  % build datasets


%% Basic identification
fprintf('[Identification] Basic identification .........\n')
pect = 0.03;                                    % gammar: remove branch thro
thro = 10^-8;                                   % decide stop basic identification
[mpfn,flag] = basicidentify(qlist,pect,thro);   % basic identification
                                                % flag == 0: wrong topo
                                                % detected
if flag == 0
    fprintf('[Basic identification] Wrong topology detected. REDUCE PECT!\n We terminate the programme to save time\n')
    quit();
end
tic;
%% select last 30 data
mpfn.q = mpfn.q(:,end-29:end);
mpfn.smat = mpfn.smat(:,end-29:end);
mpfn.dsmat = zeros(size(mpfn.smat));
mpfn.vmat = mpfn.vmat(:,end-29:end);
mpfn.thetamat = mpfn.thetamat(:,end-29:end);

%% Fine identification
fprintf('[Identification] Fine identification .........\n')
%% Fine identification SETTINGS
MAXITER = 25;             % maximum iteration time
varsigma = 0.01;          % decide whether to remove branches
xi = 0.05;                % threshold to remove branches
varphi = 1*10^-10;        % decide whether to end iterations

%% constants
[n,~] = size(mpfn.bus);
[m,~] = size(mpfn.branch);
[~,M] = size(mpfn.smat);

%% initialization
deltapqlist = deltaPQ(mpfn);
dftp = norm(deltapqlist,2);   % throshold 

%% Error evaluation 0
[gberror,wrong_add,wrong_miss] = errorevaluate(mpfn,qlist(1));
gberror = gberror*100;  % display in 100%
fprintf('[fine identification] Starts with error %e\n',dftp);
fprintf('MAPE: g: %f%%, b: %f%%; Wrong branch +/-: %d,%d\n',gberror(1),gberror(2),length(wrong_add),length(wrong_miss));

for T = 1:MAXITER
    %% Pseudo-Power Flow Calculation
    mpfn = pseudopf(mpfn);
    
    %% Jacobian matrix
    % |A B C| 
    % |D E F|
    Jmat = Jacobianmat(mpfn);

    %% Generalized reverse
    invJmat = pinv(Jmat);   
    % Delta S = Delta P,Q
    deltapqlist = deltaPQ(mpfn);
    % delta g,b,theta
    deltagbt = invJmat*deltapqlist;
    % renew g,b and theta
    mpfn = renewgbt(mpfn,deltagbt);
    
    %% calculate throshold
    deltapqlist = deltaPQ(mpfn);
    df = norm(deltapqlist,2);   % throshold
    
    %% remove wrong branches
    if abs(df - dftp)<varsigma
        mpfn = removebran(mpfn,xi);
    end
    
    %% calculate throshold
    deltapqlist = deltaPQ(mpfn);
    df = norm(deltapqlist,2);   % throshold
      
    %% Error evaluation
    [gberror,wrong_add,wrong_miss] = errorevaluate(mpfn,qlist(1));
    gberror = gberror*100;  % display in 100%
    
    %% decide whether to end iteration
    if abs(df - dftp)<varphi && T~=1
        fprintf('[fine identification] Ends with %d iteratinos, with error %e\n',T,df);
        fprintf('MAPE: g: %f%%, b: %f%%; Wrong branch +/-: %d,%d\n',gberror(1),gberror(2),length(wrong_add),length(wrong_miss));
        break;
    else
        fprintf('[fine identification] No.%d ends with error %e\n',T,df);
        fprintf('MAPE: g: %f%%, b: %f%%; Wrong branch +/-: %d,%d\n',gberror(1),gberror(2),length(wrong_add),length(wrong_miss));
    end
    dftp = df;  
end
toc;
fprintf('[Identification] Pragramme successfully ends .........\n')

