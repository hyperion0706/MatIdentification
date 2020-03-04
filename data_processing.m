function q = data_processing(pf,load,freq,err,tol)
%% import data
%pf = loadcase('case33bw');
% clear Gs on branches
pf.branch(:,5) = 0;
% run matpower ac power flow
mpopt = mpoption('verbose',0,'out.all',0);
mpf = runpf(pf,mpopt);

%% build data models
% build models, mpf, load data,period, renewable data, renewable point,
% error
%load = readLD;
pqload = pqloads(mpf,load,freq,err,tol);
% mpf list
q = pqload.qlist;
end