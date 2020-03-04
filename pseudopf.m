function mpfn = pseudopf(mpfn)
[~,M] = size(mpfn.smat);
for i = 1:M
        mpopt = mpoption('verbose',0,'out.all',0);
        mpf = runpf(mpfn.q(i),mpopt);       % run AC power flow calcation 
        mpf.bus(:,8) = mpfn.vmat(:,i);
        mpfn.q(i) = mpf;                    % renew q(i)
        mpfn.thetamat(:,i) = mpf.bus(:,9)*pi/180;  % renew theta for each 
end
end