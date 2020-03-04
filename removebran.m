function mpfn = removebran(mpfn,val)
%% constants
[m,~] = size(mpfn.branch);
[~,M] = size(mpfn.smat);
glist = mpfn.gb(:,1);
for i = 1:m
    if glist(i)<val
        mpfn.gb(i,1:2) = 0;     % clear .gb
        mpfn.branch(i,[4 11]) = 0;  % clear branch
        mpfn.branch(i,3) = Inf;
    end
end
for i = 1:M
    mpfn.q(i).branch = mpfn.branch;
end
end