function mpfn = renewgbt(mpfn,deltagbt)
%% constants
[n,~] = size(mpfn.bus);
[m,~] = size(mpfn.branch);
[~,M] = size(mpfn.smat);
id = mpfn.branch(:,11) == 1;                    % armed branches id
ref = find(mpfn.bus(:,2) == 3);% reference node
nonref = 1:n;
nonref(ref) = [];
%% renew g,b
deltg = zeros(m,1);
deltb = zeros(m,1);
deltg(id==1) = deltagbt(1:sum(id),1);
deltagbt(1:sum(id)) = [];
deltb(id==1) = deltagbt(1:sum(id),1);
deltagbt(1:sum(id)) = [];
%% renew
mpfn.gb = mpfn.gb + [deltg deltb];

%% correct unreasonable g,b
glist = mpfn.gb(:,1); blist = mpfn.gb(:,2);
glist(glist<0) = rand(1)*0.025;
blist(blist>0) = -rand(1)*0.025;

%% calculate z = y^-1
znew = 1./(glist+1j*blist);
mpfn.branch(:,3) = real(znew);
mpfn.branch(:,4) = imag(znew);
for i = 1:M
    mpfn.q(i).branch = mpfn.branch;
end

%% renew angle
for i = 1:M
    thetatp = mpfn.thetamat(nonref,i) + deltagbt(1:n-1);
    mpfn.q(i).bus(nonref,9) = thetatp*180/pi;
    mpfn.thetamat(nonref,i) = thetatp;
    deltagbt(1:n-1) = [];
end
end