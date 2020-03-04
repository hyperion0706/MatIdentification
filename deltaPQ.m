function deltapqlist = deltaPQ(mpfn)
%% constants
[n,~] = size(mpfn.bus);
[~,M] = size(mpfn.smat);
Ybus = makeYmatrix(mpfn.q(1));
deltapqlist = zeros(2*n*M,1);
for k = 1:M
    Smat = mpfn.smat(:,k);
    vlist = mpfn.vmat(:,k);
    thetalist = mpfn.thetamat(:,k);
    p = zeros(n,1);
    q = zeros(n,1);
    for i = 1:n
        p(i) = real(Smat(i));
        q(i) = imag(Smat(i));
        for j = 1:n
            p(i) = p(i) - vlist(i)*vlist(j)*(real(Ybus(i,j))...
                *cos(thetalist(i)-thetalist(j))+imag(Ybus(i,j))*sin(thetalist(i)-thetalist(j)));
            q(i) = q(i) - vlist(i)*vlist(j)*(real(Ybus(i,j))...
                *sin(thetalist(i)-thetalist(j))-imag(Ybus(i,j))*cos(thetalist(i)-thetalist(j)));
        end
    end
    mpfn.dsmat(:,k) = p+1j*q;       % renew to delta s mat
    deltapqlist(1+(k-1)*2*n:2*k*n,1) = [p;q];
end
end