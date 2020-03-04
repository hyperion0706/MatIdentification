function Jmat = Jacobianmat(mpfn)
%% constants
[n,~] = size(mpfn.bus);
[m,~] = size(mpfn.branch);
[~,M] = size(mpfn.smat);
Ybus = makeYmatrix(mpfn.q(1));
%% Jacobian matrix
    % |A B C| 
    % |D E F|
%% M matrix
miJ = zeros(n,m);
for J = 1:m
    miJ(mpfn.branch(J,1),J) = 1;
    miJ(mpfn.branch(J,2),J) = -1;
end

%% build Jacobian matrix
id = mpfn.branch(:,11) == 1;                    % armed branches id
ref = find(mpfn.bus(:,2) == 3);% reference node
nonref = 1:n;
nonref(ref) = [];
Jmat = zeros(2*n*M,sum(id)+(n-1)*M);
for u = 1:M
    %% data
    vlist = mpfn.vmat(:,u);
    thetalist = mpfn.thetamat(:,u);
    Smat = mpfn.smat(:,u);

    %% Mat A/B
    matA = zeros(n,m);
    matB = zeros(n,m);  
    for J = 1:m
        for i = mpfn.branch(J,1:2)
            to = mpfn.branch(J,1:2);
            to(to == i) = [];
            matA(i,J) = abs(miJ(i,J))*(vlist(i)^2 - ...
                vlist(i)*vlist(to)*cos(thetalist(i) - thetalist(to)));
            matB(i,J) = -abs(miJ(i,J))*vlist(i)*vlist(to)*...
                sin(thetalist(i) - thetalist(to));
        end
    end

    %% Mat E/D
    matE = -matA;
    matD = matB;

    %% mat X: related to branches
    matX1 =[matA; matD];
    matX2 =[matB; matE];

    matX1 = matX1(:,id);
    matX2 = matX2(:,id);

    matX = [matX1 matX2];

    %% Mat C/F
    % H M matrix in power flow calculation
    H0 = zeros(n);
    M0 = zeros(n);

    for i = 1:n
        for j = 1:n
            if i == j
                H0(i,j) = vlist(i)*(imag(Ybus(i,i)) + imag(Smat(i))/(vlist(i)^2))*vlist(j);
                M0(i,j) = vlist(i)*(real(Ybus(i,i)) - real(Smat(i))/(vlist(i)^2))*vlist(j);
            else
                H0(i,j) = vlist(i)*(imag(Ybus(i,j))*cos(thetalist(i)-thetalist(j))-real(Ybus(i,j))*...
                    sin(thetalist(i)-thetalist(j)))*vlist(j);         
                M0(i,j) = -vlist(i)*(-real(Ybus(i,j))*cos(thetalist(i)-thetalist(j))...
                    -imag(Ybus(i,j))*sin(thetalist(i)-thetalist(j)))*vlist(j);
            end
        end
    end
    matC = -H0(:,nonref);
    matF = -M0(:,nonref);
    matPQTM = [matC; matF];

    %% fill in Jacobian matrix
    [m1, n1] = size(matX);
    [m2, n2] = size(matPQTM);
    Jmat(1+(m1)*(u-1):(m1)*u,1:n1) = matX;
    Jmat(1+(m1)*(u-1):(m1)*u,1+n1+n2*(u-1):n1+n2*u) = matPQTM;
end
end