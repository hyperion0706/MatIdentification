function Ymat = makeYmatrix(varargin)
mpc = varargin{1};
% makeYmatrix
% get mpc
%% get mpc information
%mpc = case91;
bus = mpc.bus;
branch = mpc.branch;
N = length(mpc.bus(:,1));       % n of bus
M = length(mpc.branch(:,1));    % n of branch
baseMVA = mpc.baseMVA;

%% Taps and shifts
tap = ones(M,1);
i = find(branch(:,9));
tap(i) = branch(i,9);                           % tap
tap = tap .* exp(1j*pi/180 * branch(:, 10));    % shift
tap = conj(tap);                                % conj 

%% incidence matrix
A_b = sparse(branch(:,1:2),[1:M,1:M],[1./tap,-ones(M,1)],N+1,M); 
% i,j,v generate sparse matrix, series
A_s1 = sparse([branch(:,1) (N+1)*ones(M,1)],[1:M,1:M],[1./tap,-ones(M,1)],N+1,M);%in
A_s2 = sparse([branch(:,2) (N+1)*ones(M,1)],[1:M,1:M],[ones(M,1),-ones(M,1)],N+1,M);%out
% sh
A_sh = sparse([[1:N]',(N+1)*ones(N,1)],[1:N,1:N],[ones(N,1),-ones(N,1)],N+1,N);
% shunt
A = [A_b A_s1 A_s2 A_sh];

%% calculate series elements
stat = branch(:, 11);
yij = 1./(branch(:,3)+1j * branch(:,4)); % 1j: imaginary, yij

%% calculate shunt elements
y_lc = 1j*branch(:,5)*0.5; %% line charging susceptance
y_sh = (bus(:,5) + 1j*bus(:,6))/baseMVA;

y_list = [yij; y_lc; y_lc; y_sh];
y = spdiags(y_list,0,length(y_list),length(y_list));

%% calculate shunt consider ground
Y = (A)*y*A';
Y0 = Y(1:N,1:N);

if nargin == 1
    Ymat = Y0;
elseif nargin == 2
    if varargin{2} == 1
        Ymat = Y;
    elseif varargin{2} == 2
        Ymat = Y0;
    elseif varargin{2} == 3
        dat = find(abs(Y0-makeYbus(case91))>10^(-8));
        if isempty(dat)
        fprintf('Same as makeYbus() in MATPOWER\n')
        else
        fprintf('Disagree as makeYbus() in MATPOWER, item %d',dat)  
        end
    elseif varargin{2} == 4
        plotmatrix(real(Y))
    elseif varargin{2} == 5
        plotmatrix(imag(Y))
    end
end
