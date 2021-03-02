# MatIdentification
Electrical Engineering: A tool to identify topology and parameters of a distribution network in Matlab.

OPEN identification.mlx in MATLAB and it introduces a way to make identification on distribution networks.

Ensure MATPOWER 6.0 or later has been installed on your MATLAB.

The paper is present at
J. Zhang, Y. Wang, Y. Weng and N. Zhang, "Topology Identification and Line Parameter Estimation for Non-PMU Distribution Network: A Numerical Method," in IEEE Transactions on Smart Grid, vol. 11, no. 5, pp. 4440-4453, Sept. 2020, doi: 10.1109/TSG.2020.2979368. https://ieeexplore.ieee.org/document/9027950 

The current code shows a simplified and kernel version of the proposed model. case1github.mat show a example data and results, with simular setups in the paper. 

P, Q additional error err  = 0.01; 

V error, added by reducing pf.tol on MATPOWER, tol = 10^-5; -> 0.01%-level error

24*5 datasets, freq = 1/5, but select last 30 datasets for step 2.

pect = 0.01;                                    % gammar: remove branch thro

thro = 10^-8;                                   % decide stop basic identification

%% Fine identification SETTINGS

MAXITER = 25;             % maximum iteration time

varsigma = 0.01;          % decide whether to remove branches

xi = 0.05;                % threshold to remove branches

varphi = 1*10^-10;        % decide whether to end iterations
