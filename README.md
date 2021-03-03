# MatIdentification
## introduction
Electrical Engineering: A tool to identify topology and parameters of a distribution network in Matlab.

OPEN identification.mlx in MATLAB and it introduces a way to make identification on distribution networks. (voltage tol = 10e-8 here)

Ensure MATPOWER 6.0 or later has been installed on your MATLAB.

## paper
The paper is present at
J. Zhang, Y. Wang, Y. Weng and N. Zhang, "Topology Identification and Line Parameter Estimation for Non-PMU Distribution Network: A Numerical Method," in IEEE Transactions on Smart Grid, vol. 11, no. 5, pp. 4440-4453, Sept. 2020, doi: 10.1109/TSG.2020.2979368. https://ieeexplore.ieee.org/document/9027950 

This code is refactored after the paper accepted, and the performance is a bit different from that in the paper. The current version is more fast and clear to read, so I put it here.

## Example
The current code shows a simplified and kernel version of the proposed model. **case1github.mat** and **case1github2.mat** show example data and results, with simular setups in the paper. 
### setups
```
P, Q additional error err  = 0.01; 

V error, added by reducing pf.tol on MATPOWER, tol = 10^-5; -> 0.01%-level error

24*5 datasets, freq = 1/5;
 case1github.mat select last 30 datasets for step 2.
 case1github2.mat select all 120 dataset.
```
Other parameters:
```
pect = 0.01;                                    % gammar: remove branch thro
thro = 10^-8;                                   % decide stop basic identification

%% Fine identification SETTINGS
MAXITER = 25;             % maximum iteration time
varsigma = 0.01;          % decide whether to remove branches
xi = 0.05;                % threshold to remove branches
varphi = 1*10^-10;        % decide whether to end iterations
```

**case1github2.mat outputs. simular setups with case 1 in the paper.**
### Outputs:
```
>> identification
[Identification] Programme starts .........
[Identification] Data processing .........
[data processing] Build No. 50 dataset
[data processing] Build No. 100 dataset
[data processing] Build No. 120 dataset
[Identification] Basic identification .........
[Basic Identification] No.1 ends with error 302.440496
 Wrong branch +/-: 312,0
[Basic Identification] No.2 ends with error 150.451166
 Wrong branch +/-: 184,0
[Basic Identification] No.3 ends with error 50.923270
 Wrong branch +/-: 99,0
[Basic Identification] No.4 ends with error 30.552014
 Wrong branch +/-: 51,0
[Basic Identification] No.5 ends with error 16.561910
 Wrong branch +/-: 28,0
[Basic Identification] No.6 ends with error 10.518838
 Wrong branch +/-: 18,0
[Basic Identification] No.7 ends with error 7.397120
 Wrong branch +/-: 16,0
[Basic Identification] No.8 ends with error 2.289884
 Wrong branch +/-: 15,0
[Basic Identification] No.9 ends with error 2.049158
 Wrong branch +/-: 14,0
[Basic Identification] No.10 ends with error 0.374646
 Wrong branch +/-: 14,0
[Basic Identification] Converged after 11 iterations with error 0.000000e+00
 Wrong branch +/-: 14,0
时间已过 0.089314 秒。
[Identification] Fine identification .........
[fine identification] Starts with error 2.828044e-01
MAPE: g: 27.028450%, b: 31.582200%; Wrong branch +/-: 14,0
[fine identification] No.1 ends with error 4.114671e-01
MAPE: g: 15.028266%, b: 10.842871%; Wrong branch +/-: 14,0
[fine identification] No.2 ends with error 1.682321e-01
MAPE: g: 2.338936%, b: 3.589087%; Wrong branch +/-: 14,0
[fine identification] No.3 ends with error 2.693385e-02
MAPE: g: 0.996655%, b: 0.669106%; Wrong branch +/-: 14,0
[fine identification] No.4 ends with error 4.163705e-02
MAPE: g: 0.854492%, b: 0.738039%; Wrong branch +/-: 14,0
[fine identification] No.5 ends with error 1.015405e-02
MAPE: g: 0.826144%, b: 0.861779%; Wrong branch +/-: 0,0
[fine identification] No.6 ends with error 1.055253e-03
MAPE: g: 0.600930%, b: 0.568718%; Wrong branch +/-: 0,0
[fine identification] No.7 ends with error 1.042050e-03
MAPE: g: 0.606518%, b: 0.573085%; Wrong branch +/-: 0,0
[fine identification] No.8 ends with error 1.042050e-03
MAPE: g: 0.606678%, b: 0.573369%; Wrong branch +/-: 0,0
[fine identification] Ends with 9 iterations, with error 1.042050e-03
MAPE: g: 0.606687%, b: 0.573381%; Wrong branch +/-: 0,0
时间已过 416.990129 秒。
[Identification] Programme successfully ends .........
```
## Exist issues
1. There might be larger estimation errors, when pf.tol>10e-4. This is result from the collinear problem of voltage measurements. In simple words, the voltage records are too simular among datasets. Not algorithms' faults. The method works well when pf.tol <= 10-3, there are some cases start with **case1githubXX.mat**. Yet the algorithm fails to coverage when pf.tol>=10-2, due to basic identification fails to give rough estimations.
2. Experiences on tunning super-parameters are required, especially in case 1 and 3 in the paper. It may influence the accuray and correctness of parameter estimation and topology identification.
3. The performance of the algorithm is related to the quality of datasets. Low colinear and high accurate data is welcomed.
4. We have improved the algorithm into linear form, and the paper is under review. No speed problem later.

## About the usage of the code
You can use the code in your paper/ degree thesis as control groups, comparisons, or part of your methods. Please cite our paper.
Yet, the licence of this code does not cover the paper on IEEE trans. Smart Grid. Do not copy words, figures, or tables to your paper.
The ideas and code of this paper should not be a major part of your innovative article (contribution greater than 20%).  We do not welcome or allow the plagiarism restricted by IEEE regulations.
