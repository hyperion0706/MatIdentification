% make S matrix on buses
function [Smat,SLmat,Gmat] = makeSm(mpc)
%% import 
%mpc = case30;
bus = mpc.bus;
gen = mpc.gen;
%% 
Smat = -bus(:,3)-1j*bus(:,4);
Smat(gen(:,1)) = Smat(gen(:,1))+gen(:,2)+ 1j*gen(:,3);
Smat = Smat/mpc.baseMVA;
SLmat = (-bus(:,3)-1j*bus(:,4))/mpc.baseMVA;
Gmat = gen(:,2)/mpc.baseMVA + 1j*gen(:,3)/mpc.baseMVA;
