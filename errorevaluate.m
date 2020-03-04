function [gberror,wrong_add,wrong_miss] = errorevaluate(mpfn,qcorrect)
%% Error evaluation
correctbran = qcorrect.branch;
armedbran = find(correctbran(:,11)==1);
corry = 1./(correctbran(:,3)+1j*correctbran(:,4));
correctgb = [real(corry) imag(corry)];
% g,b error(MAPE)
pe = (mpfn.gb(armedbran,:)-correctgb(armedbran,:))./correctgb(armedbran,:);
gberror = mean(abs(pe));
% topology
nowarmedbran = find(mpfn.gb(:,1)+1j*mpfn.gb(:,2)~=0);
[~,inow,icorr] = intersect(nowarmedbran,armedbran);
wrong_add = nowarmedbran;
wrong_add(inow) = [];
wrong_miss = armedbran;
wrong_miss(icorr) = [];
end