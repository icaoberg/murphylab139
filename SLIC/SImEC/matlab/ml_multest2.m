function [p_value,both] = ml_multest2(f1,f2,Method,TestMethod,dependence,sigvalue,pvalues)

%ML_MULTEST2: multiple testing for multivariates
%   [P_VALUE,BOTH]=ML_MULTEST2(F1,F2,METHOD,TESTMETHOD,DEPEDENCE,SIGVALUE,PVALUES)
%   performs multiple testing for two multivariate samples. It can be considered as
%   a univariate test with some well-designed thresholding schemes. TESTMETHOD is the
%   univariate test method to get pvalue for each feature. If TESTMETHOD is 'ready',
%   there already exist pvalues so that the univariate test has not to be done again.
%   METHOD is the method to decide how to control false rejections. It includes:
%       'Un' - same as univariate tests
%       'BH' - the BH method (may be the best choice).
%       'Bn' - The Benferroni method
%       'Hm' - The Holm method
%   If the variables are dependent, set DEPENDENCE to 1, otherwise to 0. But if most variables
%   are positively dependent, it is better to set DEPENDENCE to 0. SIGVALUE is the
%   significance level to reject hypothesis. PVALUES is the calculated values, which is used
%   if TESTMETHOD is 'ready'.
%
%   Written by T. Zhao

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

feature_size=size(f1,2);

if dependence==1
    Cm=sum(1./[1:feature_size]);
else
    Cm=1;
end

switch(TestMethod)
case 'ready'
    p=pvalues;
    
case 'ttest2'
    for k=1:feature_size
        s1=f1(:,k);
        s2=f2(:,k);
        [h, p(k), ci]=eval(strcat(TestMethod,'(s1,s2)'));
        
    end
    
otherwise
    for k=1:feature_size
        s1=f1(:,k);
        s2=f2(:,k);
        p(k)=eval(strcat(TestMethod,'(s1,s2)'));
        
    end
end

[sorted_p index]=sort(p);
rejections=zeros(1,feature_size);

switch(Method)
case 'Un'
    p_value=sorted_p(1);
    both=[index;sorted_p<=sigvalue;sorted_p];
case 'BH'

    for(i=feature_size:(-1):1)
        if(sorted_p(i)<=(i*sigvalue/feature_size/Cm))
            rejections(1:i)=1;
            break;
        end
    end
    
    p_value=min(1,min(sorted_p*Cm*feature_size./(1:feature_size)));
    both=[index; rejections];
case 'Bn'
    t=sigvalue/feature_size;
    p_value=min(sorted_p(1)*feature_size,1);
    both=[index; sorted_p<=t; sorted_p];
case 'Hm'
    for(i=1:feature_size)
        if(sorted_p(i)>sigvalue/(feature_size-i+1))
            rejections=rejections+1;
            rejections(i:feature_size)=0;
            break;
        end
    end
    
    p_value=min(sorted_p(1)*feature_size,1);
    both=[index; rejections; sorted_p];
end