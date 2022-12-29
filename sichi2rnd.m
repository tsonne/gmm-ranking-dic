function y = sichi2rnd(a,b,Nrow,Ncol)
% scaled inverse chi squared random numbers
% scaled inverse chi squared (a,b)
% a = degrees of freedom (aka nu)
% b = scale paramter (aka tau^2, s^2)
% ts 2015-02-19
y = (b*a)./chi2rnd(a,Nrow,Ncol);