
clear all %#ok<CLALL>
clc
N=30; % Number of search agents

Function_name='F1'; 
T=100;
runs =1;
[lb,ub,dim,fobj]=Get_F(Function_name);
for iter = 1:runs
    [Loc_cssma,Bfit_cssma,CNVG_cssma]=EISMA(N,T,lb,ub,dim,fobj) ;   
end



