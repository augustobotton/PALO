T0 = readtable('dados.csv');
T1 = table2array(T0);
h=T1(1,:);
v=T1(2,:);
a=width(h);
T=zeros(a);


for i=1:(a)
    [T(i),Tm,p,rho,ainf,M,mu,Pr,Kn,d,Re,R]=atm_padrao(h(i),v(i),1.5,10);
end

plot(h,v)