clear all;
close all;
clc;
q = 8;
T = 2;
f0 = 20;
f = f0 + ((0:q-1)/T);
dt = (1/max(f))/10;
t = 0:dt:T;
E = 10;
eps=0;
m = sqrt(eps/2);
D = (1-eps)/2;
d=sqrt(D);
si = zeros(q,length(t));
for i = 1:q
    si(i,:) = sqrt(2*E/T)*cos(2*pi*f(i)*t);
end
s_cos = zeros(q,length(t));
s_sin = zeros(q,length(t));
f_sin = zeros(q,length(t));
f_cos = zeros(q,length(t));
 
for i = 1:q
    f_cos(i,:) = sqrt(2/T)*cos(2*pi*f(i)*t);
    f_sin(i,:) = sqrt(2/T)*sin(2*pi*f(i)*t); 
    s_cos(i,:) = sqrt(2*E/T)*f_cos(i,:);
    s_sin(i,:) = sqrt(2*E/T)*f_sin(i,:);
end
SNRdB = 1:14;
SNRdBtheor = 1:14;
SNRtheor = 10.^(SNRdBtheor/10);
nErrMax = 50;
for nSNR = 1:length(SNRdB)
    nErr = 0;
    ntest = 0;
    SNR = 10^(SNRdB(nSNR)/10);
    sigma2=(sum(sum(si.^2)))/(q*2*SNR);
    sigma=sqrt(sigma2);
    while nErr < nErrMax
        x = normrnd(m,d);
        y = normrnd(m,d);
        Mu = sqrt((x^2)+(y^2));
        XX = q*rand + 1;
        i = floor(XX);
        tetta = rand();
        r = Mu*cos(2*pi*tetta)*s_cos(i,:) + Mu*sin(2*pi*tetta)*s_sin(i,:) + sigma*randn(1,length(t));
        for a = 1:(q)
            r_c(a) = sum(r.*f_cos(a,:))*dt;
            r_s(a) = sum(r.*f_sin(a,:))*dt;
        end
        [index, i_result] = max((r_s.^2) + (r_c.^ 2));
        if (i_result)~= i
             nErr = nErr + 1;
             clc
             disp([SNRdB(nSNR), ntest, nErr, nErrMax]);
        end 
         ntest = ntest + 1;
    end 
    Pe_experement(nSNR) = nErr/ntest;
end 
Pe_theory = zeros(1, length(SNRdBtheor));
for m=1:q-1
   koef = nchoosek(q-1,m);
   mult = koef*((-1)^(m+1));
   mult2=mult*(1./(1+m+m*(1-eps)*SNRtheor));
   mult3 = mult2.*exp(((-1*m*eps*SNRtheor)./(1+m+m*(1-eps)*SNRtheor)));
   Pe_theory = Pe_theory + mult3;
end
    figure(1);
    axis('square');
    semilogy(SNRdBtheor, Pe_theory, 'black', SNRdB, Pe_experement,'red.-','MarkerSize',10);
    xlabel('SNRdB')
    ylabel('Pe')
    grid on;
