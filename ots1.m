clear all;
close all;
clc;
%Исходные данные для 8 варианта
q = 8;
T = 2;%мкс
f0 = 20;%МГц
Ns = 10;
f = f0 + (0:q-1)/(T);
dt = (1/max(f))/Ns;
t = 0:dt:T;
E = 15;
nfig = 1;
%Сигнал шум соотношение
SNRdB = 5:13;  
SNRdBtheor = 5:13;
SNRtheor = 10.^(SNRdBtheor/10);
 
 
Petheor = zeros(1, length(SNRdBtheor));
for l=1:q+1
   
    Petheor = Petheor + nchoosek(q+1, l) * (-1)^(l+1) * (1/(1+l))*exp(-l/(l+1)*SNRtheor);
    
end
 
%базисные ф-и
 
disp(size(SNRdBtheor));
ss = zeros(q,length(t));
sc = zeros(q,length(t));
ps = zeros(q,length(t));
pc = zeros(q,length(t));
 
 
 
for i = 1:q
    ps(i,:) = sqrt(2/T)*sin(2.*pi.*f(i)*t); 
    pc(i,:) = sqrt(2/T)*cos(2.*pi.*f(i)*t);
 
    ss(i,:) = sqrt(E)*ps(i,:);
    sc(i,:) = sqrt(E)*pc(i,:);
end
 
 
p = [ps;pc];
l = sum(pc(2,:).*pc(2,:))*dt;
R = zeros(2*q, 2*q);
 
for i = 1:2*q
    for k = i:2*q
        R(i,k)= sum(p(i,:).*p(k,:))*dt;
        R(k,i) = R(i,k);
    end
end
figure(nfig);
nfig = nfig + 1;
bar3(R);
 
si = zeros(q,length(t));
for i = 1:q
    si(i,:) = sqrt(2*E/T)*cos(2 * pi * f(i) .* t);
end
 
 
%функции для сигнала и базисных функци
nErrMax = 20;
for nSNR = 1:length(SNRdB)
    nErr = 0;
    nRun = 0;
    SNR = 10^(SNRdB(nSNR)/10);
    sigma =  sqrt(sum(sum((si).^2))/(2*q*SNR));
    
 %цикл моделирования при одном значении отношения синал шум
    while nErr < nErrMax
        
        i = floor(q*rand)+1;
        phas = 2*pi*rand;
        
        %вычисление сигнала на выходе канала
        r = sin(phas)*ss(i,:) + cos(phas)*sc(i,:) + sigma*randn(1,length(t));
 
       for k = 1:q
            rc(k) = sum(r.*pc(k,:))*dt;
            rs(k) = sum(r.*ps(k,:))*dt;
       end
       
       %формирование решения
       
        [resi, ind] = max(rs.^2 + rc.^ 2);
        if (ind)~= i
             nErr = nErr + 1;
             disp([SNRdB(nSNR),nErr,nErrMax,nRun]);
        end
         nRun = nRun + 1;
    end
  
    Pexp(nSNR) = nErr/nRun;
end
    figure(nfig);
    nfig = nfig+1;
    semilogy(SNRdBtheor, Petheor, 'b', SNRdB, Pexp,'ro-');
    grid on;
