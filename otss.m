clear all;
close all;
clc;
%Исходные данные для 8 варианта
q = 8;
T = 2;%мкс
f0 = 20;%МГц
f = f0 + ((0:q-1)/T);%набор центральных частот
dt = (1/max(f))/10;
t = 0:dt:T;
E = 10;%энергия принятого сигнала

%формула сигналов si(t),если 0<t<T
si = zeros(q,length(t));
for i = 1:q
    si(i,:) = sqrt(2*E/T)*cos(2 * pi * f(i).* t);
end
%Строим базисные функции и sin и cos составляющие сигналов по формуле
s_cos = zeros(q,length(t));
s_sin = zeros(q,length(t));
f_sin = zeros(q,length(t));
f_cos = zeros(q,length(t));
 
for i = 1:q
    f_cos(i,:) = sqrt(2/T)*cos(2.*pi.*f(i)*t);
    f_sin(i,:) = sqrt(2/T)*sin(2.*pi.*f(i)*t);
    s_cos(i,:) = sqrt(E)*f_cos(i,:);
    s_sin(i,:) = sqrt(E)*f_sin(i,:);
end

%Вычисляем отношение сигнал/шум
SNRdB = 1:15; %набор значений отношения сигнал/шум
SNRdBtheor = 1:13;
SNRtheor = 10.^(SNRdBtheor/10);
 
nErrMax = 20;%назначаем максимальное число ошибок nErrMax в пределах 20..50..100

%Цикл по значениям отношения сигнал/шум в дБ
for nSNR = 1:length(SNRdB)
    nErr = 0;%начальное значение счетчика ошибок
    ntest = 0;%начальное значение счетчика испытаний
    SNR = 10^(SNRdB(nSNR)/10);%Y
    sigma2=(sum(sum(si.^2)))/(q*2*SNR);
    sigma=sqrt(sigma2);
    
 %цикл моделирования при одном значении отношения синал шум
    while nErr < nErrMax
        %моделирование передатчика и канала
        %i = floor(q*rand)+1;%случайно выбираем i на интервале 0,1...,q-1
        %здесь функция rand отвечает за равномерно распределенные случайные числа
        X = q*rand + 1;
        i = floor(X);
        
        %Значение случайной фазы тоже будем задавать через рандом
        tetta = rand;
        
        %Вычисляем сигнал на выходе канала
        r = cos(2*pi*tetta)*s_cos(i,:) + sin(2*pi*tetta)*s_sin(i,:) + sigma*randn(1,length(t));
        
        %Моделирование приемника
        %Для i=0,1,..,q-1 вычисляем значение r_c(i) и r_s(i)
        for a = 1:(q)
            r_c(a) = sum(r.*f_cos(a,:))*dt;
            r_s(a) = sum(r.*f_sin(a,:))*dt;
        end
       
       %Формируем решения
        [index, i_result] = max((r_s.^2) + (r_c.^ 2));
        
        %Фиксация результата
        if (i_result)~= i
             nErr = nErr + 1;%Увеличиваем счетчик числа испытаний
             clc
disp([SNRdB(nSNR), ntest, nErr, nErrMax]);%Отображаем массив на экране
% disp([SNRdB(nSNR),nErr,nErrMax,ntest]);

        end 
         ntest = ntest + 1;
    end %Конец цикла при одном значении сигнал/шум
    
    %Вычисление эксперементальной вероятности ошибки
    Pe_experement(nSNR) = nErr/ntest;%деления числа произошедших ошибок на число испытаний
end %Конец цикла по значениям отношения сигнал/шум в дБ

%Вычисляем теоритическое значение вероятности ошибки
Pe_theory = zeros(1, length(SNRdBtheor));
for l=1:q+1
   % вычисления делаем по действиям
   %1 вычисляем биномиальный коэффициент С
   koef = nchoosek(q+1, l);
   %2 Умножаем коэффициент на (-1) в степени (l+1)
   mult = koef*(-1)^(l+1);
   %3 умножаем получившееся значение на отношение 1/1+l
   mult2=mult*(1/(1+l));
   %4 Умножаем на экспоненту с аргументом ((-1/1+l)*(E/N0))
   mult3 = mult2*exp(-l/(l+1)*SNRtheor);
   % Находим сумму на промежутке от 0 до q-1
   Pe_theory = Pe_theory + mult3;
end
    figure(1);
    %nfig = nfig+1;
    axis('square');
    semilogy(SNRdBtheor, Pe_theory, 'black', SNRdB, Pe_experement,'red.-','MarkerSize',10);
    %xlim([0 13]);
    %ylim([10^(-7) 1]);
    title('nErrMax=20');
    xlabel('SNRdB')
    ylabel('Pe')
    grid on;
