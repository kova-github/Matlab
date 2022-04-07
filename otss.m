clear all;
close all;
clc;
%�������� ������ ��� 8 ��������
q = 8;
T = 2;%���
f0 = 20;%���
f = f0 + ((0:q-1)/T);%����� ����������� ������
dt = (1/max(f))/10;
t = 0:dt:T;
E = 10;%������� ��������� �������

%������� �������� si(t),���� 0<t<T
si = zeros(q,length(t));
for i = 1:q
    si(i,:) = sqrt(2*E/T)*cos(2 * pi * f(i).* t);
end
%������ �������� ������� � sin � cos ������������ �������� �� �������
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

%��������� ��������� ������/���
SNRdB = 1:15; %����� �������� ��������� ������/���
SNRdBtheor = 1:13;
SNRtheor = 10.^(SNRdBtheor/10);
 
nErrMax = 20;%��������� ������������ ����� ������ nErrMax � �������� 20..50..100

%���� �� ��������� ��������� ������/��� � ��
for nSNR = 1:length(SNRdB)
    nErr = 0;%��������� �������� �������� ������
    ntest = 0;%��������� �������� �������� ���������
    SNR = 10^(SNRdB(nSNR)/10);%Y
    sigma2=(sum(sum(si.^2)))/(q*2*SNR);
    sigma=sqrt(sigma2);
    
 %���� ������������� ��� ����� �������� ��������� ����� ���
    while nErr < nErrMax
        %������������� ����������� � ������
        %i = floor(q*rand)+1;%�������� �������� i �� ��������� 0,1...,q-1
        %����� ������� rand �������� �� ���������� �������������� ��������� �����
        X = q*rand + 1;
        i = floor(X);
        
        %�������� ��������� ���� ���� ����� �������� ����� ������
        tetta = rand;
        
        %��������� ������ �� ������ ������
        r = cos(2*pi*tetta)*s_cos(i,:) + sin(2*pi*tetta)*s_sin(i,:) + sigma*randn(1,length(t));
        
        %������������� ���������
        %��� i=0,1,..,q-1 ��������� �������� r_c(i) � r_s(i)
        for a = 1:(q)
            r_c(a) = sum(r.*f_cos(a,:))*dt;
            r_s(a) = sum(r.*f_sin(a,:))*dt;
        end
       
       %��������� �������
        [index, i_result] = max((r_s.^2) + (r_c.^ 2));
        
        %�������� ����������
        if (i_result)~= i
             nErr = nErr + 1;%����������� ������� ����� ���������
             clc
disp([SNRdB(nSNR), ntest, nErr, nErrMax]);%���������� ������ �� ������
% disp([SNRdB(nSNR),nErr,nErrMax,ntest]);

        end 
         ntest = ntest + 1;
    end %����� ����� ��� ����� �������� ������/���
    
    %���������� ����������������� ����������� ������
    Pe_experement(nSNR) = nErr/ntest;%������� ����� ������������ ������ �� ����� ���������
end %����� ����� �� ��������� ��������� ������/��� � ��

%��������� ������������� �������� ����������� ������
Pe_theory = zeros(1, length(SNRdBtheor));
for l=1:q+1
   % ���������� ������ �� ���������
   %1 ��������� ������������ ����������� �
   koef = nchoosek(q+1, l);
   %2 �������� ����������� �� (-1) � ������� (l+1)
   mult = koef*(-1)^(l+1);
   %3 �������� ������������ �������� �� ��������� 1/1+l
   mult2=mult*(1/(1+l));
   %4 �������� �� ���������� � ���������� ((-1/1+l)*(E/N0))
   mult3 = mult2*exp(-l/(l+1)*SNRtheor);
   % ������� ����� �� ���������� �� 0 �� q-1
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
