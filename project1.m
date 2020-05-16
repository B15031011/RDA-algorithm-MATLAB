clc;clear all;close all;
%% LSJ
%% ��������
% �״����
C = 2.9979e8;                                                               % ����
f0 = 5.3e9;                                                                 % ����Ƶ��
lambda=C/f0;                                                                % ����
% ����ƽ̨
Vr= 7062;                                                                   % ƽ̨�ٶ�
% ������
Kr = -0.72135e12;                                                           % ��Ƶ��
Tr = 41.75e-6;                                                              % �������ʱ��
Fr = 32.317e6;                                                              % ���������Ƶ��
dt=1/Fr;                                                                    % ���������ʱ����
% ��λ��
Ka = 1733;                                                                  % ��λ���Ƶ��
Ba = 1256.98;                                                               % ��λ�����Ƶ�� PRF
fc = -6900;                                                                 % ����������Ƶ��
ds = 1/Ba;                                                                  % ��λ�����ʱ����
% ��������
t0 = 6.5956e-3;                                                             % ���ݴ���ʼʱ��  %
%load('C:\Users\Dell\Desktop\ң���źŴ���\����ҵ\EXTRACTED_DATA\CDdata1.mat');
load('EXTRACTED_DATA\CDdata1'); 
data1=data;
load('EXTRACTED_DATA\CDdata2');                                            %ѡ����һ�����ݵ������
data2=data;
data=[data1;data2];
data=double(data);
%figure,imagesc(abs(data))
[Nslow,Nfast] = size(data);
Rmin=t0*C/2;   
Rmax=t0*C/2+Nfast/2*dt*C;
%% ȷ���ο�����
%Rref=(t0+Nfast/2*dt)*C/2;                                                   % �ο�����ѡ�ڳ�����������
Rref=(Rmax+Rmin)/2;
%% ����ĩβ����
Za = 800;                                                                   % ��λ������
Zr = ceil(Tr/dt);                                                           % ����������
data = cat(2,data,zeros(Nslow,Zr));                                         % ��������
data = cat(1,zeros(Za,Nfast+Zr),data);                                      % ��λ����
Nslow = Nslow+Za;                                                           %��λ�����
Nfast = Nfast+Zr;                                                           %���������
figure,imagesc(abs(data));axis image;set(gcf,'Color','w');
title('Time domain��raw data after zero-padding');
xlabel('Range direction');ylabel('Azimuth direction��');

%% ʱ���ᡢƵ��������
tr = t0+(Nfast-1)/Fr;                                                       %������ʱ����
Ffast = ((0:Nfast-1)-Nfast/2)/Nfast*Fr;                                     % ������Ƶ����[-Fr/2,Fr/2)
ts = ((0:Nslow-1)-Nslow/2)*ds;                                              % ��λ��ʱ����
Fslow = fc+((0:Nslow-1)-Nslow/2)/Nslow*Ba;                                  % ��λ��Ƶ����

%% ���뵽��λ�����������
S0=data.*exp(-1i*2*pi*fc*(ts'*ones(1,Nfast)));                              % ������������Ƶ������
%% �ڶ�άƵ����ʵ��SRC�;���ѹ��
t=t0:1/Fr:t0+(Nfast-1)/Fr;
xl=t*C/2;                                                                   % x���ע ������
yl=ts*Vr;                                                                   % y���ע ��λ��
t_m=repmat(t,Nslow,1);                                                      %ʱ�����
ts_m=repmat(ts',1,Nfast);
R0=t_m*C/2;
Fslow_m = repmat(Fslow',1,Nfast);                                           % ��λƵ�ʾ���
Ffast_m = repmat(Ffast,Nslow,1);                                            % ����Ƶ�ʾ���
I=ones(Nslow,Nfast);
D=(I-C^2*Fslow_m.^2/(4*Vr^2*f0^2)).^1/2;                                    %��άƵ���е��㶯���� �ο���ʽ��5.30��
Ksrc=2*Vr^2*f0^3*D.^3./(C*R0.*Fslow_m.^2);                                  %�ο���ʽ��6.22����б�ӽ��� ��ʼ�����Ƶ�ʱ仯
%% ����ѹ��/���ξ���ѹ������ѡ��
Refr=exp(1i*pi*Ffast_m.^2/Kr);                                             %������ƥ���˲����������ο���ʽ��3.35�������ڼ򵥵�б�ӽ�
%Refr=exp(1i*pi*Ffast_m.^2/Kr-1i*pi*Ffast_m.^2./Ksrc);                      %���ξ���ѹ����SRC�˲����;���ѹ���˲����ϲ�����˲���,�ο���6.28��
Sr=ifty(fty(S0).*Refr);                                                     %S0������fft+���뷽��ѹ�����ο���6.3��
figure;            
imagesc(xl,yl,abs(Sr));
title('SRC result in time domain');
xlabel('Range direction');ylabel('Azimuth direction');
S1=ftx(Sr);                                                                 %��λ��FFT������ʱ�򣬷�λƵ��
figure;
imagesc(xl,yl,abs(S1));
title('SRC result in range-doppler domain');                                        
xlabel('Range direction');ylabel('Azimuth direction');
%% WRC
%% RCMC �����㶯����
delta_R = lambda^2*R0.*Fslow_m.^2/(8*Vr^2);                                 %�ο���ʽ��6.11��
%delta_R = R0*((1-D)./D)                                                    %��б�ӽǵ���� �ο���ʽ��6.25��

%��һ���ھ�����������о����ֵ���� ÿһ�����б�ࣨR0�������ž����ŵĲ�ͬ���ı䣨����
num_range = C/(2*Fr);                                                       % һ�����������Ԫ����Ӧ�ĳ���
delta_R_num = delta_R./num_range;                                           % ÿһ����λ��Ƶ�ʣ���RCM��Ӧ�ľ��������Ԫ��

R = 8;                                                                      % sinc��ֵ�˳���

for p = 1:Nslow
    for q = 1:Nfast
        delta_R_p = delta_R_num(p,q);
        R_p = q + delta_R_p;
        R_p_zheng = ceil(R_p);                                              %����ȡ��
        ii = (R_p-(R_p_zheng-R/2):-1:R_p-(R_p_zheng+R/2-1));
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);                               %��ֵ��һ��
        %ii��sinc��ֵ���̵ı���
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
        % ����S_rdֻ��������ȡֵ���ҷ�Χ���ޡ���˲�ֵ��Ҫ��������ȡֵ����߽����⡣������ѭ����λ��˼��
        if (R_p_zheng-R/2) > Nfast                                          %ȫ����
            ll = (R_p_zheng-R/2-Nfast:1:R_p_zheng+R/2-1-Nfast);
        else
            if (R_p_zheng+R/2-1) > Nfast                                    %��������
                ll_1 = (R_p_zheng-R/2:1:Nfast);
                ll_2 = (1:1:R_p_zheng+R/2-1-Nfast);
                ll = [ll_1,ll_2];
            else
                if(R_p_zheng+R/2-1) < 1                                     %ȫ���磨�����ܷ�����������Ҫ���ǣ�
                   ll = (R_p_zheng-R/2+Nfast:1:R_p_zheng+R/2-1+Nfast);
                else
                    if (R_p_zheng-R/2) < 1                                  % ��������
                        ll_1 = (R_p_zheng-R/2+Nfast:1:Nfast);
                        ll_2 = (1:1:R_p_zheng+R/2-1);
                        ll = [ll_1,ll_2];
                    else
                        ll = (R_p_zheng-R/2:1:R_p_zheng+R/2-1);
                    end                    
                end
            end
        end   
        rcmc_S = S1(p,ll);
        S2(p,q) = sum( rcmc_sinc.*rcmc_S );
    end
end
        
%%������FFT ������λ��� IFFT
% 
% G = exp(1i*4*pi*Ffast_m.*delta_R/C);                                        %��λ�˷��� ��6.12��
% S2 = S1.*G;
% 
figure;
imagesc(xl,yl,abs(S2));
title('RCMC result in range-doppler domain');                                        
xlabel('Range direction');ylabel('Azimuth direction');


%% JCWN
%% ��λѹ��
Ka = 2*Vr^2/lambda./R0;                                                     %��λ��Ƶ�� ��5.3��
H = exp(-1i*pi*Fslow_m.^2./Ka);                                             %��λ��ƥ���˲��� ��6.16�� 
S3 = S2.*H;                                                                 %��λѹ��(6.17)
s_ac = iftx(S3);                                                            %��λ��IFFT
J=abs(s_ac);
figure;            
image(xl,yl,J);
title('Azimuth compression result');
xlabel('Range direction');ylabel('Azimuth direction');
%% �����任
v=100;
r=mat2gray(double(J));
S=log(1+v*r)/(log(v+1));
S=flipud(S);
%S=imcrop(S,[0,800,2048,1535]);
S=imcrop(S,[0,800,2048,3071]);
figure; 
imshow(S,[]);
xlabel('Logarithmic transformation'); 
%imwrite(S,'4.jpg');

