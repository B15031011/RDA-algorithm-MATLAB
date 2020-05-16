clc;clear all;close all;
%% LSJ
%% 参数设置
% 雷达参数
C = 2.9979e8;                                                               % 光速
f0 = 5.3e9;                                                                 % 中心频率
lambda=C/f0;                                                                % 波长
% 飞行平台
Vr= 7062;                                                                   % 平台速度
% 距离向
Kr = -0.72135e12;                                                           % 调频率
Tr = 41.75e-6;                                                              % 脉冲持续时间
Fr = 32.317e6;                                                              % 距离向采样频率
dt=1/Fr;                                                                    % 距离向采样时间间隔
% 方位向
Ka = 1733;                                                                  % 方位向调频率
Ba = 1256.98;                                                               % 方位向采样频率 PRF
fc = -6900;                                                                 % 多普勒中心频率
ds = 1/Ba;                                                                  % 方位向采样时间间隔
% 其他参数
t0 = 6.5956e-3;                                                             % 数据窗开始时间  %
%load('C:\Users\Dell\Desktop\遥感信号处理\大作业\EXTRACTED_DATA\CDdata1.mat');
load('EXTRACTED_DATA\CDdata1'); 
data1=data;
load('EXTRACTED_DATA\CDdata2');                                            %选其中一个数据导入就行
data2=data;
data=[data1;data2];
data=double(data);
%figure,imagesc(abs(data))
[Nslow,Nfast] = size(data);
Rmin=t0*C/2;   
Rmax=t0*C/2+Nfast/2*dt*C;
%% 确定参考距离
%Rref=(t0+Nfast/2*dt)*C/2;                                                   % 参考距离选在成像区域中心
Rref=(Rmax+Rmin)/2;
%% 数据末尾补零
Za = 800;                                                                   % 方位向补零数
Zr = ceil(Tr/dt);                                                           % 距离向补零数
data = cat(2,data,zeros(Nslow,Zr));                                         % 距离向补零
data = cat(1,zeros(Za,Nfast+Zr),data);                                      % 方位向补零
Nslow = Nslow+Za;                                                           %方位向点数
Nfast = Nfast+Zr;                                                           %距离向点数
figure,imagesc(abs(data));axis image;set(gcf,'Color','w');
title('Time domain：raw data after zero-padding');
xlabel('Range direction');ylabel('Azimuth direction）');

%% 时间轴、频率轴设置
tr = t0+(Nfast-1)/Fr;                                                       %距离向时间轴
Ffast = ((0:Nfast-1)-Nfast/2)/Nfast*Fr;                                     % 距离向频率轴[-Fr/2,Fr/2)
ts = ((0:Nslow-1)-Nslow/2)*ds;                                              % 方位向时间轴
Fslow = fc+((0:Nslow-1)-Nslow/2)/Nslow*Ba;                                  % 方位向频率轴

%% 对齐到方位向多普勒中心
S0=data.*exp(-1i*2*pi*fc*(ts'*ones(1,Nfast)));                              % 搬移至多普勒频率中心
%% 在二维频域中实现SRC和距离压缩
t=t0:1/Fr:t0+(Nfast-1)/Fr;
xl=t*C/2;                                                                   % x轴标注 距离向
yl=ts*Vr;                                                                   % y轴标注 方位向
t_m=repmat(t,Nslow,1);                                                      %时间矩阵
ts_m=repmat(ts',1,Nfast);
R0=t_m*C/2;
Fslow_m = repmat(Fslow',1,Nfast);                                           % 方位频率矩阵
Ffast_m = repmat(Ffast,Nslow,1);                                            % 距离频率矩阵
I=ones(Nslow,Nfast);
D=(I-C^2*Fslow_m.^2/(4*Vr^2*f0^2)).^1/2;                                    %二维频域中的徙动因子 参考公式（5.30）
Ksrc=2*Vr^2*f0^3*D.^3./(C*R0.*Fslow_m.^2);                                  %参考公式（6.22）大斜视角下 初始距离调频率变化
%% 距离压缩/二次距离压缩（可选择）
Refr=exp(1i*pi*Ffast_m.^2/Kr);                                             %距离向匹配滤波器函数，参考公式（3.35）适用于简单低斜视角
%Refr=exp(1i*pi*Ffast_m.^2/Kr-1i*pi*Ffast_m.^2./Ksrc);                      %二次距离压缩，SRC滤波器和距离压缩滤波器合并后的滤波器,参考（6.28）
Sr=ifty(fty(S0).*Refr);                                                     %S0距离向fft+距离方向压缩，参考（6.3）
figure;            
imagesc(xl,yl,abs(Sr));
title('SRC result in time domain');
xlabel('Range direction');ylabel('Azimuth direction');
S1=ftx(Sr);                                                                 %方位向FFT，距离时域，方位频域
figure;
imagesc(xl,yl,abs(S1));
title('SRC result in range-doppler domain');                                        
xlabel('Range direction');ylabel('Azimuth direction');
%% WRC
%% RCMC 距离徙动矫正
delta_R = lambda^2*R0.*Fslow_m.^2/(8*Vr^2);                                 %参考公式（6.11）
%delta_R = R0*((1-D)./D)                                                    %大斜视角的情况 参考公式（6.25）

%法一：在距离多普勒域中距离插值运算 每一个最近斜距（R0）都随着距离门的不同而改变（？）
num_range = C/(2*Fr);                                                       % 一个距离采样单元，对应的长度
delta_R_num = delta_R./num_range;                                           % 每一个方位向频率，其RCM对应的距离采样单元数

R = 8;                                                                      % sinc插值核长度

for p = 1:Nslow
    for q = 1:Nfast
        delta_R_p = delta_R_num(p,q);
        R_p = q + delta_R_p;
        R_p_zheng = ceil(R_p);                                              %向上取整
        ii = (R_p-(R_p_zheng-R/2):-1:R_p-(R_p_zheng+R/2-1));
        rcmc_sinc = sinc(ii);
        rcmc_sinc = rcmc_sinc/sum(rcmc_sinc);                               %插值归一化
        %ii是sinc插值过程的变量
        % g(x)=sum(h(ii)*g_d(x-ii)) = sum(h(ii)*g_d(ll));
        % 由于S_rd只有整数点取值，且范围有限。因此插值中要考虑它的取值溢出边界问题。——用循环移位的思想
        if (R_p_zheng-R/2) > Nfast                                          %全右溢
            ll = (R_p_zheng-R/2-Nfast:1:R_p_zheng+R/2-1-Nfast);
        else
            if (R_p_zheng+R/2-1) > Nfast                                    %部分右溢
                ll_1 = (R_p_zheng-R/2:1:Nfast);
                ll_2 = (1:1:R_p_zheng+R/2-1-Nfast);
                ll = [ll_1,ll_2];
            else
                if(R_p_zheng+R/2-1) < 1                                     %全左溢（不可能发生，但还是要考虑）
                   ll = (R_p_zheng-R/2+Nfast:1:R_p_zheng+R/2-1+Nfast);
                else
                    if (R_p_zheng-R/2) < 1                                  % 部分左溢
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
        
%%法二：FFT 线性相位相乘 IFFT
% 
% G = exp(1i*4*pi*Ffast_m.*delta_R/C);                                        %相位乘法器 （6.12）
% S2 = S1.*G;
% 
figure;
imagesc(xl,yl,abs(S2));
title('RCMC result in range-doppler domain');                                        
xlabel('Range direction');ylabel('Azimuth direction');


%% JCWN
%% 方位压缩
Ka = 2*Vr^2/lambda./R0;                                                     %方位调频率 （5.3）
H = exp(-1i*pi*Fslow_m.^2./Ka);                                             %方位向匹配滤波器 （6.16） 
S3 = S2.*H;                                                                 %方位压缩(6.17)
s_ac = iftx(S3);                                                            %方位向IFFT
J=abs(s_ac);
figure;            
image(xl,yl,J);
title('Azimuth compression result');
xlabel('Range direction');ylabel('Azimuth direction');
%% 对数变换
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

