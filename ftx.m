
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward FFT w.r.t. the first variable方位向fft %
%fft(X)是把每个列作为向量，得到每列元素的傅里叶变换
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fs=ftx(s)
 fs=fftshift(fft(fftshift(s)));


