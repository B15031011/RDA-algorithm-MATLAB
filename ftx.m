
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward FFT w.r.t. the first variable��λ��fft %
%fft(X)�ǰ�ÿ������Ϊ�������õ�ÿ��Ԫ�صĸ���Ҷ�任
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fs=ftx(s)
 fs=fftshift(fft(fftshift(s)));


