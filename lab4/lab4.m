clear;
I=imread('lena.bmp');

fftI=fft2(I);
sfftI=fftshift(fftI);
RRfdp=real(sfftI);
IIfdp=imag(sfftI);
a=sqrt(RRfdp.^2+IIfdp.^2);
a=(a-min(min(a)))/(max(max(a))-min(min(a)))*255;
subplot(2,2,1),imshow(real(a));
b=angle(fftI);
subplot(2,2,2),imshow(real(b));
theta=30;
RR1=a*cos(theta);
II1=a*sin(theta);

fftI1=RR1+1i.*II1;
C=ifft2(fftI1)*255;
subplot(2,2,3),imshow(real(C));

mm=150;
RR2=mm*cos(angle(fftI));
II2=mm*sin(angle(fftI));
fftI2=RR2+1i.*II2;
D=ifft2(fftI2);

subplot(2,2,4),imshow(real(D));
subplot(2,2,1),title('幅值谱');
subplot(2,2,2),title('相位谱');
subplot(2,2,3),title('幅值谱重构图像');
subplot(2,2,4),title('相位谱重构图像');
