clear all;
close all;

Bw = 1e6; %Bandwidth of the lowpass filter
F0 = 1.5e6;%Bandpass filter's center frequency
Fc = 0.5*Bw;%cutoff frequency of the initial lowpass filter prototype

[B,A] = ellip(6,2,40,2*pi*Fc,'s');%Matlab's built-in elliptical filter design
                                  % 6thorder, 2 dB passband ripple, 40 db rejection

p = roots(A);%Find roots of output-side polynomial (the poles)
pp = p + 1j*2*pi*F0;%Replicate copies of the poles shifted up and down by the 
                    % %bandpass filter's center frequency       
pn = p -1j*2*pi*F0;
p = [pp; pn];

z = roots(B);%Find roots of input-side polynomial (the zeros)
zp = z + 1j*2*pi*F0;%Replicate copies of the zeros shifted up and down by the 
                    %bandpass filter's center frequency
zn = z -1j*2*pi*F0;
z = [zp; zn]; 

Fs = 10*max(imag(p));%10 times larger than largest imaginary part of any pole
L = 1e4;
dF = Fs/L;
lenFreq = L*dF;
F = -0.5*Fs:dF:0.5*Fs;
Ts = 1/Fs;

s = j*2*pi*F;
Z = prod(s-z);
P = prod(s-p);
H = Z./P;

figure(1);
grid on;
plot(F,H);
axis([-10^7 10^7 -80 100]);
xlabel('Frequency');
ylabel('Frequency Response');
title('System Frequency Response with center frequency 1.5MHz');

H_ = fftshift(H);
h_ = ifft(H_); %do inverse fft
h = Fs.*h_; %mucltiply result by Fs; result = h(t)
h_real = real(h); %take real part of h(t)

Tk = L*Ts;%This is a given
t = 0:Ts:Tk; 
figure(2);
grid on;
plot(t,h_real);
axis([0 4e-5 -1.5e8 1.5e8]);
xlabel('Frequency');
ylabel('Frequency Response');
title('System Frequency Impulse Response');

t_u = (0:Ts:2*Tk);% ' is very important to unit step function
impulse = t_u==0;
unitstep = t_u>=0;
Fo = 1.5e6;%Fo is not F0
x = exp(j*2*pi*Fo*t_u);
figure(3);
%subplot(2,1,1);
plot(t_u,x);
title('System Frequency Input Response')
%subplot(2,1,2);
%plot(t_u,x);%This does start at 0
%axis([-Tk Tk -2 2]);

%y = conv(1,h_real);
sf0 = j*2*pi*F0;
Zf0 = prod(sf0-z);
Pf0 = prod(sf0-p);
Hf0 = Zf0./Pf0;
y = Hf0.*x;
%y = conv(1,H);
yss = conv(x,h_real);
%yss = conv(x,H);
t_1 = linspace(0,Tk,length(y));
t_x = linspace(0,Tk,length(yss));
figure(4);
subplot(2,1,1);
plot(t_1,real(y),'r');
plot(t_x,real(yss),'b');
subplot(2,1,2);
plot(t_1,imag(y),'r');
plot(t_x,imag(yss),'b');
%xlabel('');
figure(5);subplot(2,1,1);plot(t_1,real(y),'r');subplot(2,1,2);plot(t_x,real(yss),'b');
