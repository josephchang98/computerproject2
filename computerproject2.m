Bw = 1e6; %Bandwidth of the lowpass filter
F0 = 1.5e6;%Bandpass filter's center frequency
Fc = 0.5*Bw;%cutoff frequency of the initial lowpass filter prototype
%t = 0:1e9;
%x = exp(j*2*pi*F0*t);%input signal 
%y(t)=x(t)*h(t)

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

%H = z/p; %H(f)
%don't know how many values we want in frequency array
F_f = linspace(-1e10, 1e10, 10000);
H = ones(1, length(F_f));
% equation 2 in pdf
for freq = F_f
    for num = 1:length(z)
        H(num) = H(num)*(((1j*2*pi*freq)-z(num))/((1j*2*pi*freq)-p(num)));
    end
end

figure(1);
plot(F_f,H);
xlabel('Frequency');

dF = abs(F_f(1)-F_f(2));
L = length(F_f);
lenFreq = L*dF;
Fs = 10*max(imag(p)); %10 times larger than largest imaginary part of any pole
Ts = 1/Fs;

F1 = -0.5*Fs:dF:-0.5*Fs+.5*lenFreq*dF;
F2 = 0.5*Fs-.5*lenFreq*dF:dF:0.5*Fs-dF;
Freq = [F1,F2];
%Freq = [-0.5*Fs (-0.5Fs + dF) (-0.5Fs + 2dF) (-0.5Fs + .5*lenFreq*dF) (0.5*Fs - .5*lenFreq*dF) (0.5*Fs - 2*dF) (0.5*Fs - dF)]; 
%do inverse fftshift on H
H_ = fftshift(H);
h_ = infft(H'); %do inverse fft
h = Fs.*h_; %multiply result by Fs; result = h(t)
h_real = real(h); %take real part of h(t)

Tk = L*Ts;
t = 0:Ts:Tk; %need to check if starting value is 0
figure(2)
plot(t, h_real);

%not sure if t_u is what we need for x; i think t is what we use.
t_u = (-5:0.01:5)';% ' is very important to unit step function
impulse = t_u==0;
unitstep = t_u>=0;
%ramp = t_u.*unitstep;
%quad = t_u.^2.*unitstep;
x = exp(j*2*pi*F0*t_u);%F0 here needs to be a different variable
x_u = x.*unitstep;
figure(3);
%plot(t_u,[impulse unitstep ramp quad]);
subplot(2,1,1);
plot(t_u,real(x_u));
subplot(2,1,2);
plot(t_u,imag(x_u));

Hfo = 1;
for num = 1:length(z)
    Hfo = Hfo*(((1j*2*pi*F0)-z(num))/((1j*2*pi*F0)-p(num)));
end

y = Hfo.*x_u;

%again, don't think t_u is the correct time array
figure(4);
subplot(2,1,1);
plot(t_u,real(y));
subplot(2,1,2);
plot(t_u,imag(y));