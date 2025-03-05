data = readmatrix('scope_41.csv','NumHeaderLines',2);

t = data(:,1);
u = data(:,2);
y = data(:,3);
y_z = data(:,4);


 figure(1);
 subplot(3,1,1);
 plot(t,u);xlabel('timp [s]');ylabel('Semnal de intrare u [V]')
 subplot(3,1,2)
 plot(t,y,'g');xlabel('timp [s]');ylabel(' Semnal iesire y [V]')
 subplot(3,1,3)
 plot(t,y_z,'r');xlabel('timp [s]');ylabel(' Semnal iesire y_z [V]')
%%
%identificare folosind frecventa de rezonanta te
wr = pi/(t(445)-t(433))
Ay = (max(y)-min(y))/2
Au = (max(u)-min(u))/2
Mr = Ay/Au
K =mean(y)/mean(u)
val = roots([-4 0 4 0 -1/(Mr)^2])
z = min(val(val>0))
wn  = wr/sqrt(1-2*z^2)


dy=(y(2)-y(1))/(t(2)-t(1))
A = [0,1; ...
    -(wn^2),-2*z*wn];
B = [0;K*wn^2];
C = [1 0];
D = 0;
H = ss(A,B,C,D);


[y_calc,~,~] = lsim(H,u,t,[y(1),dy]);

figure(2) 
plot(t,[y,y_calc])

error1 = norm(y-y_calc)/norm(y-mean(y))*100
%Raspunsul in frecventa a sistemului fara zerouri
%Calculam w-urile si modulul la frecvente joase
%%
w1 = pi/(t(141)-t(100));
w2 = pi/(t(205)-t(176));
w3 = pi/(t(251)-t(229));
w4 = pi/(t(292)-t(273));
w5 = pi/(t(327)-t(311));

M1 = (y(145)-y(99))/(u(141)-u(100));
M2 = (y(206)-y(179))/(u(205)-u(176));
M3 = (y(255)-y(232))/(u(251)-u(229));
M4 = (y(295)-y(276))/(u(292)-u(273));
M5 = (y(330)-y(314))/(u(327)-u(311));

% %Calculam w-uri si modulul la frecvente medii
w6 = pi/(t(388)-t(374));
w7 = pi/(t(416)-t(403));
w8 = pi/(t(489)-t(478));

M6 = (y(392)-y(378))/(u(388)-u(374));
M7 = (y(419)-y(406))/(u(416)-u(403));
M8 = (y(493)-y(481))/(u(489)-u(478));

%Calculam w-uri si modulul la frecvente inalte
w9 = (pi/(t(661)-t(653))+pi/(t(662)-t(653)))/2
w10 = pi/(t(949)-t(943));

M9 = (y(665)-y(658))/(u(661)-u(653));
M10 = (y(954)-y(948))/(u(949)-u(943));

%Calculam faza
ph1 = -w1*(t(100)-t(97));
ph2 =-w2*(t(179)-t(176));
ph3 = -w3*(t(232)-t(229));
ph4 = -w4*(t(276)-t(273));
ph5 = -w5*(t(314)-t(311));
ph6 = -w6*(t(378)-t(374));
ph7 = -w7*(t(406)-t(403));
ph8= -w8*(t(481)-t(478));
phr = -wr*(t(437)-t(433));
ph9 = -w9*(t(658)-t(653));
ph10 = -w10*(t(948)-t(943));


w = [w1 w2 w3 w4  w5  w6  wr   w9 w10]
M = [M1 M2 M3 M4 M5 M6  Mr   M9  M10]
ph = rad2deg([ph1,ph2,ph3,ph4,ph5,ph6,phr,ph9,ph10])

figure(3)
semilogx(w,20*log10(M),'r*') ;grid on;hold on;
bode(H);grid on;hold on;
semilogx(w,ph,'b*');grid on;hold on;

panta = 20*log10(M10/M9)/(w10/w9)*10
%%
%Identificare parametrica ARMAX sistem fara zero
dt = 1.0e-05
datay = iddata(y,u,dt)
datay_z=iddata(y_z,u,dt)

Marmax_y = armax(datay,[2 2 2 0])
Marmax_y = pem(datay,Marmax_y)
[Num,Den] = tfdata(Marmax_y,'v')
sys = idss(Marmax_y)
[A,B,C,D] = tf2ss(Num,Den)

figure(4)
subplot(2,1,1);resid(Marmax_y,datay,'corr',7)
subplot(2,1,2);compare(datay,Marmax_y)

%conditii initiale
x1 = [y(1);0]
x2 = A'*x1+C'*u(1)

y_armax_calc = dlsim(A',C',B',D,u,[2.8125,-2.160]);
figure(5)
plot(t,[y,y_armax_calc]);xlabel('timp [s]');ylabel('Amplitudine [V])');title('Simulare model ARMAX pentru sistemul de ordin 2 fara zero')

H_armax = tf(Num,Den,dt)

error_armax = norm(y_armax_calc-y)/norm(y-mean(y))*100
%%
%Identifacare parametrica iv sistem fara zero
Te = t(2)-t(1);
datay = iddata(y,u,Te);
M_iv = iv4(datay,[2 1 1])

figure(5);
subplot(2,1,1);resid(M_iv,datay,'corr',5);
subplot(2,1,2);compare(M_iv,datay);hold off;

x1 = [y(1);0];
x2 = A'*x1+C'*u(1);

[num,den] = tfdata(M_iv,'v');
H_iv =tf(num,den,Te)
[A,B,C,D] = tf2ss(num,den);
y_sim = dlsim(A',C',B',D,u,[2.8125 -2.1620]);

figure(6)
plot(t,[y,y_sim]);xlabel('timp [s]');ylabel('Amplitudine [V])');title('Simulare model IV pentru sistemul de ordin 2 fara zero')

error_iv = norm(y_sim-y)/norm(y-mean(y))*100
%%
%Identificare parametrica armax pentru sistemul cu zero
datayz = iddata(y_z,u,1.0e-5)
M_armax = armax(datayz,[2 1 1 0])
M_armax = pem(datayz,M_armax)

figure(7);
subplot(2,1,1);resid(M_armax,datayz,'corr',5);
subplot(2,1,2);compare(M_armax,datayz);

[num,den] = tfdata(M_arx,'v');
Hz_armax = tf(num,den,1.0e-05)
[A,B,C,D] = tf2ss(num,den);

x1 = [y(1);0];
x2 = A'*x1+C'*u(1);

y_sim = dlsim(A',C',B',D,u,[2.8125 -2.3412]);
figure(8);plot(t,[y,y_sim]);xlabel('timp [s]');ylabel('Amplitudine [V])');title('Simulare model ARMAX pentru sistemul de ordin 2 cu zero')

error_z_armax = norm(y_sim-y_z)/norm(y_z-mean(y_z))*100
%%
%Identificare parametrica iv sistem cu zero
M_ivz = iv4(datayz,[2 1 0])
figure(9);

subplot(2,1,1);resid(M_ivz,datayz,'corr',5);
subplot(2,1,2);compare(M_ivz,datayz);

[num,den] = tfdata(M_ivz,'v');
Hz_iv= tf(num,den,1.0e-05)
[A B C D] = tf2ss(num,den);

x1 = [y(1);0];
x2 = A'*x1+C'*u(1);

yz_sim = dlsim(A',C',B',D,u,[2.8125 -2.2984]);
figure(10);plot(t,[y_z,yz_sim]);xlabel('timp [s]');ylabel('Amplitudine [V])');title('Simulare model IV pentru sistemul de ordin 2 cu zero')

error_z_iv = norm(yz_sim-y_z)/norm(y_z-mean(y_z))*100