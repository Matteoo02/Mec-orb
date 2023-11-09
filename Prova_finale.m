%Prova finale


r = [-7663.5213,-6485.4986,-2201.1930];
v = [3.5150,-2.9160,-3.8140];
mu = 398600;
 [a,e,in,OM,w,theta] = car2kep(r,v,mu);
theta_tot = [0:pi/180:2*pi];

rp1 = a*(1-e)
ra1 = a*(1+e)

%  [Xp,Yp,Zp,X,Y,Z] = plotOrbit(a,e,deg2rad(in),deg2rad(OM),deg2rad(w),theta_tot);
% Terra3d;
% plot3(X,Y,Z);
% hold on

%dati orbita finale:

a2 = 13200;
e2 = 0.3860;
i2 = 1.4840;
OM2 = 2.7570;
w2 = 0.9111;
theta2 = (0.2903);

rp2 = a2*(1-e2)
ra2 = a2*(1+e)


% [Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,i2,OM2,w2,theta_tot);
% Terra3d;
% plot3(X2,Y2,Z2);



%Passaggio 1: Cambio forma %figure 1
%Passaggio 2: cambio w %figure 2
%Passaggio 3: cambio piano %figure 3


%% Disegno cambio forma con w uguale


r = [-7663.5213,-6485.4986,-2201.1930];
v = [3.5150,-2.9160,-3.8140];
mu = 398600;
 [a,e,in,OM,w,theta] = car2kep(r,v,mu);
theta_tot = [0:pi/180:2*pi];

 [Xp,Yp,Zp,X,Y,Z] = plotOrbit(a,e,deg2rad(in),deg2rad(OM),deg2rad(w),theta_tot);
Terra3d;
plot3(X,Y,Z);
hold on
h=plot3(nan,nan,nan,'or');

a2 = 13200;
e2 = 0.3860;
    set(h,'XData',X(177),'YData',Y(177),'ZData',Z(177));
    drawnow


[Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,deg2rad(in),deg2rad(OM),deg2rad(w),theta_tot);
Terra3d;
plot3(X2,Y2,Z2);

% Disegno cambio forma con w diverso

% figure;
% a = 13200;
% e = 0.3860;
% 
% [Xp,Yp,Zp,X,Y,Z] = plotOrbit(a,e,deg2rad(in),deg2rad(OM),deg2rad(w),theta_tot);
% Terra3d;
% plot3(X,Y,Z);
% hold on
% 
% w4 = 0.9111;
% 
% 
% [Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,deg2rad(in),deg2rad(OM),w4,theta_tot);
% Terra3d;
% plot3(X2,Y2,Z2);

% cambio piano

figure;

[Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,deg2rad(in),deg2rad(OM),deg2rad(w),theta_tot);
Terra3d;
plot3(X2,Y2,Z2);


i6 = 1.4840;
OM6 = 2.7570;

% w2 da calcoli su quaderno
[Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,i6,OM,deg2rad(w),theta_tot);
Terra3d;
plot3(X2,Y2,Z2);
hold on


i6 = 1.4840;
OM6 = 2.7570;
% parametri cambio piano, caso D_i > 0, D_OM > 0
i_i = deg2rad(42.4233); i_f = i6;    D_i = i_f-i_i;
OM_i = deg2rad(OM); OM_f = OM6;  D_OM = OM_f-OM_i; 
w = deg2rad(w);
alpha = acos(cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(D_OM));
% Calcolo u_i
cos_ui = (cos(alpha)*cos(i_i)-cos(i_f))/(sin(alpha)*sin(i_i));
sin_ui = sin(i_f)*sin(D_OM)/sin(alpha);
u_i = atan2(sin_ui,cos_ui);
% Calcolo u_f
cos_uf = (cos(i_i)-cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
sin_uf = sin(i_i)*sin(D_OM)/sin(alpha);
u_f = atan2(sin_uf,cos_uf);
% Calcolo parametri orbite
th_man_1 = u_i-w;
th_man_2 = th_man_1;
w2 = u_f - th_man_2 + 2*pi
figure()
[Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,(i_i),(OM_i),w,theta_tot);
Terra3d;
plot3(X2,Y2,Z2);

[Xp2,Yp2,Zp2,X2,Y2,Z2] = plotOrbit(a2,e2,(i_f),(OM_f),w2,theta_tot);
plot3(X2,Y2,Z2);


%prova mia dei dati (mi torna uguale al quaderno ma non al risultato
%precedente)

u1_2 = asin(sin(pi-i_f)*(sin(D_OM)/sin(alpha)));
theta1_2 = u1_2-w;
u2_2 = asin(sin(i_i)*(sin(D_OM)/sin(alpha)));
theta2_2 = theta1_2;
w2_2 = (u2_2-theta2_2)+2*pi
% hold on
% [Xp4,Yp4,Zp4,X4,Y4,Z4] = plotOrbit(a2,e2,(i_f),(OM_f),w2_2,theta_tot);
% plot3(X4,Y4,Z4);

%% Calcolo Delta V

%Dati

r = [-7663.5213,-6485.4986,-2201.1930];
v = [3.5150,-2.9160,-3.8140];
mu = 398600;
[a,e,in,OM,w,theta] = car2kep(r,v,mu);

rp1 = a*(1-e);
ra1 = a*(1+e);

%dati orbita finale:

a2 = 13200;
e2 = 0.3860;
i2 = 1.4840;
OM2 = 2.7570;
w2 = 0.9111;
theta2 = (0.2903);

rp2 = a2*(1-e2);
ra2 = a2*(1+e);

rat = ra1;
rpt = rp2;
at = (rat+rpt)/2;
et = (rat-rpt)/(rat+rpt);
% manovra cambio forma


Dv1 = sqrt(2*mu*(1/ra1-1/(2*at)))-sqrt(2*mu*(1/ra1-1/(2*a)));
Dv2 = sqrt(2*mu*(1/rp2-1/(2*a2)))-sqrt(2*mu*(1/rp2-1/(2*at)));

% manovra cambio piano

i_i = deg2rad(in); i_f = i2;    D_i = i_f-i_i;
OM_i = deg2rad(OM); OM_f = OM2;  D_OM = OM_f-OM_i; 
w_r = deg2rad(w);
alpha = acos(cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(D_OM));
% Calcolo u_i
cos_ui = (cos(alpha)*cos(i_i)-cos(i_f))/(sin(alpha)*sin(i_i));
sin_ui = sin(i_f)*sin(D_OM)/sin(alpha);
u_i = atan2(sin_ui,cos_ui);
% Calcolo u_f
cos_uf = (cos(i_i)-cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
sin_uf = sin(i_i)*sin(D_OM)/sin(alpha);
u_f = atan2(sin_uf,cos_uf);
% Calcolo parametri orbite
th_man_1 = u_i-w_r;
th_man_2 = th_man_1;
w_plane = u_f - th_man_2 + 2*pi;

p2 = a2*(1-e2^2);
V_th = sqrt(mu/p2)*(1+e2*cos(th_man_1));
Dv_plane = 2*V_th*sin(alpha/2);

%cambio anomalia pericentro
deltaw = w2-w_plane+2*pi;
Dv_per = 2*sqrt(mu/p2)*e2*sin(deltaw/2);
theta_A = deltaw/2
theta_B = 2*pi-deltaw/2

%dv totale
Dv_tot = Dv1+Dv2+Dv_plane+Dv_per

%Calcoli tempi:
t1 = timeOfFlight(a,e,deg2rad(theta),pi,mu)
t_t = pi*sqrt(at^3/mu)
t2 = timeOfFlight(a2,e2,0,th_man_1,mu)
t3 = timeOfFlight(a2,e2,th_man_1,theta_B,mu)


delta_t = t1+t_t+t2+t3
heur = delta_t/3600








