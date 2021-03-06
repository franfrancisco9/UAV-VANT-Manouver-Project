clear; clc; close all;

%Conversoes
deg = pi/180; %conversao graus em rad
kn=0.514444444; %conversao kt em m/s

%===============INICIALIZACAO DE VARIAVEIS=====================

h = 100;             
aa0 = 8.79*deg;      
gg0 = 0*deg;
g = 9.81;
tt0 = aa0 + gg0;    
u0 = 29.2*kn;        
flaps = 0*deg;
t_sim = 1000;

th0 =58;             
de0 = -7.44*deg; 
da0 = 0.58*deg; 
dr0 = -0.01*deg;
Teng = 0.14;        

demax = 30*deg;
demin = -30*deg;
damax = 30*deg; 
drmax = 30*deg; 
flapmax = 40*deg;
m =22.3 ;
Ix =1.548 ;
Iy =2.841 ;
Iz =3.828 ;
Ixz =0.1 ;


S =0.90 ; 
b =3.000; 
c =0.300; 
aamax =18.00*deg;


xu = -0.1468; 
xw = 0.5034; 
zu = -1.4226; 
zw = -2.1822; 
zwp = -0.0036; 
zq = -1.1256; 
mu = 0.0000; 
mw = -1.0989; 
mq = -1.8976; 
mwp = -0.0829;

ybb = -0.1159; 
lbb = -30.5785; 
nbb = 10.5991; 
yp = 0.0000; 
lp = -12.2557; 
np = -0.9814; 
yr = 0.0000; 
lr = 4.1257; 
nr = -0.2519;

xde = 0.000; 
zde = 0.711; 
mde = -16.356; 
xdf = -0.563; 
zdf = -1.688; 
mdf = -0.661; 
xdt = 3.146; 
zdt = 0.000; 
mdt = 0.000;

Lda = -36.160; 
Nda = 0.000; 
Ydr = -0.006; 
Ldr = 0.000; 
Ndr = 10.894;

%===================PONTO 1: ANALISE MODELO========================
%Matriz da dinamica em Espaco de Estados
A = [ybb,yp+sin(tt0),yr-1,g*cos(tt0)/u0, 0;
    lbb + Ixz/Ix*nbb,lp + Ixz/Ix*np,lr + Ixz/Ix*nr,0, 0; 
    nbb + Ixz/Iz*lbb,np + Ixz/Iz*lp,nr + Ixz/Iz*lr,0, 0;
    0,1,tan(tt0),0, 0;
    0 0 1/(cos(tt0)) 0 0];

B = [0, Ydr;
    Lda + Ixz/Ix*Nda, Ldr + Ixz/Ix*Ndr;
    Nda + Ixz/Iz*Lda, Ndr + Ixz/Iz*Ldr;
    0, 0;
    0, 0];
    
%uu = [de; dt] 

C = eye(5);
D = zeros(5,2);
 
damp(A);

%==============PONTO 2: SAE=========================================
% FUGOIDE
% Realimentacao de u 
%[num,den] = ss2tf(A,B,C,D,2);  
%fun_tranfer = tf(num(1,:),den);
%rlocus(fun_tranfer)
%[k,POLES] = rlocfind(fun_tranfer)

% PERIODO CURTO
% Realimentacao de q 
[num,den] = ss2tf(A,B,C,D,1);  
fun_tranfer = tf(num(3,:),den);
% rlocus(-fun_tranfer)
% [k,POLES] = rlocfind(-fun_tranfer)

% Ganho da Realimentacao de u e de q
k_sae = [0 0 0.4 0 0];

Aaf = A-B(:,2)*k_sae;
damp(Aaf)


%==============PONTO 3: CONTROLO======================================
% SEM INTEGRADOR
% umax = 2.5;
% wmax = 1; 
% qmax = 8*deg;
% ttmax = 3.5*deg;
% hmax =  5;
% 
% demax = 5*deg;   
% dtmax = 0.15;
% 
% Q = diag([1/(umax)^2 1/(wmax)^2 1/(qmax)^2 1/ttmax^2 1/hmax^2]);
% R = diag([1/demax^2 1/dtmax^2]);
% 
% K_lqr = lqr(A,B,Q,R)
% damp(A-B*K_lqr)
% 
% Ky = [K_lqr(1,1), K_lqr(1,5);
%       K_lqr(2,1), K_lqr(2,5)];
% 
%   
% Kc = [0, K_lqr(1,2), K_lqr(1,3), K_lqr(1,4), 0;
%       0, K_lqr(2,2), K_lqr(2,3), K_lqr(2,4), 0];
%   
% sel_u_h = [1 0 0 0 0;
%            0 0 0 0 1];


% COM INTEGRADOR
% umax = 1;
% wmax = 1; 
% qmax = 8*deg;
% ttmax = 3.5*deg;
% hmax = 1;
% demax_lqr = 5*deg;   
% dtmax = 0.15;
% 
% 
% Q = diag([1/(umax)^2 1/(wmax)^2 1/(qmax)^2 1/ttmax^2 1/hmax^2 0.1 0.1]);
% R = diag([1/demax_lqr^2 1/dtmax^2]);



%Com Sensores
% umax = 1;
% wmax = 1; 
% qmax = 8*deg;
% ttmax = 3.5*deg;
% hmax =  5;
% demax_lqr = 1*deg;   
% dtmax = 0.01;
% 
% Q = diag([1/(umax)^2 1/(wmax)^2 1/(qmax)^2 1/ttmax^2 1/hmax^2 0.1 0.01]);
% R = diag([1/demax_lqr^2 1/dtmax^2]);

%Com estimador
max_bb = 15*0.5*deg; 
max_p = 0.5;
max_r = 0.5;
max_phi = 30*0.5*deg;
max_psi = 45*0.5*deg;
max_da = 30*0.5*deg;
max_dr = 30*0.5*deg;

Q11 = 1/(max_bb^2);
Q22 = 1/(max_p^2);
Q33 = 1/(max_r^2);
Q44 = 1/(max_phi^2);
Q55 = 1/(max_psi^2);
R11 = 1/(max_da^2);
R22 = 1/(max_dr^2);

Q = 10* diag([Q11 Q22 Q33 Q44 Q55 5 5]);
R = 5*diag([R11 R22]);

 
A_lqr_int  = [ybb,yp+sin(tt0),yr-1,g*cos(tt0)/u0, 0, 0, 0;
    lbb + Ixz/Ix*nbb,lp + Ixz/Ix*np,lr + Ixz/Ix*nr,0, 0, 0, 0; 
    nbb + Ixz/Iz*lbb,np + Ixz/Iz*lp,nr + Ixz/Iz*lr,0, 0, 0, 0;
    0,1,tan(tt0),0, 0, 0, 0;
    0 0 1/(cos(tt0)) 0 0, 0, 0;
    1 0 0 0 0 0 0;
    0 0 0 1 0 0 0];

B_lqr_int = [0, Ydr;
    Lda + Ixz/Ix*Nda, Ldr + Ixz/Ix*Ndr;
    Nda + Ixz/Iz*Lda, Ndr + Ixz/Iz*Ldr;
    0, 0;
    0, 0;
    0, 0;
     0, 0];

K_lqr = lqr(A_lqr_int, B_lqr_int, Q, R);
damp(A_lqr_int - B_lqr_int*K_lqr)

Ky = [K_lqr(1,1), K_lqr(1,4);
      K_lqr(2,1), K_lqr(2,4)];

  
% Kc = [0, K_lqr(1,2), K_lqr(1,3), K_lqr(1,4), 0;
%       0, K_lqr(2,2), K_lqr(2,3), K_lqr(2,4), 0];
  
Kc_control =[K_lqr(1,2), K_lqr(1,3), K_lqr(1,5);
             K_lqr(2,2), K_lqr(2,3), K_lqr(2,5)];
  
K_int = [K_lqr(1,6), K_lqr(1,7);
         K_lqr(2,6), K_lqr(2,7)];
  
sel_u_h = [1 0 0 0 0 ;
           0 0 0 1 0 ];
       
       
%estimador_kalman       
C_kalman = [0 1 0 0 0;
            0 0 1 0 0;
            0 0 0 1 0;
            0 0 0 0 1];
        
D_kalman = zeros(4,2);



%==============PONTO 4: SENSORES E ATUADORES=========================


atuadores_vmax = 60 * deg; 
atuadores_t = 40 * 0.001; 
atuadores_f = 100; 

pe_upper = 110000 * 100; 
pe_lower = 45000 * 100;  
pe_f = 10;
pe_res = 0.02 * 100; 
pe_rms = 0.36; 


pd_upper = 2 * 1000; 
pd_lower = -2 * 1000;
pd_upper_v = 5; 
pd_lower_v = 0;
pd_t = 10 * 0.001; 

a_upper = 4*g; 
a_lower = -4*g;
a_res = 0.4 * 0.001 * g; 
a_rms = 3 * 0.001 *g; 

% Girosc?pio de Raz?o Angular
g_upper = 300; 
g_lower = -300; 
g_upper_v = 4.3; 
g_lower_v = 0.7; 
g_rms = 4.4; 
g_gain_gyro = (g_upper_v-g_lower_v)/(g_upper-g_lower);
g_gain = g_gain_gyro * (5-0)/(g_upper_v-g_lower_v);
g_offset = 5 - g_upper * g_gain;

% Magnet?metro
m_upper = 8; 
m_lower = -8 ;
m_res = 0.005 ; 
m_rms = 0.015 ; 
mag_lisboa_y = 0.44 ;
mag_lisboa_z = 0.44 ; 

% Sonar
s_upper = 7.5;
s_lower = 0.2;
s_res = 0.01; 
s_f = 10; 

% GPS
gps_f = 5; 
gps_p_res = 0.5; 
gps_p_rms = 2.5;
gps_v_res = 0.01; 
gps_v_rms = 0.1; 

% Conversor A/D
ad_bits = 12; 
ad_upper_v = 5; 
ad_lower_v = 0; 
ad_quantization = (ad_upper_v-ad_lower_v)/(2^(ad_bits) - 1);
ad_rmd = (1.5 * ad_quantization);
ad_f = atuadores_f; 
%==============PONTO 5: SIMULACAO=========================         
         
R_est= diag([(g_rms)^2 (g_rms)^2 (m_rms)^2 (m_rms)^2]);
Q_est= 100* eye(5);

%Vari?veis de Seguimento

global passo;
global vec_anterior;
global j;
global Influencia;  
global Vector_Base;

j=1; passo=2;
Influencia=100/16;
vec_anterior=[0 0];
Vector_Base=[0 0];


%Defini??o do Percurso

global Path

passo=1;
frac = 16;
%AB
x0N=0*u0/frac:passo:400*u0/frac;
x0E=zeros(1, max(size(x0N)));

%BC
r = 100*u0/frac;
x = 100*u0/frac;
y = 400*u0/frac;

theta = linspace(pi, 0, 100);
xCirc = r * cos(theta) + x;
yCirc = r * sin(theta) + y;

%CD-SemiCirc1
r1 = 100*u0/frac;
x1 = 300*u0/frac;
y1 = 400*u0/frac;
n1 = r1*pi; %perimetro calculado de modo a ter um ponto por metro de trajet?ria

theta1 = linspace(-pi, 0, 100); %perimetro aproximadamente = 4720
xCirc1 = r1 * cos(theta1) + x1;
yCirc1 = r1 * sin(theta1) + y1;

%CD-SemiCirc2
r2 = 95*u0/frac;
x2 = 305*u0/frac;
y2 = 400*u0/frac;
n2 = r2*pi*u0;

theta2 = linspace(0, pi, 100);
xCirc2 = r2 * cos(theta2) + x2;
yCirc2 = r2 * sin(theta2) + y2;

%DE
x1N=400*u0/frac:-passo:0*u0/frac;
x1E=zeros(1, max(size(x0N)));
x1E(1,:) = 210*u0/frac;

%EF
r3 = 105*u0/frac;
x3 = 105*u0/frac;
y3 = 0*u0/frac;
n3 = r3*pi;

theta3 = linspace(0, -pi, 100);
xCirc3 = r3 * cos(theta3) + x3;
yCirc3 = r3 * sin(theta3) + y3;

%Percurso completo
Path = cat(1, [x0E' x0N'], cat(1,[xCirc' yCirc']), cat(1,[xCirc1' yCirc1']), cat(1,[xCirc2' yCirc2']), [x1E' x1N'], cat(1,[xCirc3' yCirc3']));
plot(Path(:,1), Path(:,2))

