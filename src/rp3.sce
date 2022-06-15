/*
---------------
--UAV: flight condition: 4
h=100 m; aa0 =8.79 deg; gg0=0 deg; u0=29.2 kn; flaps=0 deg.
inputs: th0 =58(%); de0= -7.44 deg; da0 =0.58 deg; dr0=-0.01 deg;
Teng =0.14 s; demax =+30/- -30 deg; damax =30 deg; drmax =30 deg; flapmax =40 deg
inertial data:
m=22.3 kg; Ix =1.548 kg.m^2; Iy =2.841 kg.m^2; Iz =3.828 kg.m^2; Ixz =0.1 kg.m^2
wing data: S=0.90 m^2; b=3.000 m; c=0.300 m; aamax =18.00 deg
derivatives (no units or SI units):
xu xw zu zw zwp zq mu mw mq mwp
-0.1468 0.5034 -1.4226 -2.1822 -0.0036 -1.1256 0.0000 -1.0989 -1.8976 -0.0829
ybb lbb nbb yp lp np yr lr nr
-0.1159 -30.5785 10.5991 0.0000 -12.2557 -0.9814 0.0000 4.1257 -0.2519
xde zde mde xdf zdf mdf xdt zdt mdt
0.000 0.711 -16.356 -0.563 -1.688 -0.661 3.146 0.000 0.000
Lda Nda Ydr Ldr Ndr
-36.160 0.000 -0.006 0.000 10.894
-----------------
*/

clear;
close all;

// Constantes e conversões
kn = 0.5144444444444444;
deg = %pi/180; 
g = 9.81;

// dados_iniciais
h = 100;             // m 
aa0 = 8.79*deg;      // alfa0_rad
gg0 = 0*deg;         // gama0_rad
tt0 = aa0 + gg0;     // teta0_rad
u0 = 29.2*kn;        // em m/s
flaps = 0*deg;

// inputs : 
th0 =58;             // %
de0 = -7.44*deg; 
da0 = 0.58*deg; 
dr0 = -0.01*deg;
Teng = 0.14;         //s 
demax = 30*deg;
demin = -30*deg;
damax = 30*deg; 
drmax = 30*deg; 
flapmax = 40*deg;

// inertial data :
m =22.3 /*kg*/; 
Ix =1.548 /*kg .m ^2*/; 
Iy =2.841 /*kg .m ^2*/; 
Iz =3.828 /*kg .m ^2*/; 
Ixz =0.1 /*kg .m ^2*/;

// wing data :
S =0.90 /*m ^2*/; 
b =3.000 /*m*/; 
c =0.300 /*m*/; 
aamax =18.00*deg;

// derivadas ( sem unidades ou unidades SI):

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

//==============================================================================
// =============================== RP1 =========================================
//============================== Ponto 1 =======================================
//==============================================================================

// Matriz A
A = [ybb,yp+sin(tt0),yr-1,g*cos(tt0)/u0;
    lbb + Ixz/Ix*nbb,lp + Ixz/Ix*np,lr + Ixz/Ix*nr,0; 
    nbb + Ixz/Iz*lbb,np + Ixz/Iz*lp,nr + Ixz/Iz*lr,0;
    0,1,tan(tt0),0];

// Matriz B
B = [0, Ydr;
    Lda + Ixz/Ix*Nda, Ldr + Ixz/Ix*Ndr;
    Nda + Ixz/Iz*Lda, Ndr + Ixz/Iz*Ldr;
    0, 0];
    
// Matriz C
C = [1, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1];

// Matriz D
D = zeros(4,2);


// Valores Próprios
valores_proprios = spec(A);
[wn,z] = damp(valores_proprios);

// Dados para analisar qualidades de voo
T2 = log(2)/valores_proprios(4);
tau = 1/-valores_proprios(1);
disp(T2, tau)

//==============================================================================
//==============================================================================


//==============================================================================
// =============================== RP2 =========================================
//============================== Ponto 2 =======================================
//==============================================================================

// Função Tranferência 
H = syslin('c', A, B, C, D);
h = ss2tf(H)
disp(h(3,2))

//Rootlocus
//clf();
//evans(h(3, 2),10)
//sgrid('red')

// k retirado do gráfico
k_sae = [0 0 2.849 0];

Aaf = A-B(:,2)*k_sae;
Caf =[0, 0, 0, 0;
0, 0, 0, 0;
0, 0, 1, 0;
0, 0, 0, 0];
valores_proprios_f = spec(Aaf);
[wn_f,z_f] = damp(valores_proprios_f);
tauf = 1/-valores_proprios_f(1);

//disp(valores_proprios, valores_proprios_f)
//disp(wn, wn_f)
//disp(z, z_f) 

//==============================================================================
// =============================== RP2 =========================================
//============================== Ponto 3 =======================================
//==============================================================================

//controlo ótimo
//Método de Bryson

//Valores escolhidos para Q
max_bb = 15*0.5*deg; 
max_p = 0.5;
max_r = 0.5;
max_phi = 30*0.5*deg;
//Valores escolhidos para R
max_da = 30*0.5*deg;
max_dr = 30*0.5*deg;

//Dados da Matriz Q
Q11 = 1/(max_bb^2);
Q22 = 1/(max_p^2);
Q33 = 1/(max_r^2);
Q44 = 1/(max_phi^2);
//Dados da Matriz R
R11 = 1/(max_da^2);
R22 = 1/(max_dr^2);

// Matriz Q e R
Q = diag([Q11 Q22 Q33 Q44]);
R = diag([R11 R22]);

// Ganho K
K1=-lqr(H,Q,R);

// Valores Próprios, w_n e fator de amortecimento
valores_proprios_c=spec(A-B*K1);
[wn_c,z_c] = damp(valores_proprios_c);
disp(valores_proprios, valores_proprios_f, valores_proprios_c, K1)
disp(wn, wn_c)
disp(z, z_c)

A1 = A-B*K1;
C1 = [1, 0, 0, 0;
    0, 0, 0, 1;]
D1 = zeros(2,2);
G = -C1*inv(A1)*B;
//disp(G)
F = inv(G); // ganho estático para seguimento de referência

//gráfico
bb_ref = 0;
phi_ref = max_phi;
x0 = zeros(4,1);

select_bb_phi = [1, 0, 0, 0;
                0, 0, 0, 1];
Ky = [K1(1,1), K1(1,4);
    K1(2,1), K1(2,4)];
Kc = [0, K1(1,2), K1(1,3), 0;
    0, K1(2,2), K1(2,3),  0];

// xcos("LQR.zcos")

// H1 = syslin('c', A1, B, C1, D1);
// h1 = ss2tf(H1)
// disp(h1)
// t = 0: 0.01 : 10;
// y1 = csim('step', t, h1(1,1));
// figure
// plot (t, y1)

//Integrativo
Q =  diag([Q11 Q22 Q33 Q44 5 5]);
// Sensore:
Q =  diag([Q11*1000 Q22*10 Q33*10 Q44*1000 500 50]);
R = diag([R11*10000 R22*100000])
A_int = [ybb,yp+sin(tt0),yr-1,g*cos(tt0)/u0, 0, 0;
    lbb + Ixz/Ix*nbb,lp + Ixz/Ix*np,lr + Ixz/Ix*nr,0, 0, 0; 
    nbb + Ixz/Iz*lbb,np + Ixz/Iz*lp,nr + Ixz/Iz*lr,0, 0, 0;
    0,1,tan(tt0),0, 0, 0;
    1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0];

// Matriz B com integrativos
B_int = [0, Ydr;
    Lda + Ixz/Ix*Nda, Ldr + Ixz/Ix*Ndr;
    Nda + Ixz/Iz*Lda, Ndr + Ixz/Iz*Ldr;
    0, 0;
    0, 0;
    0, 0];
C_int = [1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0
    0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1];
D_int =  zeros(6,2);
H_int = syslin('c', A_int, B_int, C_int, D_int);
K_int = -lqr(H_int , Q, R)
disp(spec(A_int - B_int *K_int))
valores_proprios_int=spec(A_int-B_int*K_int);
[wn_int,z_int] = damp(valores_proprios_int);
disp(valores_proprios, valores_proprios_f, valores_proprios_int, K1)
disp(wn, wn_int)
disp(z, z_int)

Ky = [K_int(1,1), K_int(1,4);
        K_int(2,1), K_int(2,4)];

    
Kc = [0, K_int(1,2), K_int(1,3), 0;
      0, K_int(2,2), K_int(2,3), 0];


    
Kint = [K_int(1,5), K_int(1,6);
            K_int(2,5), K_int(2,6)];



    
// xcos("int.zcos")    

//==============================================================================
// =============================== RP3 =========================================
//============================== Ponto 4 =======================================
//==============================================================================
// Atuadores
atuadores_vmax = 60 * deg; // rad/s
atuadores_t = 40 * 0.001; // s
atuadores_f = 100; // Hz

// Pressao Estática
pe_upper = 110000 * 100; // Pa
pe_lower = 45000 * 100;  // Pa
pe_f = 10; // Hz
pe_res = 0.02 * 100; // Pa
pe_rms = 0.36; // m

// Pressao Dinâmica
pd_upper = 2 * 1000; // Pa
pd_lower = -2 * 1000;// Pa
pd_upper_v = 5; // V
pd_lower_v = 0; // V
pd_t = 10 * 0.001; // s

// Aceleração
a_upper = 4*g; // m/s^2
a_lower = -4*g;// m/s^2
a_res = 0.4 * 0.001 * g; // m/s^2
a_rms = 3 * 0.001 *g; // m/s^2

// Giroscópio de Razão Angular
g_upper = 300;  // º/s
g_lower = -300; // º/s
g_upper_v = 4.3; // V
g_lower_v = 0.7; // V
g_rms = 4.4; // º/s
g_gain = (g_upper_v-g_lower_v)/(g_upper-g_lower);
g_offset = 5 - g_upper * g_gain;

// Magnetómetro
m_upper = 8; // gauss
m_lower = -8; // gauss
m_res = 0.005; // gauss
m_rms = 0.015; // gauss

// Sonar
s_upper = 7.5; // m
s_lower = 0.2; // m
s_res = 0.01; // m
s_f = 10; // Hz

// GPS
gps_f = 5; // Hz
gps_p_res = 0.5; // m
gps_p_rms = 2.5; // m
gps_v_res = 0.01; // m/s
gps_v_rms = 0.1; // m/s

// Conversor A/D
ad_bits = 12; // bits
ad_upper_v = 5; // V
ad_lower_v = 0; // V
ad_quantization = (ad_upper_v-ad_lower_v)/(2^(ad_bits) - 1);
ad_rmd = 1.5 * ad_quantization;
ad_f = atuadores_f; // Hz
