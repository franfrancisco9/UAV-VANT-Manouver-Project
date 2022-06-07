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

// Matriz A
// Matriz A
A = [ybb,yp+sin(tt0),yr-1,g*cos(tt0)/u0, 0;
lbb + Ixz/Ix*nbb,lp + Ixz/Ix*np,lr + Ixz/Ix*nr,0, 0; 
nbb + Ixz/Iz*lbb,np + Ixz/Iz*lp,nr + Ixz/Iz*lr,0, 0;
0,1,tan(tt0),0, 0;
0,0,1/cos(tt0),0,0];
// Cálculo dos valores para o relatório
valores_proprios = spec(A);
[wn,z] = damp(valores_proprios);
T2 = log(2)/valores_proprios(5);
tau = 1/-valores_proprios(4);

re = real(valores_proprios(5,1))
if re < 0
    disp("Espiral estável")
    disp(valores_proprios(5,1))
elseif T2 < 12 then
    disp("T2 Espiral menor que 12");
    disp(T2);
end

//modo rolamento
Tp = 1 / wn(2,1) // tem que ser menor 1
if 1 < Tp then
    disp("Tp Rolamento maior que 1");
    disp(Tp);
else
    disp("Tp Rolamento menor que 1");
    disp(Tp);
end


//modo rolamento holandes
zz_rh = z(3,1) // maior do que 0.19 (imposto pelo projeto maior que 0.6)
if zz_rh < 0.6 then
    disp("Amortecimento RH menor que 0.6");
    disp(zz_rh);
else
    disp("Amortecimento RH maior que 0.6");
    disp(zz_rh);
end

wn_rh = wn(3,1) // maior do que 0.5
if wn_rh < 1 then
    disp("Frequencia RH menor que 1");
    disp(wn_rh);
else
    disp("Frequencia RH maior que 1");
    disp(wn_rh);
end

zz_wn_rh = zz_rh * wn_rh // maior do que 0.35
if zz_wn_rh < 0.35 then
    disp("Freq*Amort RH menor que 0.35");
    disp(zz_wn_rh);
else
    disp("Freq*Amort RH maior que 0.35");
    disp(zz_wn_rh);
end
