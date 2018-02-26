clc;
clear;

syms xs ys zs real;
syms r1 r2 r3 r4 r5 r6 r7 r8 r9 real;   % 旋转矩阵
syms x y z real;
syms theta psi real;

%{    未知数    %}
syms dx dy dz real;      % 相机到imu的偏移量
syms itR itP itH real;   % 相机与imu之间的安置角

R = [r1 r2 r3; r4 r5 r6; r7 r8 r9];
bx = x - xs;    by = y - ys;    bz= z- zs;

Rrph = ang2Rotation(itR, itP, itH) * R; % (R * Rx)-1
Vxyz = Rrph * ([bx by bz]' - R' * [dx, dy, dz]');

fx = atan(Vxyz(1, 1)/Vxyz(2, 1)) - theta;
fy = atan(sqrt(Vxyz(2, 1)^2 + Vxyz(1, 1)^2)/Vxyz(3, 1)) - psi;

x11 = diff(fx, dx);
x12 = diff(fx, dy);
x13 = diff(fx, dz);
x14 = diff(fx, itR);
x15 = diff(fx, itP);
x16 = diff(fx, itH);
y11 = diff(fy, dx);
y12 = diff(fy, dy);
y13 = diff(fy, dz);
y14 = diff(fy, itR);
y15 = diff(fy, itP);
y16 = diff(fy, itH);

foc = fopen('d20180204_resection.txt', 'w');
fprintf(foc, '%s\n%s\n%s\n%s\n%s\n%s\n\n', ccode(x11), ccode(x12), ccode(x13), ccode(x14), ccode(x15), ccode(x16));
fprintf(foc, '%s\n%s\n%s\n%s\n%s\n%s\n\n', ccode(y11), ccode(y12), ccode(y13), ccode(y14), ccode(y15), ccode(y16));
fprintf(foc, '%s\n%s\n', ccode(fx), ccode(fy));
fclose(foc);
open('d20180204_resection.txt');
%}

%{
% 后方交会
a11 = diff(fx, xs);
a12 = diff(fx, ys);
a13 = diff(fx, zs);
a14 = diff(fx, ro);
a15 = diff(fx, pt);
a16 = diff(fx, he);

a21 = diff(fy, xs);
a22 = diff(fy, ys);
a23 = diff(fy, zs);
a24 = diff(fy, ro);
a25 = diff(fy, pt);
a26 = diff(fy, he);

ca11 = ccode(a11);
ca12 = ccode(a12);
ca13 = ccode(a13);
ca14 = ccode(a14);
ca15 = ccode(a15);
ca16 = ccode(a16);
ca21 = ccode(a21);
ca22 = ccode(a22);
ca23 = ccode(a23);
ca24 = ccode(a24);
ca25 = ccode(a25);
ca26 = ccode(a26);

cfx = ccode(fx);
cfy = ccode(fy);

foc = fopen('d20160103_panoEquation.txt', 'w');
fprintf(foc, '%s\n%s\n%s\n%s\n%s\n%s\n\n', ca11, ca12, ca13, ca14, ca15, ca16);
fprintf(foc, '%s\n%s\n%s\n%s\n%s\n%s\n\n', ca21, ca22, ca23, ca24, ca25, ca26);
fprintf(foc, '%s\n%s\n', cfx, cfy);
fclose(foc);
open('d20160103_panoEquation.txt');
%}