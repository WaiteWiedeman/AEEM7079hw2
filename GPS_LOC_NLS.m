clc; clear; close all; 

data = load("gps_meas1.mat"); % load data
rho = data.rho; % measured pseudorange
sigm = 5; % standard deviation 5 meters
maxiter = 7;

% SV position coordinates [m]
SV_5 = [15.764733 -1.592675 21.244655]*1e6;
SV_13 = [6.057534 -17.186958 19.396689]*1e6;
SV_18 = [4.436748 -25.771174 -1.546041]*1e6;
SV_22 = [-9.701586 -19.687467 15.359118]*1e6;
SV_26 = [23.617496 -11.899369 1.492340]*1e6;
SV_27 = [14.540070 -12.201965 18.352632]*1e6;

% Weight
W = 1/(sigm^2)*eye(length(rho));

% Initial Conditions
xc = [0 0 0 0]';
xe = zeros(maxiter+1,4); xe(1,:)=xc';

% Nonlinear Least Squares
for i = 1:maxiter
    rho1e = sqrt((SV_5(1)-xc(1))^2 + (SV_5(2)-xc(2))^2 + (SV_5(3)-xc(3))^2) + xc(4);
    rho2e = sqrt((SV_13(1)-xc(1))^2 + (SV_13(2)-xc(2))^2 + (SV_13(3)-xc(3))^2) + xc(4);
    rho3e = sqrt((SV_18(1)-xc(1))^2 + (SV_18(2)-xc(2))^2 + (SV_18(3)-xc(3))^2) + xc(4);
    rho4e = sqrt((SV_22(1)-xc(1))^2 + (SV_22(2)-xc(2))^2 + (SV_22(3)-xc(3))^2) + xc(4);
    rho5e = sqrt((SV_26(1)-xc(1))^2 + (SV_26(2)-xc(2))^2 + (SV_26(3)-xc(3))^2) + xc(4);
    rho6e = sqrt((SV_27(1)-xc(1))^2 + (SV_27(2)-xc(2))^2 + (SV_27(3)-xc(3))^2) + xc(4);

    H1 = [ -[(SV_5(1)-xc(1)) (SV_5(2)-xc(2)) (SV_5(3)-xc(3))]/(rho1e-xc(4)) 1];
    H2 = [ -[(SV_13(1)-xc(1)) (SV_13(2)-xc(2)) (SV_13(3)-xc(3))]/(rho2e-xc(4)) 1];
    H3 = [ -[(SV_18(1)-xc(1)) (SV_18(2)-xc(2)) (SV_18(3)-xc(3))]/(rho3e-xc(4)) 1];
    H4 = [ -[(SV_22(1)-xc(1)) (SV_22(2)-xc(2)) (SV_22(3)-xc(3))]/(rho4e-xc(4)) 1];
    H5 = [ -[(SV_26(1)-xc(1)) (SV_26(2)-xc(2)) (SV_26(3)-xc(3))]/(rho5e-xc(4)) 1];
    H6 = [ -[(SV_27(1)-xc(1)) (SV_27(2)-xc(2)) (SV_27(3)-xc(3))]/(rho6e-xc(4)) 1];

    H = [H1; H2; H3; H4; H5; H6];

    dy = [rho(1)-rho1e;rho(2)-rho2e;rho(3)-rho3e;rho(4)-rho4e;rho(5)-rho5e;rho(6)-rho6e];
    dx = inv(H'*inv(sigm^2)*H)*H'*inv(sigm^2)*dy;

    xc = xc + dx;
    xe(i+1,:) = xc';

    sig3 = diag(inv(H'*inv(sigm^2)*H))'.^(0.5)*3;

end
