clc
clear all
% 3D lift calculations
data2D = [
        -5.0000 -0.2446 0.0140;
        -4.0000 -0.1465 0.0091;
        -3.0000 -0.0401 0.0073;
        -2.0000 0.0658 0.0059;
        -1.0000 0.1717 0.0049;
        0.0000 0.2737 0.0043;
        1.0000 0.4058 0.0045;
        2.0000 0.5143 0.0050;
        3.0000 0.6167 0.0057;
        4.0000 0.7194 0.0066;
        5.0000 0.8201 0.0078;
        6.0000 0.9193 0.0092;
        7.0000 1.0129 0.0112;
        8.0000 1.1027 0.0134;
        9.0000 1.1844 0.0165;
        10.0000 1.2533 0.0201;
        11.0000 1.2865 0.0252;
        12.0000 1.2763 0.0332;
        13.0000 1.2329 0.0475;
        14.0000 1.1635 0.0720;
        15.0000 1.0951 0.1052;
        ];
    
alpha = data2D(:,1); % Angle of attack
Clift = data2D(:,2); % 2D coefficient of lift
CDrag = data2D(:,3); % 2D coefficient of drag

e = 0.9; % Span efficiency factor
AR = 16.5; % Aspect Ratio

%linear region
alpha_Upper = 5;
alpha_Lower = 4;

%Find index of alpha values
x1=find(data2D(:,1)==alpha_Lower);
x2=find(data2D(:,1)==alpha_Upper);

%initial estimate of alpha
a_0 = (Clift(x2)-Clift(x1))/(alpha(x2)-alpha(x1));

%estimate alpha for finite wing
a = (a_0)/(1 + ((57.3 * a_0)/(pi * e * AR)));

%find zero lift angle of attack
sign = find(Clift>0,1);
Alpha0 = (Clift(sign) - Clift(sign - 1))/(alpha(sign) - alpha(sign - 1)) * (alpha(sign));

%use alpha to calculate cl
CL_3D = a .* ( alpha - Alpha0);
 
% 3D coefficient of drag
D_i = (CL_3D).^2 ./ (pi*e*AR);

%calculate total drag
D = CDrag + D_i;

%find min drag alpha value
MinDragAlpha = alpha(find(D == min(D)));

%define variables
%oswald's efficiency factor
e0 = 1.78 * (1 - 0.045 .* AR.^(0.68)) - 0.64;
k = 1 / (pi * e0 * AR);

%For CFD
dataCFD=[
    -5,-0.32438,0.044251;
    -4,-0.21503,0.033783;
    -3,-0.10081,0.028627;
    -2,0.010503,0.025864;
    -1,0.12155,0.024643;
    0,0.24163,0.025099;
    1,0.34336,0.025635;
    2,0.45256,0.02766;
    3,0.56037,0.030677;
    4,0.66625,0.034855;
    5,0.76942,0.040403;
    6,0.86923,0.04759;
    7,0.96386,0.057108;
    8,1.0441,0.070132;
    9,1.0743,0.090921;
    10,1.0807,0.11193;
    11,1.0379,0.13254;
    12,1.034,0.15645;
    ];
alphaCFD = dataCFD(:,1); % Angle of attack
CliftCFD = dataCFD(:,2); % 3D coefficient of lift
CDragCFD = dataCFD(:,3); % 3D coefficient of drag

%plot and compare the results
figure(1)
plot(alpha, CL_3D);
hold on
plot(alphaCFD, CliftCFD);
hold on
plot(alpha, Clift);
hold off

%mark each plot
title('Compare Results of Calculated and CFD Data');
legend('3D Calculated Data', 'CFD Data', '2D Calculated Data');