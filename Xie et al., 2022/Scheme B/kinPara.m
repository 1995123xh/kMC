
%% get kinetic parameters

function [predic, KIE, KIE_13C] =kinPara(T)
R=8.314;
predic.CC=zeros(10,1);
CCpara=[1.071e20,-1.04,379942;  7.32325e24,-2.58,375564;  5.498e23/2,-1.91,371360;...
    1.353e33/4,-4.74,379810;  3.261e27,-3.19,373005; 2.718e26,-2.77,365689;
    5.694e30,-4.05,362698; 1.476e28,-3.11,359168; 1.6618e34,-5.08,357637;
    1.66423e38,-6.36,351299]; %parameters for CC reactions ordered by degree;...
%from 1-1, 1-2, 1-3, 1-4, 2-2, 2-3, 2-4, 3-3, 3-4, 4-4
predic.CC=CCpara(:,1).*T.^CCpara(:,2).*exp(-CCpara(:,3)/R/T);
%predic.CC(:)=6.64671E23*T^-1.99*exp(-375755/R/T);       %rates in s^-1 RMG-Py rate rules
predic.CA=1.58849e20*(T)^-0.76*exp(-442154/8.314/T); %
predic.CN=0; %4.68781e8*T^0.52*exp(-349516.61/R/T);
predic.CO=2.63224e18*T^-0.64*exp(-351864.4/R/T);
predic.CS=4.3545e13*T^0.56*exp(-305100.73/R/T);
predic.CAs=5.1379e25*T^-2.86*exp(-325726/R/T);    %for the secondary bond on the phenyl sidechains; reference: octylbenzene

DDH=zeros(5,4);
A_H=zeros(5,4);

DDH(1,1)=384; % 1 deg carbon
DDH(1,2)=328; % 2 deg carbon
DDH(1,3)=316; % 3 deg carbon

DDH(2,1)=359; % 1 deg C-A
DDH(2,2)=359+(DDH(1,2)-DDH(1,1));
DDH(2,3)=359+(DDH(1,3)-DDH(1,1));

DDH(3,1)=427; % 1 deg C-O
DDH(3,2)=427+(DDH(1,2)-DDH(1,1));
DDH(3,3)=427+(DDH(1,3)-DDH(1,1));

DDH(4,1)=377; % 1 deg C-S
DDH(4,2)=377+(DDH(1,2)-DDH(1,1));
DDH(4,3)=377+(DDH(1,3)-DDH(1,1));
DDH(:,4)=DDH(:,3);

DDH(5,1)=82.4; %secondary KIE Ni et al., 2011; Apply octane data to everything
DDH(5,2)=8.5;   % this value is suspicious (even though it comes from Ni et al. Might need recheck.
DDH(5,3)=DDH(5,2);
DDH(:,4)=DDH(:,3);

A_H(:)=1.20;   %pre-exponential factor
A_H(5,1)=1.06;
A_H(5,2)=1.02;
A_H(5,3)=A_H(5,2);
A_H(5,4)=A_H(5,2);

KIE = 1./(A_H.*exp(-DDH.*4.184/8.314/T));


% For carbon isotopes from table 2 in Tang et al., 2001
DDC=zeros(5,4);  %four for the 4 types of bonds, and 3 for carbon degrees; 5th row is secondary isotope effects
A_13C=zeros(5,4);

DDC(1,1)=48.853; % DEA of 1 deg carbon, in cal/mol
A_13C(1,1)=1.03;
DDC(1,2)=38.763; % 2 deg carbon
A_13C(1,2)=1.02;
DDC(1,3)=35.261; % 3 deg carbon
A_13C(1,3)=1.023;

DDC(2,1)=51.562; % 1 deg C-A
A_13C(2,1)=1.034;
DDC(2,2)=51.562+38.763-48.853; % 2 deg C-A 2021-08 update
A_13C(2,2)=1.018;
DDC(2,3)=51.562+35.261-48.853; % 3 deg C-A
A_13C(2,3)=1.018;

DDC(3,1)=51.539; % 1 deg C-O
A_13C(3,1)=1.034;
DDC(3,2)=51.539+38.763-48.853; % 2 deg C-O 2021-08 update
A_13C(3,2)=1.025;
DDC(3,3)=51.539+35.261-48.853; % 3 deg C-O
A_13C(3,3)=1.025;

DDC(4,1)=34.962; % 1 deg C-S
A_13C(4,1)=1.026; % 1 deg C-S
DDC(4,2)=34.962+38.763-48.853; % 2 deg C-S
A_13C(4,2)=1.022; % 2 deg C-S
DDC(4,3)=34.962+35.261-48.853; % 3 deg C-S
A_13C(4,3)=1.022; % 2 deg C-S

DDC(5,1)=4.69;   %equation 11b in Tang et al., 2005;
A_13C(5,1)=1.003;
DDC(5,2)=5.34;
A_13C(5,2)=1.00;
DDC(5,3)=DDC(5,2);
A_13C(5,3)=A_13C(5,2);

DDC(6,1)=47.442;   %for C-R bonds, Tang et al., 2005
A_13C(6,1)=1.039;
DDC(6,2)=28.314;
A_13C(6,2)=1.018;
DDC(6,3)=DDC(5,2);
A_13C(6,3)=A_13C(5,2);

DDC(:,4)=DDC(:,3);
A_13C(:,4)=A_13C(:,3);
KIE_13C=1./(A_13C.*exp(-DDC.*4.184./8.314./T));


end