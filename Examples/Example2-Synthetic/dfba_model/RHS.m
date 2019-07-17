function [lb,ub] = RHS( t,y,INFO )

% This subroutine updates the upper and lower bounds for the fluxes in the
% exID arrays in main. The output should be two matrices, lb and ub. The lb matrix
% contains the lower bounds for exID{i} in the ith row in the same order as
% exID. The same is true for the upper bounds in the ub matrix.
% Infinity can be used for unconstrained variables, however, it should be 
% fixed for all time. 

% Y1 = Biomass EColi (gDW/L)
% Y2 = Carbon (mmol/L)
% Y3 = Nitrogen (mmol/L)
% Y4 = Oxygen (mmol/L)
% Y5 = Lipids (mmol/L)
% Y6 = Ethanol (mmol/L)
% Y7 = Oxidation product (mmol/L)
% Y8 = Penalty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
if isfield(INFO,'param')
    vmaxc = INFO.param(1);
    Kc = INFO.param(2);
    Kie = INFO.param(3);
    vmaxn = INFO.param(4);
    Kn = INFO.param(5);
    vmaxo = INFO.param(6);
    Ko = INFO.param(7);
    vATPm = INFO.param(8);
else
    vmaxc = 1.5;
    Kc = 0.05;
    Kie = 15;
    vmaxn = 0.25;
    Kn = 0.5;
    vmaxo = 2;
    Ko = 1.2;
    vATPm = 0.18;
end

% Carbon
lb(1,1) = 0;
if (y(2)<0)
    ub(1,1) = 0;
else
    ub(1,1) = vmaxc*y(2)/(Kc + y(2))*1/(1+y(6)/Kie);
end

% Nitrogen
lb(1,2) = 0;
if (y(3)<0)
    ub(1,2) = 0;
else
    ub(1,2) = vmaxn*y(3)/(Kn + y(3));
end

% Oxygen
lb(1,3) = 0;
if (y(4)<0)
    ub(1,3) = 0;
else
    ub(1,3) = vmaxo*y(4)/(Ko + y(4));
end

% ATP
lb(1,4) = vATPm;
ub(1,4) = vATPm;
end