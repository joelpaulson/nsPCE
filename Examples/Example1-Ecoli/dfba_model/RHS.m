%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez and Kai Höffner                                %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., Höffner, K. and Barton, P. I. (2014).                      %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. BMC Bioinformatics, 15:409                                    % 
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lb,ub] = RHS( t,y,INFO )

% This subroutine updates the upper and lower bounds for the fluxes in the
% exID arrays in main. The output should be two matrices, lb and ub. The lb matrix
% contains the lower bounds for exID{i} in the ith row in the same order as
% exID. The same is true for the upper bounds in the ub matrix.
% Infinity can be used for unconstrained variables, however, it should be 
% fixed for all time. 

% Y1 = Volume (L)
% Y2 = Biomass EColi (gDW/L)
% Y3 = Glucose (mmol/L)
% Y4 = Xylose (mmol/L)
% Y5 = O2 (mmol/L)
% Y6 = Ethanol (mmol/L)
% Y7 = Penalty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Extract parameters [see "Generalized derivatives of dynamic systems with
%%% a linear program embedded", equations 25-27 or "Dynamic Flux Balance
%%% Modeling of MicrobialCo-Cultures for Efficient Batch Fermentation
%%% ofGlucose and Xylose Mixtures"]
p = INFO.param;
vgmax = p(1); % maximum uptake rate of glucose
Kg = p(2); % saturation constant for glucose
vzmax = p(3); % maximum uptake rate of xylose
Kz = p(4); % saturation constant for xylose
Kig = p(5); % inhibition constant for glucose
vo = p(6); % uptake rate for oxygen (do not use MM because we assume extracelluar oxygen conc controlled)

% EColi
% Glucose
% if (y(3)<0)
%     lb(1,1) = 0;
% else
%     lb(1,1) = -(10.5*y(3)/(0.0027 + y(3)))*1/(1+y(6)/20);
% end
% ub(1,1) = 0;
if (y(3)<0)
    lb(1,1) = 0;
else
    lb(1,1) = -(vgmax*y(3)/(Kg + y(3)))*1/(1+y(6)/20);
end
ub(1,1) = 0;

% Xylose
% if (y(4)<0)
%     lb(1,2) = 0;
% else
%     lb(1,2) = -(6*y(4)/(0.0165 + y(4)))*1/(1+y(6)/20)*1/(1+y(3)/0.005);
% end
% ub(1,2) = 0;
if (y(4)<0)
    lb(1,2) = 0;
else
    lb(1,2) = -(vzmax*y(4)/(Kz + y(4)))*1/(1+y(6)/20)*1/(1+y(3)/Kig);
end
ub(1,2) = 0;

% O2
% lb(1,3) = -(15*y(5)/(0.024 + y(5)));
% ub(1,3) = 0;
lb(1,3) = -vo;
ub(1,3) = 0;

% Ethanol
lb(1,4) = 0;
ub(1,4) = Inf;
end