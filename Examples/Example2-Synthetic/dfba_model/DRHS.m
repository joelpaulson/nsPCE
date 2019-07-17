function dy = DRHS(t, y, INFO)

% Y1 = Biomass EColi (gDW/L)
% Y2 = Carbon (mmol/L)
% Y3 = Nitrogen (mmol/L)
% Y4 = Oxygen (mmol/L)
% Y5 = Lipids (mmol/L)
% Y6 = Ethanol (mmol/L)
% Y7 = Oxidation product (mmol/L)
% Y8 = Penalty
X = y(1);
C = y(2);
N = y(3);
O = y(4);
L = y(5);
E = y(6);
COX = y(7);

%% Update bounds and solve for fluxes
[flux,penalty] = solveModel(t,y,INFO);
%%

% set fluxes to zero if penalty is > 0
if penalty(1) > 1e-6
    flux = zeros(1,7);
end

% growth
mu = flux(1,1);
% lipid
vlip = flux(1,2);
% fermentation
vferm = flux(1,3);
% carbon
vc = flux(1,4);
% nitrogen
vn = flux(1,5);
% oxygen
vo = flux(1,6);
% oxidation products
vcox = flux(1,7);

% get uncertain coefficients from INFO
if isfield(INFO,'param')
    SOXOX = INFO.param(19);
    SOXFerm = INFO.param(20);
else
    SOXOX = 1.0;
    SOXFerm = 2.0;
end


%% Dynamics
dy = zeros(8,1);
dy(1) = mu*X;
dy(2) = -vc*X;
dy(3) = -vn*X;
dy(4) = -vo*X;
dy(5) = vlip*X;
dy(6) = vferm*X;
dy(7) = (SOXOX*vcox + SOXFerm*vferm)*X;
dy(8) = penalty(1);
end