% MASTER FILE
%
% CoalmineEnergy - Comparing the energy needed to run the system vs enerygy
%                   output.
% 
% Inputs:
%    Lp - Length of injection pipe (m)
%    rho - Density of fluid (kg/m^3)
%    D - Injection pipe diameter (meters)
%    mu - viscosity of water (kg/ms)
%    Q_flow - Volumetric flow rate (m^3/s)
%    epsilon -  Roughness height coefficient of the injection well, based
%               on material used (meters but will usually be given in mm)
%    kperm - Permeability of fracture (m^2)
%    h - Thickness of fracture (m)
%    Cp - heat capacity of fluid (J/kg°C)    
%    K_d - Effective thermal dispersivity of fracture ( 10^-5 to 10^-7)
%    ncyc - Total number of cycles 
%    years - Number of years system is run
%    K_r - Molecular heat diffusivity of rock (10^-7) 
%    Tinj - Temperature of injected fluid (kelvin K)
%    Taq - Initial aquifer temperature K
%    Ly - dimensional vertical length of aquifer (m)
%    Lx - dimensional horizontal length of aquiferm (m)
%
%
%--------------------------------------------------------------------------
% Pressure change in whole system whcih needs to be worked agaisnt when
% pumping in.
%
% [deltaP, Re] = pumpinAquifer2D(Lp,rho,D,mu,Q_flow,epsilon,Lx,...
%                                       kperm,h)
%
% Outputs:
%    deltaP - Sum of pressure change in injection well deltaPwell and down
%    fracture deltaPfracture. Units Pa or kg/ms^2.
%               deltaPwell - Pressure difference between pressure at top of
%               well adn pressure at bottom. Derived from the Darcy Weisbach
%               equation and Moody chart. Assume injection well is 3D cylinder
%               deltaPfracture - Pressure diference between injection well
%               and end of fracture. Derived from Darcy's Law. Assume
%               Fracture is 2D.
%    
%    Re = Reynold's number.
%
%--------------------------------------------------------------------------
% The thermal energy added and retrieved to the system in J/s for a 2D 
% channel and isothermal fracture.
%
% [netQ_heat] = netheatTransferAquifer2D(rho,Cp,Q_flow,K_d,ncyc,...
%                                                years, K_r, Tinj, Taq ,...
%                                                h, Ly, Lx)
%
%
% Outputs:
%    Q_heat - net heat transferred into system J/s 
%
%-------------
% ASSUMES velocity in fracture is the same as the velocity in the injection
% pipe
%
%-------------
% Other m-files required: colebrook.m, flowVelocity.m, solveCoalmineRobin.m
% MAT-files required: pumpinAquifer2D.m, 

%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 10/10/2022; Last revision: 
% Version: R2022a

%------------- BEGIN CODE -------------------------------------------------

%% Parameter selection

%-------------
% Injection Pipe
% Length of injection pipe (m): Coal min, i = 4 is Eavor
Lp = [100, 500, 1000, 4500];

% Injection pipe diameter (m)
D = 0.25;

% Roughness height coefficient of the injection well, based on material ...
% used (meters but will usually be given in mm)
epsilon = 1.5*10^(-6); %PVC

%-------------
% Fluid properties 
% Density of fluid (kg/m^3)
rho = 1000;

% Viscosity of fluid (kg/ms)
mu = [10^(-3), 0.28*10^(-3)]; % @ 20degrees, 100 degrees

% Heat capacity of fluid (J/kg°C) 
C_p = 4200;

%-------------
% Aquifer properties
% Thickness of fracture (m)
h = [1 2];

% Permeability of fracture (m^2)
kperm = [10^(-7),10^(-11), 10^(-13)]; % Highly fractured rock, Netherland...
%                                       case study, Oil reservoir

% Effective thermal dispersivity of fracture (m^2/s)
K_d = [10^(-7), 10^(-6), 10^(-5)];

% Molecular heat diffusivity of rock (m^2/s)
K_r = 10^(-7);

% Dimensional vertical length of aquifer (m)
Ly = 6;

% Dimensional horizontal length of aquifer (m)
Lx = 3000;

% Initial aquifer temperature (°C)
Taq = [20, 180]; 

%-------------
% System Properties
% Volumetric flow rate (m^3/s)
Q_flow = [2.778*10^(-4),10^(-7), 0.0583]; %1 tonne per hr, Netherland case study

% Total number of cycles
ncyc = [1, 10, 20, 50];

% Number of years system operates for
years = [0.2, 1, 10, 20];

% Temperature of injected fluid (°C)
Tinj = [40, 60, 80]; 


%% Coal Mine scenario

%-------------
% Choice of parameter

% Pipe length
iLp = 4;

% Viscosity
imu = 1;

% Thickness of fracture
ih = 1;

% Permeability 
ikperm = 1;

% Effective dispersivity
iKd = 1;

% Initial Aquifer Temp
iTaq = 1;

% Flow rate
iQf =1;

% Number of cycles
incyc = 1;

% Years
iyear = 1;

% Injection Temp
iTinj = 1;


%-------------
% Velocity in fracture
U0 = 4*Q_flow(iQf) * (1/pi) * (1/(D^2));

% Rate of heat transferred to system J/m^2s
[Q_heat, time] = netheatTransferAquifer2D(rho,C_p,U0,K_d(iKd),...
                                        ncyc(incyc),years(iyear), K_r, ...
                                        Tinj(iTinj), Taq(iTaq) ,h(ih),...
                                        Ly, Lx);
% Net Q_heat J/s

netQ_heat = Lx * h(ih) * Q_heat; % Unsure about this wetted area

%-------------
% Work done by system at each second
[deltaP, Re] = pumpinAquifer2D(Lp(iLp),rho,D,mu(imu),Q_flow(iQf),epsilon,...
                                Lx,kperm (ikperm),h(ih));
WorkdonePersecond = deltaP * Q_flow(iQf) *ones(1,length(time)); %Joules/second

WdPressure = cumsum(WorkdonePersecond); %Joules



%------------- END OF CODE ------------------------------------------------
%
%
%
%
%
%
%
%
%
%
%
%