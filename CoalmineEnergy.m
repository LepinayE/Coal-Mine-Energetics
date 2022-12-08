% MASTER FILE
%
% CoalmineEnergy - Comparing the energy needed to run the system vs energy
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
%
%
%--------------------------------------------------------------------------
% The thermal energy added and retrieved to the system in J/s for a 2D 
% channel and isothermal fracture.
%
%[Q_heatext,Q_heatsystem,Q_heatavg, time] = netheatTransferAquifer2D(rho,...
%                                   C_p,U0,K_d,ncyc, years, K_r, Tinj,...
%                                   Taq, h, Ly, Lx)
%
%
% Outputs:
%    Q_heatsystem - net heat transferred into and out of system J/s 
%    Q_heatavg - net averaged heat transferred into and out of system
%    time - t_vec from system solver
%    Q_heatext - Heat retrieved during extraction J/s
%
%
%
%-------------
% ASSUMES velocity in fracture is the same as the velocity in the injection
% pipe
%
%-------------
% Other m-files required: colebrook.m, flowVelocity.m, solveCoalmineRobin.m
% MAT-files required: pumpinAquifer2D.m, netheatTransferAquifer2D.m

%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 10/10/2022; Last revision: 3/11/2022
% Version: R2022a

%------------- BEGIN CODE -------------------------------------------------

%% Parameter selection

%-------------
% Injection Pipe
% Length of injection pipe (m): Coal min, i = 3 is Eavor-Lite, i=4 is Eden
Lp = [100, 1000, 2400, 4500];

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
h = [1, 2 ];

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
Lx = 1000;

% Initial aquifer temperature (°C)
Taq = [20, 75, 180]; 

%-------------
% System Properties
% Volumetric flow rate (m^3/s)
Q_flow = [10^(-3),10^(-2), 0.00973,0.008,10^(-7), 0.0583]; 
% i = 2 is Eavor-lite
% i = 3 is Eden
%1 tonne per hr, Netherland case study

% Total number of cycles
ncyc = [1, 10, 20, 50];

% Number of years system operates for
years = [0.5, 1, 10, 20];

% Temperature of injected fluid (°C)
Tinj = [25, 40, 60, 80]; 

%% Pumping cost vs kperm changing
%-------------
% Choice of parameter

% Pipe length
iLp = 1;

% Viscosity
imu = 1;

% Thickness of fracture
ih = 1;

% Effective dispersivity
iKd = 2;

% Flow rate
iQf =2;


j = 1;

for kperm = logspace(-13,-10,20);

    k = kperm;
% Work done by system at each second
    [deltaP, Re] = pumpinAquifer2D(Lp(iLp),rho,D,mu(imu),Q_flow(iQf),epsilon,...
                                    Lx, kperm,h(ih));
    WorkdonePersecond = deltaP * Q_flow(iQf);% *ones(1,length(time)); %Joules/second
    
    %WdPressure = cumsum(WorkdonePersecond); %Joules

    Wdpresperm(j) = WorkdonePersecond;

    j = j+1 ;
end

figure(10)

kperm = logspace(-13,-10,20);
semilogx(kperm,Wdpresperm,'--.', "MarkerSize",10,'LineWidth',3, "Color",[0 0.4470 0.7410])
grid on

xlabel("Permeability (m^2)",'FontSize', 16)
ylabel (" Workdone Pumping into System (J/s)", 'FontSize', 16)

 

%% Coal Mine scenario

%-------------
% Choice of parameter

% Pipe length
iLp = 1;

% Viscosity
imu = 1;

% Thickness of fracture
ih = 1;

% Permeability 
ikperm = 3;

% Effective dispersivity
iKd = 2;

% Initial Aquifer Temp
iTaq = 1;

% Flow rate
iQf =1;

% Number of cycles
incyc = 1;

% Years
iyear = 1;

% Injection Temp
iTinj = 2;

%-------------
% Velocity in fracture
U0 = Q_flow(iQf) * (1/h(ih)) * (1/(Lx));

% Rate of heat retrieved from system J/m^2s
[Q_heatext,Q_heatsystem,Q_heatavg, time] = netheatTransferAquifer2D(rho,C_p,U0,K_d(iKd),...
                                        ncyc(incyc),years(iyear), K_r, ...
                                        Tinj(iTinj), Taq(iTaq) ,h(ih),...
                                        Ly, Lx);

% Net Q_heat retrieved J/s

netQ_heatsystem =-1* Lx * h(ih) * Q_heatsystem; % 1- extTf at each time step
netQ_heatavg = -1* Lx * h(ih) * Q_heatavg; % 1- average of extracted temp 

% Difference between injected temp and aquifer temp
injtime = length(time)/ncyc(incyc);
netInjQ_heat = Q_flow(iQf) * rho* C_p *(Tinj(iTinj) - Taq(iTaq)); 
netInjQ_heat = netInjQ_heat *ones(1,injtime);

%-------------
% Work done by system at each second
[deltaP, Re] = pumpinAquifer2D(Lp(iLp),rho,D,mu(imu),Q_flow(iQf),epsilon,...
                                Lx,kperm (ikperm),h(ih));
WorkdonePersecond = deltaP * Q_flow(iQf) *ones(1,length(time)); %Joules/second

WdPressure = cumsum(WorkdonePersecond); %Joules




%%
%-------------
% Percentage of Energy needed to pump into system compared to Energy
% retrieved
CostHeatvPressavg =  WdPressure(end) / netQ_heatavg; % Average heat retrieved
display(CostHeatvPressavg)

CostHeatvPresssys =  WdPressure ./ netQ_heatsystem;


% Heat injected vs Aquifer temp
CostInjvPres = WdPressure(end) / netInjQ_heat(1); 
display(CostInjvPres)


%-------------
%Graph

figure(5)
plot(time, netQ_heatsystem, '.b', 'DisplayName', 'Thermal Energy injected and extracted from system (Watts)')

legend show
xlabel('Time (yrs)')
ylabel ('Thermal Power (kW)')

%Changing y labels 
yt = yticks;
yt = yt/1000;
yticklabels(yt);




figure(3)

plot(WdPressure,  'DisplayName','Workdone pumping (J)')
legend show

figure(1)

plot(netQ_heatsystem, 'DisplayName', 'Net Heat injected and retrieved from system')

hold on 

plot(netInjQ_heat, 'DisplayName', 'Net Heat Transferred during Injection')
grid on 
legend show

figure(2)
plot(CostHeatvPresssys)

%% Net Thermal Power (Heat flux - pressure work per kperm)

Q_flow = [0.01,0.001];

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