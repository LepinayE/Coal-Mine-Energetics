% NetThermalPower.m - File compares the power demand of pumping into
% aquifer based on permeability and with the power output of heat storage.
% Graph produced show permeability 
%
% 
%
% Inputs: Common aquifer parameters with focus on for loops over
%    kperm - Permeability of aquifer (m^2)
%    Qflow - Volumetric flow rate (m^3/s)
%
% Outputs:
%    SECTION 1:Demand and Output of each system
%              - Pumping cost for different permeability
%              - Typical output of a specific heat storage system
%              - Section 1 graph
%    SECTION 2: Net Thermal Power for different flow rate and deltaT
%               figure(2)
%
% Subfunctions: none
% MAT-files required: pumpinAquifer2D.m, netheatTransferAquifer2D.m
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 06/12/2022; Last revision: 
% Version: R2022a

%------------- BEGIN CODE -------------------------------------------------


% NEED TO SET iQf and kperm FOR EACH SECTION

%% Establishing Parameters

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
h = [1, 2];

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
Taq = [10, 10]; 

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
years = [1, 10, 20];

% Temperature of injected fluid (°C)
Tinj = [20, 60, 100]; 

%% Parameter Selection

% Pipe length
iLp = 1;

% Viscosity
imu = 1;

% Thickness of fracture
ih = 1;

% Effective dispersivity
iKd = 2;

% Number of cycles
incyc = 1;

% Years
iyear = 1;

% Permeability 
kperm = logspace(-13,-12,20);

%% SECTION 1: Demand and Output of each system
 
%-------------
% Pumping cost for different permeability

% Flow rate 
iQf = 2; %0.01 m^3/s

j = 1;

for k = kperm ;


% Work done by system at each second
    [deltaP, Re] = pumpinAquifer2D(Lp(iLp),rho,D,mu(imu),Q_flow(iQf),epsilon,...
                                    Lx, k,h(ih));
    WorkdonePersecond = deltaP * Q_flow(iQf);% *ones(1,length(time)); %Joules/second
    
    %WdPressure = cumsum(WorkdonePersecond); %Joules

    Wdpressvperm(j) = 2*WorkdonePersecond;

    j = j+1 ;
end

%-------------
% Power output of Heat storage system

% Initial Aquifer Temp
iTaq = 2;

% Injection Temp
iTinj = 3;

% Velocity in fracture
U0 = Q_flow(iQf) * (1/h(ih)) * (1/(Lx));
    
% Rate of heat retrieved from system J/m^2s
[Q_heatext,Q_heatsystem,Q_heatavg, time] = netheatTransferAquifer2D(rho,...
                                            C_p,U0,K_d(iKd),ncyc(incyc),...
                                            years(iyear), K_r,...
                                            Tinj(iTinj), Taq(iTaq) ,h(ih),...
                                            Ly, Lx);

netQ_heatavg = -1* Lx * h(ih) * Q_heatavg;

%% SECTION 1: Graph 

figure(1)

% Pumping Graph vs Permeability
semilogx(kperm,Wdpressvperm, '-', 'LineWidth',3, "Color",[18 62 116]/255 )
    
grid on


%Changing y labels 

yt = yticks;
yt = yt/1000000;
yticklabels(yt);
ylim([0 2000000])


xlim([10^(-13) 1])

% To do to get split
%break_axis('axis', 'x', 'position', 0.5)
%xlim([10^(-13) 10^(-10)])
%xticks(0)


hold on 

% Thermal Power output of system and Heat pump power usage
kperm2 = logspace(-13, -10,30);
Qh = netQ_heatavg.*ones(1,length(kperm2));

semilogx(kperm2,Qh, '-', 'LineWidth',3, "Color",[104 172 229]/255 )
hold on 

semilogx(kperm2,Qh/3, '--', 'LineWidth',3, "Color",[242 32 35]/255 )



ylabel (" Workdone (MW)", 'FontSize', 16,'Interpreter','latex')
xlabel("Permeability",'FontSize', 16,'Interpreter','latex')
legend("Workdone pumping fluid through aquifer","Typical power produced ..." + ...
    "using heat storage","Workdone using a heat pump for same heat produced")
legend("location","east")

%% SECTION 2: Net Thermal Power for different flow rate and deltaT

% Flow rate 0.001
iQf = 1;

% Permeability 
kperm = logspace(-13,-10,20);

j = 1;

for k = kperm ;

    % Work done by system at each second
    [deltaP, Re] = pumpinAquifer2D(Lp(iLp),rho,D,mu(imu),Q_flow(iQf),epsilon,...
                                    Lx, k,h(ih));
    WorkdonePersecond = deltaP * Q_flow(iQf);% *ones(1,length(time)); 

    WdprespermSlow(j) = WorkdonePersecond; %Joules/second

    j = j+1 ;
end


j = 1;

% Flow rate 0.01
iQf = 2;

for k = kperm ;

    % Work done by system at each second 
    [deltaP, Re] = pumpinAquifer2D(Lp(iLp),rho,D,mu(imu),Q_flow(iQf),epsilon,...
                                    Lx, k,h(ih));
    WorkdonePersecond = deltaP * Q_flow(iQf);% *ones(1,length(time));

    WdprespermFast(j) = 2*WorkdonePersecond; %Joules/second

    j = j+1 ;
end

 


figure(2)

%c = [[0 0.4470 0.7410]; [0.6350 0.0780 0.1840]; [0 1 0.7410]; [0.8500 0.3250 0.0980];];
c = uint8([[18 62 116]; [242 32 35]; [104 172 229]; [237 103 62];]);

for i = [1,2] % for loop over temperature

    % Initial Aquifer Temp
    iTaq = i;
    
    % Injection Temp
    iTinj = i;

    % Loop over Velocity in fracture
    for iQf = [1,2];
    
        U0 = Q_flow(iQf) * (1/h(ih)) * (1/(Lx));
    
        % Rate of heat retrieved from system J/m^2s
        [Q_heatext,Q_heatsystem,Q_heatavg, time] = netheatTransferAquifer2D(...
                                                rho,C_p,U0,K_d(iKd),...
                                                ncyc(incyc),years(iyear), K_r, ...
                                                Tinj(iTinj), Taq(iTaq) ,h(ih),...
                                                Ly, Lx);
        
        % Net Q_heat retrieved J/s
        
        %netQ_heatsystem =-1* Lx * h(ih) * Q_heatsystem; % 1- extTf at each time step
        netQ_heatavg = -1* Lx * h(ih) * Q_heatavg; % 1- average of extracted temp 
        
        if iQf == 1
    
            Wdslow = netQ_heatavg - WdprespermSlow(:);
        elseif iQf == 2
    
            Wdfast = netQ_heatavg - WdprespermFast(:);
        end
    
    end

% Workdone Overall

    if i == 1
        dT = Tinj(i) - Taq(i);
        
    elseif i == 2
        dT = Tinj(i) - Taq(i);

    end

    txt = ['Slow, Temperature change = ',num2str(dT) , ' $^{\circ}$C'];
    txt2 = ['Fast, Temperature change = ',num2str(dT) , ' $^{\circ}$C'];

    semilogx(kperm,Wdslow, '-', "MarkerSize",10,'LineWidth',3, "Color",c(i,:),'DisplayName', txt)
    
    hold on 
    
    semilogx(kperm,Wdfast,'-', "MarkerSize",10,'LineWidth',3, "Color",c(i+2,:), 'DisplayName', txt2 )
    
    hold on

end

xlabel("Permeability ($m^2$)",'FontSize', 16,'Interpreter','latex')
ylabel (" Net Power output of system (MW)", 'FontSize', 16)

%Change yticks to be in MW
yt = yticks;
yt = yt/1000000;
yticklabels(yt);

legend('location', 'southeast','Interpreter','latex', 'FontSize', 14)
grid on 



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

