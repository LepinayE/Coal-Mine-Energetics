function [deltaP, Re] = pumpinAquifer2D(Lp,rho,D,mu,Q_flow,epsilon,Lx,...
                                        kperm,h)

%pumpinAquifer - Pressure change experienced when pumping down a well and
%into a horizontal fracture. All in 3D cylindrical coordinate
%
%Syntax:  [deltaP, Re] = pumpinAquifer(L,rho,D,mu,Q,epsilon,r_frac,...
%                                        hydraulicK,h)
% Inputs:
%    Lp - Length of injection pipe (m)
%    rho - Density of fluid (kg/m^3)
%    D - Injection pipe diameter (meters)
%    mu - viscosity of water (kg/ms)
%    Q_flow - Volumetric flow rate (m^3/s)
%    epsilon -  Roughness height coefficient of the injection well, based
%               on material used (meters but will usually be given in mm)
%    Lx - Linear extent of fracture (m)
%    kperm - Permeability of fracture (m^2)
%    h - Thickness of fracture (m)
%
% Outputs:
%    deltaP - Sum of pressure change in injection well deltaPwell and down
%    fracture deltaPfracture. Units Pa or kg/ms^2.
%               deltaPwell - Pressure difference between pressure at top of
%               well and pressure at bottom. Derived from the Darcy Weisbach
%               equation and Moody chart. Assume injection well is 3D cylinder
%
%               deltaPfracture - Pressure diference between injection well
%               and end of fracture. Derived from Darcy's Law. Assume
%               Fracture is 2D (Height y and length x). 
%    
%    Re = Reynold's number.
%%%% ASSUME SPEED FROM PIPE IS SAME DOWN FRACTURE
%
% Other m-files required: colebrook.m
%
% See also: CoalmineEnergy.m

%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 04/10/2022; Last revision: 03/11/2022
% Version: R2022a

addpath '/Users/lepinay/Desktop/Aquifer Matlab/Energetical Efficiency'
%------------- BEGIN CODE -------------------------------------------------
    % Pressure Change down Injection Well
    
    % Reynold's number 
    Re = (4 *rho *Q_flow ) /(mu * pi *D) ;
    
    % Finding the Darcy friction coefficient f_d
    
    if Re > 4000 % Assumes Flow is turbulent
        
        relrough = epsilon / D; % Relative roughness coefficient
        
        f_d = colebrook(Re,relrough); % Colebrook-White Equation
    
    elseif Re <= 4000 % Assumes Flow is Laminar
        
        f_d = 64/Re ;
    
    end

    deltaPwell= (8 * Lp * f_d * rho * (Q_flow^2))/((pi^2) * D^5);

    % Velocity

    U0 = Q_flow * (1/h) * (1/Lx);
    %-------------
    % Pressure Change down Fracture
    
    % Parameters

    deltaPfracture =  U0 * Lx * mu / kperm; % Assume linear flow

    %-------------
    % Outputs
    deltaP = deltaPwell + deltaPfracture ;
end
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
