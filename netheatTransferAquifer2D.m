function [Q_heatsystem,Q_heatavg, time] = netheatTransferAquifer2D(rho,...
                                            C_p,U0,K_d,ncyc,years, K_r, ...
                                            Tinj, Taq ,h, Ly, Lx)


% heatTransferAquifer2D - Calculates the Thermal Energy added and retrieved
%                           to the system in J/s for a 2D channel and 
%                           isothermal fracture.
%
% Syntax:  [Q_heat] = heatTransferAquifer(rho,c_p,Q)
%
% Inputs:
%    rho - density of fluid (kg/m^3)
%    Cp - heat capacity of fluid (J/kg°C)
%    Q_flow - Volumetric flow rate (m^3/s)
%    U0 - Velocity in fracture (m/s) assumed uniform and linear
%    K_d - Effective thermal dispersivity of fracture ( 10^-5 to 10^-7)
%    ncyc - Total number of cycles 
%    years - Number of years system is run
%    K_r - Molecular heat diffusivity of rock (10^-7) 
%    Tinj - Temperature of injected fluid
%    Taq - Initial aquifer temperature
%    h - Height of fracture (model is symmetrical and takes h/2)
%    Ly - dimensional vertical length of aquifer
%    Lx - dimensional horizontal length of aquifer
%
% Outputs:
%    Q_heatsystem - net heat transferred into and out of system J/s 
%    Q_heatavg - net averaged heat transferred into and out of system
%    time - t_vec from system solver
%
%
%
% Other m-files required: solveCoalmineRobin.m
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 31/10/2022; Last revision: 3/11/2022
% Version: R2022a

%------------- BEGIN CODE -------------------------------------------------
    

    % Solve Aquifer temp
    
    results = solveCoalmineRobin(K_d,ncyc,years,K_r, U0, h, Ly, Lx);
    
    %-------------
    % Finding extracted Temperature

    onofflogic = ischange(results.velocity); %Logical 1 at abrupt changes

    onoffIndex = find(onofflogic); % Time index of changes in velocity
    onoffIndex(2:2:end) = onoffIndex(2:2:end) -1; % End of extraction is ...
    % one index before start of injection
   

    extTf = zeros (size(results.t_vec));
    for n = 1 : ncyc % Fracture temp at x=0
        
        extTf(onoffIndex((2*n)-1):onoffIndex(2*n)) = ...
                    results.frac(1,onoffIndex((2*n)-1):onoffIndex(2*n));
       
        % Averages Extraction temp
        extTfavg(n) = mean (extTf(onoffIndex((2*n)-1):onoffIndex(2*n))); 
    end
    
    
    %-------------
    % Change in Temperature between injected temperature and extracted
    % temperature 
    
    deltaTsystem = extTf - 1; % Since T_inj non dim is 1 
    
    deltaTavg = extTfavg -1;

    % Dimensionalise

    deltaTsystem = (Tinj - Taq) *deltaTsystem; 
    
    deltaTavg = (Tinj - Taq)* deltaTavg; 

    

    %-------------
    % Outcome
    Q_heatsystem = rho *C_p * U0 * deltaTsystem; %Units J/sm^2
    Q_heatavg = rho *C_p * U0 * deltaTavg;
    time = results.t_vec;

end
%------------- END OF CODE ------------------------------------------------
%