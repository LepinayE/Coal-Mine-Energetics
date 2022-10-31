function [results] = solveAquifertemp2(K_d,ncyc,years,K_r, Q_flow, Tinj,...
                                        Taq , h, Ly, Lx, rho, C_p)

% solveAquifertemp2 - solves the Advection-Diffusion equation in
% the fracture and the heat equation in the surroudning rock.
% Model contains dispersion, diffusion, step function flow velocity 
% and constant input flux. 
%
% Inputs:
%   K_d - Effective thermal dispersivity of fracture ( 10^-5 to 10^-7)
%   ncyc - Total number of cycles 
%   years - Number of years system is run
%   K_r - Molecular heat diffusivity of rock (10^-7)
%   Q_flow - Volumetric flow rate (m^3/s)
%   Tinj - Temperature of injected fluid (Kelvin)
%   Taq - Initial aquifer temperature (°C)
%   h - Height of fracture (model is symmetrical and takes h/2)
%   Ly - Dimensional vertical length of aquifer
%   Lx - Dimensional horizontal length of aquifer
%   rho - Density of water kg/m^3
%   C_p - Specific heat capacity of fluid J/kgK
%
% Outputs:
%    results - structure with NON-DIM fields:
%                       .frac - T_f(X,t), fracture temperature
%                       .rock = T_r(X,Y,t), rock temperature
%                       .velocity = U, flow velocity
%                       .dt - time step
%                       .t_vec - vectorised time s.t t_vec(k) = k*dt
%                       .X_vec - vectorised horizontal length space
%                       .dX - horizontal step
%                       .Y_vec - vectorised vertical length space
%                       .dY - vertical step
%                       .Qhat - input heat flux
%
% Other m-files required: flowVelocity.m
%
% See also: CoalmineEnergy.m
addpath '/Users/lepinay/Desktop/Aquifer Matlab/FluxbcStepAquifer'
%--------------------------------------------------------------------------
% Non-dimentional 2D Cyclic Aquifer mode
% 
% Non-dimentional variables: 
% y = b/2 * Y - Vertical height
% t = b^2/(4*K) * \tau - Time from start of injection
% u = U0* U - Injection and extraction flow velocity
% x = ((b^2 * U0)/(4*K) )  * X - Horizontal Length
% Pe = (bU0)^2 /4 * K_r * K_d - Peclet number
% T_r = (T_inj - T_aq) \hat{T_r} + T_aq - Rock temperature
% T_f = (T_inj - T_aq) \hat{T_f} + T_aq - Fracture temperature
%
%------------------------
% Fracture model T_f(x,t) (drop hat)
% T_t + U T_X = (1/Pe) T_XX + df/dY(at Y=1) - Advection-Diffusion Eq
%
% For injection FDM forward: 
% T_f(i,k+1) = T_f(i,k ) - U(k)*(dt/dX) (T_f(i,k) -T_f(i-1,k)) ...
%               + (dt/dY)*(T_r(j+1,k,i) - T_r(j,k,i)) ...
%               + (dt/dX^2)*(1/Pe)*(T_f(i+1,k) -2 T_f(i,k) + T_f(i-1,k))
%
% For extraction FDM backward:
% T_f(i,k+1) = T_f(i,k ) + U(k)*(dt/dX) *(T_f(i+1,k)- T_f(i,k)) ...
%               + (dt/dY)*(T_r(j+1,k,i) - T_r(j,k,i)) ...
%               +(dt/dX^2)*(1/Pe)*(T_f(i+1,k) -2 T_f(i,k) + T_f(i-1,k))
%
%------------------------
% Rock Model T_r(y,t,x)
% T_t = T_YY 
% Using FDM:
% T_r(j, k+1,i) = T_r(j,k,i) + dt/(dY)^2 *(T_r(j+1,k,i) + ...
%                   T_r(i,j-1,k) -2 T_r(i,j,k) )
%
%------------------------
% Boundary conditions
%
% T_r(Y = 1,t,x) = T_f(x,t) - Instantaneous equilibrium between rock and
% fracture
% T_r(inf,t,X) = 0 - Far field rock temp held at aquifer initial temperature
% T_r(Y,t,inf) = 0
% 
% (T_f)_ XXX(inf,t) = 0 - No loss of heat and fluid at end of fracture
% 
% During injection: constant input flux
% Qhat = U * ( T_f + (T_aq/(T_inj - T_aq)) ) - (1/Pe)* ((T_f)_X) at X=0
% 
% During extraction: purely advective
% ((T_f)_X) = 0 at X = 0
%
%%------------------------
% Initial conditions
% T_r(y,0,x) = 0 - Initial (non-dim) temperature of aquifer
% T_f(x,0) = 0
%
%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 31/08/2022; Last revision: 6/09/2022
% Version: R2022a

%------------- BEGIN CODE -------------------------------------------------

    % Model Set up
    % Thermal Injection flux Q_thermal
    Q_thermal = rho*C_p * (Tinj -Taq) * (Q_flow /(h * Lx));

    % Speed
    U0 = Q_flow/(Lx * h);

    % Non-dim input flux
    Qhat = Q_thermal / ((Tinj-Taq)*U0*rho*C_p); 

    % Peclet number
    Pe = ((h^2)* (U0^2))/(4*K_d *K_r) ; 

    %-------------
    % Vertical height of domain
    LY = (2/h) *Ly;  % Non-dim length of cap rock cap
    NY = 151; % Number of elements
    Y_vec = linspace(0, LY, NY);  % Discretise Yspace starting from 0 ...
    %       and incrementing to Ly with length of vector Ny (LY/(NY -1))
    dY = Y_vec(2) - Y_vec(1); % Concequence of linspace 
    
    % Horizontal length of domain
    LX = ((4*K_r)/(U0*h^2)) *Lx; % Non-dim length of fracture
    NX = 151;  % Number of elements
    X_vec = linspace(0, LX, NX);  % Discretise Xspace starting from 0 ...
    %                       and incrementing to Lx with length of vector Nx
    dX = X_vec(2) -X_vec(1); 

    %-------------
    % Selecting correct time step to satisfy all conditions
    % d\tau < dX/(3) - Advection condition
    % d\tau < ((dY)^2)/2 - Courant ondition for rock
    % d\tau <((dX)^2)Pe/2 - Courant ondition for fracture

    minPe = ((h^2)* (U0^2))/(4*3*10^(-5) *K_r); % To keep dt constant ...
    %                       over experiments with different dispersivity
    
    if ( (dY^2)/ 4 <= dX/(3) ) && ( (dY^2)/ 4 <= ((dX)^2)* (minPe/ 4) )
    
        dt = ((dY^2)/ 4) ;
    
    elseif ( ((dX)^2)* (minPe/ 4) <= (dY^2)/ 4 ) && ...
            ( ((dX)^2)* (minPe/ 4) <= dX/(3) )
    
        dt  = ((dX)^2)* (minPe/ 4);
    else
        dt = dX/(3);
    
    end
    t_vec = 0:dt:(years+dt) ; % Discretise time such that t_vec(end) = years
    period = years/ncyc; % Period of oscillations

    %-------------
    % Flow velocity
    U = flowVelocity(t_vec,period); % Step function over each cycles

    %-------------
    % Satisfying ICs
    T_r = zeros (length(Y_vec), length(t_vec), length(X_vec)); % Rock ...
    %                                           Temp T_r(Y,tX)
    T_f = zeros (length(X_vec), length(t_vec)); % Fracture T_f(X,t)


    % FDM Matrix for rock T_r
    beta = dt/ ((dY)^2);

    A0 = (1 - 2* beta) * ones(length(Y_vec),1); 
    Ap1 = beta *  ones(length(Y_vec)-1,1);
    Am1 = beta *  ones(length(Y_vec)-1,1);
    
    A = diag(A0) + diag(Ap1, 1) + diag (Am1, -1); % 2D array
    
    A(1,:) = 0 ; % To impose BCs
    A(end,:) = 0; 
    
    A = repmat(A, 1, 1, length(X_vec)); % Creates 3D array by repeatedly...
    %                                     stacking A in the 3rd directions
    
    %-------------
    % Courant condition check 
    alpha2 = (1/Pe) * (dt/(dX)^2);
    
    if alpha2 >= 0.5
       msg = 'Courant condition not satisfied for fracture';
         error(msg)
    end

    %-------------
   % FDM Solver
    % Ghost variable at (i = 1, j =1)
    for k = 1 : length(t_vec) -1 
    
        % BC at end of domain T_xxx = 0
        T_f(end, k ) = 3 * T_f(end-1, k)- 3* T_f(end-2, k) ...
                        + T_f(end-3, k) ;
        %Far field BC for rock temp
        T_r(end,k,:) = 0;
        T_r(:, k , end) = 0;

        %-------------
        % Loop over Injection and Extraction cases to build Matrix
        if U(k) >= 0 % Injection
           
            % Injection BC for T_f 
            T_f(1,k) = ( (Pe* dX)/ (( U(k)* Pe* dX ) + 1) )* ( Qhat - ( (U(k)* Taq)/ (Tinj - Taq) ) +( T_f(2,k)/ (Pe* dX)) );
            

            %-------------
            % FDM Matrix input for fracture T_f
            alpha1 = U(k) * (dt/(dX) );

            C0 = (1 - alpha1 - 2*alpha2 ) * ones(length(X_vec),1) ;
            Cm1 = (alpha1 + alpha2) * ones(length(X_vec)-1,1);
            Cp1 = alpha2* ones(length(X_vec)-1,1);

            %-------------
            % Build FDM Matrix for fracture T_f
            C = diag(C0) + diag(Cp1,1)+ diag(Cm1,-1);
            C(end,:) = 0; % To impose far field BC
            C(1,:) = 0 ; % To impose BC during injection

        elseif U(k) < 0 % Extraction

            % Extraction BC for T_f
            T_f(1,k) = T_f(2,k);

            %-------------
            % FDM Matrix input for fracture T_f
            alpha1 = U(k) * (dt/(dX) );

            C0 = (1 + alpha1 - 2*alpha2 ) * ones(length(X_vec),1);
            Cp1 = (- alpha1 + alpha2) * ones(length(X_vec)-1,1);
            Cm1 = alpha2* ones(length(X_vec)-1,1);
        
            %-------------
            % Build FDM Matrix for fracture T_f
            C = diag(C0) + diag(Cp1,1)+ diag(Cm1,-1);
            C(end,:) = 0; % To impose far field BC
            C(1,:) = 0 ; % To impose BC during extraction
           
        end % Ends if loop to switch between injection and extraction

      
        %-------------
        % Solve fracture equation
         
        M1 = squeeze(T_r(2,k,:) - T_r(1,k,:)); % Heat loss to rock
        
        T_fnext = C*T_f(:,k) + (dt/dY)*(M1); % Find Temp at next time step
        
        T_f(2:length(X_vec) - 1, k+1) = T_fnext(2:length(X_vec) - 1);


        %-------------

        % Solve rock equation 
        T_rnext = pagemtimes(A,T_r( :, k,:)) ;
        
        T_r( 2:length(Y_vec) - 1, k+1,2:length(X_vec) - 1) = ...
            T_rnext( 2:length(Y_vec) - 1, 1, 2:length(X_vec) - 1);
        
       % relate fracture temp to rock temp at boundary y = b/2
        T_r(1,k+1,:) = T_f(:,k+1);

    end

    T_r = shiftdim(T_r, 2); % Shift dimensions such that T_r(x,y,t)

    % Removing ghost variables
    T_f = T_f(2:end-1,:); 

    T_r = T_r(2:end-1,:,:);
    %-------------
    % Checking errors
    
     fractsum = sum(T_f,2);
     if fractsum(end) > 0.0001
         msg = 'Fracture Temperature changes past domain size';
         error(msg)
     end
   
     rockYsum = T_r(2,:,:);
     rockYsum = squeeze(rockYsum);
     rockYsum = sum(rockYsum, 2);
     if rockYsum(end) > 0.001
         msg = 'Rock Temperature changes past Y-domain size';
         error(msg)
     end

     rockXsum = T_r(:,2,:);
     rockXsum = squeeze(rockXsum);
     rockXsum = sum(rockXsum,2);
     if rockXsum(end) > 0.001
         msg = 'Rock Temperature changes past X-domain size';
         error(msg)
     end

     

    %-------------
    % Outputs

    results.frac = T_f;
    results.rock = T_r;

    % Additional info
    results.velocity = U;
    results.dt = dt;
    results.t_vec = t_vec;
    results.X_vec = X_vec;
    results.dX = dX;
    results.Y_vec = Y_vec;
    results.dY = dY;
    results.Qhat = Qhat;
end

%------------- END OF CODE ------------------------------------------------

