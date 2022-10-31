function [results] = solveCoalmineRobin(K_d,ncyc,years,K_r, U0, h, Ly,...
                                        Lx)

% solveCoalmineRobin - solves the Advection-Diffusion equation in
% the fracture and the heat equation in the surroudning rock with Robin
% boundary condition during injection. Model contains dispersion, 
% diffusion, step function flow velocity and constant input flux.
% FDM explicit scheme requires courant conditions for stability.
%
%
% Inputs:
%   K_d - Effective thermal dispersivity of fracture ( 10^-5 to 10^-7)
%   ncyc - Total number of cycles 
%   years - Number of years system is run
%   K_r - Molecular heat diffusivity of rock (10^-7)
%   U0 - Fluid speed in fracture (10^-5) 
%   h - Height of fracture (model is symmetrical and takes b/2)
%   Ly - Dimensional vertical length of aquifer
%   Lx - Dimensional horizontal length of aquifer
%    
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
%                      
%
% Other m-files required: flowVelocity.m
% Subfunctions: none
% Other files: HeatconsAquiferRobin.m, simpleAquifertempGraphRobin.m
%
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
% EXPLICIT 
% Fracture model T_f(x,t) (drop hat)
% T_t + U T_X = (1/Pe) T_XX + df/dY(at Y=1) - Advection-Diffusion Eq
%
% For injection FDM forward time, backward convection, CS: 
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
% T_r(inf,t,X)_YY = 0 - Far field rock temp held at aquifer initial temperature
% T_r(Y,t,inf) = 0
% 
% (T_f)_ X(inf,t) = 0 - No loss of heat and fluid at end of fracture
% 
% During injection: constant input flux (backward difference)
% U * ( T_f + -1 ) = (1/Pe)* ((T_f)_X) at X=0
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
% Date: 14/10/2022; Last revision: 31/10/2022
% Version: R2022a

addpath '/Users/lepinay/Desktop/Aquifer Matlab'/FluxbcStepAquifer
%------------- BEGIN CODE -------------------------------------------------

    % Model Set up

    Pe = ((h^2)* (U0^2))/(4*K_d *K_r) ; % Peclet number

    %-------------
    % Vertical height of domain
    LY = (2/h) *Ly;  % Non-dim length of cap rock cap
    NY = 251; % Number of elements
    Y_vec = linspace(0, LY, NY);  % Discretise Yspace starting from 0 ...
    %       and incrementing to Ly with length of vector Ny (LY/(NY -1))
    dY = Y_vec(2) - Y_vec(1); % Concequence of linspace 
    
    % Horizontal length of domain
    LX = ((4*K_r)/(U0*h^2)) *Lx; % Non-dim length of fracture
    NX = 251;  % Number of elements
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
    T_r = zeros(length(Y_vec), length(t_vec), length(X_vec)); % Rock ...
    %                                           Temp T_r(Y,t,X)
    T_f = zeros(length(X_vec), length(t_vec)); % Fracture T_f(X,t)


    % FDM Matrix for rock T_r A(y,t)
    beta = dt/ ((dY)^2);

    A0 = (1 - 2* beta) * ones(length(Y_vec),1); 
    Ap1 = beta *  ones(length(Y_vec)-1,1);
    Am1 = beta *  ones(length(Y_vec)-1,1);
    
    A = diag(A0) + diag(Ap1, 1) + diag (Am1, -1); % 2D array (Y by Y)

  
    
    % Impose  far field BCs
    A(end,end-1) = 0; 
    A(end,end) = 1;
    
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

   for k = 1 : length(t_vec) -1 
    
 
        %-------------
        % Loop over Injection and Extraction cases to build Matrix
        if U(k) >= 0 % Injection
           
            % Injection BC for T_f 
            T_f(1,k) = ( (Pe* dX)/((U(k)* Pe* dX) + 1) )...
                            * ( U(k)  + (T_f(2,k)/ (Pe* dX)) );

            %-------------
            % FDM Matrix input for fracture T_f
            alpha1 = U(k) * (dt/(dX) );

            C0 = (1 - alpha1 - 2*alpha2 ) * ones(length(X_vec),1) ;
            Cm1 = (alpha1 + alpha2) * ones(length(X_vec)-1,1);
            Cp1 = alpha2* ones(length(X_vec)-1,1);

            %-------------
            % Build FDM Matrix for fracture T_f
            C = diag(C0) + diag(Cp1,1)+ diag(Cm1,-1);

            %-------------
            % Impose far field BC
            C(end,end) = 1 - alpha1 - alpha2;
            C(end, end-1) = alpha1 + alpha2;

            

        elseif U(k) < 0 % Extraction

            % Extraction BC for T_f
            T_f(1,k) = T_f(2,k); % SHOUDL MAYBE BE THE OTHER WAY AROUND

            %-------------
            % FDM Matrix input for fracture T_f
            alpha1 = U(k) * (dt/(dX) );

            C0 = (1 + alpha1 - 2*alpha2 ) * ones(length(X_vec),1);
            Cp1 = (- alpha1 + alpha2) * ones(length(X_vec)-1,1);
            Cm1 = alpha2* ones(length(X_vec)-1,1);
        
            %-------------
            % Build FDM Matrix for fracture T_f
            C = diag(C0) + diag(Cp1,1)+ diag(Cm1,-1);

            %-------------
            % Impose far field BC
            C(end,end) = 1 - alpha1; 
            C(end,end-1) =  alpha1;
            
           
        end % Ends if loop to switch between injection and extraction

      
        %-------------
        % Solve fracture equation
         
        M1 = squeeze(T_r(2,k,:) - T_r(1,k,:)); % Heat loss to rock
        
        T_fnext = C*T_f(:,k) + (dt/dY)*(M1); % Find Temp at next time step
        
        T_f(1:length(X_vec) , k+1) = T_fnext(1:length(X_vec));


        %-------------

        % Solve rock equation 
        T_rnext = pagemtimes(A,T_r( :, k,:)) ;
        
        T_r( 1:length(Y_vec) , k+1,1:length(X_vec)) = ...
            T_rnext( 1:length(Y_vec), 1, 1:length(X_vec));
        
       % Relate fracture temp to rock temp at boundary y = b/2
        T_r(1,k+1,:) = T_f(:,k+1);
          
    
    end % End of FDM solver

    T_r = shiftdim(T_r, 2); % Shift dimensions such that T_r(x,y,t)


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
   
end

%------------- END OF CODE ------------------------------------------------

