function getMohrCoulomb ( E, nu, phi, cohesion, psi, sigt, dEps1, p_initial, iteration)
% *************************************************************************
% *                          MOHR-COULOMB MODEL                            *
% *************************************************************************
% Matlab code to simulate drained and undrained triaxial compression tests 
% with the original Mohr-Coulomb model. This code 
% was written in the simplest way for educational reasons and can only 
% simulate loading paths.
% 
% 
% INPUT PARMATERS
% E: Young's modulus
% nu: Poisson ratio
% phi: friction angle of soil in angle
% cohesion: cohesion of soil
% psi: dilatancy angle of soil
% sigt: maximum tension allowed in soil
% porosity: porosity of soil
% iteration: number of iterations
% dEps1: axial strain increment 
%
% OTHER PARAMETERS
% Nphi: (1 + sin phi) / (1 - sin phi)
% Npsi: (1 + sin psi) / (1 - sin psi)
%
% To run code, click on run or type the command line in the command window.
% You can change the parameters in the 'run' by clicking on the arrow,
% select a command line, righ click and edit. 

%
% example of command line:
% getMohrCoulomb(10^7, 0.2, 30, 1000, 0, 0, -0.0001, 0, 10)

%% ************************ ALLOCATE MEMORY *******************************    
       
    EpsV        = zeros(iteration,1);
    Eps1        = zeros(iteration,1);
    
    Sig         = zeros(6,1);
    Eps         = zeros(6,1); 
    dEps        = zeros(6,1); 
    De          = zeros(6,6);     
   
%%*************************************************************************
% *                         INITIALIZATION                                *
% *************************************************************************
    % Initial index
    a = 1;   
    
    % Initial state
    Sig(1,1) = p_initial;
    Sig(2,1) = p_initial;
    Sig(3,1) = p_initial;
    Nphi = (1 + sind(phi)) / (1 - sind(phi));
    Npsi = (1 + sind(psi)) / (1 - sind(psi));
    
    % get elastic moduli
    K = E / 3 / (1 - 2 * nu);
    G = E / 2 / (1 + nu);
    a1 = K + 4/3*G;
    a2 = K - 2/3*G;
     
    %build elastic stiffness matrix
    for x=1:6
        for y=1:6
            if x<=3
                if y<=3
                    if x==y
                        De(x,y)=a1; 
                    else
                        De(x,y)=a2;
                    end
                end
            end
            
            if x>3
                if y>3
                    if x==y
                        De(x,y)=2*G;
                    end
                end
            end
        end 
    end
 
%%*************************************************************************
% *                         LOADING                                       *
% *************************************************************************
    
    % Loading condition
    dEps(2) = dEps1;
    dEps(1) = - nu * dEps1;
    dEps(3) = - nu * dEps1;
 
%{
    Eps_1 = zeros(6,1); 
    Eps_2 = zeros(6,1);  
    
    % OLD:
    Eps_1(1) = 0.000378924;
    Eps_1(2) = -0.00052;
    Eps_1(3) = 0;
    
    Eps_2(1) = 0.000378954;
    Eps_2(2) = -0.00052002;
    Eps_2(3) = 0;
    
    dEps = Eps_2 - Eps_1;
    
    Sig(1, 1) = 0.0553322;
    Sig(2, 1) = -3463.94;
    Sig(3, 1) = -692.776;

    % NEW:
    Eps_1(1) = 4.84113e-06;
    Eps_1(2) = -0.00052;
    Eps_1(3) = 0;
    
    Eps_2(1) = 4.8415e-04;
    Eps_2(2) = -0.00052002;
    Eps_2(3) = 0;
    
    dEps = Eps_2 - Eps_1;
    
    Sig(1, 1) = -1390.65;
    Sig(2, 1) = -5764.33;
    Sig(3, 1) = -1431;    
%}    
    while a < iteration
        
        % Assume elastic guess
        SigTrial = Sig(:, a) + De * dEps;
        
        % Rotate to get pricipal stresses
        % Rotate x and y with z direction being the into the page
        % Working in principal stress to get stress
        sx = SigTrial(1);
        sy = SigTrial(2);
        sz = SigTrial(3);
        sxy = SigTrial(4);
        
        mohr_c = 0.5 * (sx + sy);
        mohr_r = 0.5 * sqrt((sx - sy)^2 + (2 * sxy)^2);
        
        s1 = mohr_c - mohr_r;
        s2 = mohr_c + mohr_r;

        if (s1 > sz)
            sigma1 = sz;
            sigma2 = s1;
            sigma3 = s2;
            mohr_flag = 2;
        elseif (s2 < sz)
            sigma1 = s1;
            sigma2 = s2;
            sigma3 = sz;
            mohr_flag = 3;
        else
            sigma1 = s1;
            sigma2 = sz;
            sigma3 = s2;
            mohr_flag = 1;
        end
%{
s = [s1; s2; sz]
%}
        disp(['Step ', num2str(a)]);
        sig = [sigma1; sigma2; sigma3];
 
        % Mohr-Coulomb failure criteria
        % shear
        yield_shear = sigma1 - sigma3 * Nphi + 2 * cohesion * sqrt(Nphi); 
        disp(['yield_shear: ', num2str(yield_shear)]);
        
        % tension
        sigt_max = cohesion / (tand(phi) + 0.0001);
        if sigt > sigt_max
            sigt = sigt_max;
        end
        yield_tension = sigt - sigma3;
        disp(['yield_tension: ', num2str(yield_tension)]);

        % Get hsig to know whether it is failing in shear or tension
        alphaP = sqrt(1 + Nphi^2) + Nphi;
        sigmaP = sigt * Nphi - 2 * cohesion * sqrt(Nphi);
        h = sigma3 - sigt + alphaP * (sigma1 - sigmaP);

        % Shear potential function gs (non-associated flow rule) - not used
        gs = sigma1 - sigma3 * Npsi;
    
        % Plastic part
        if (yield_tension < 0)
            if (h <= 0)
                % failure in shear and tension
                lambda_s = yield_shear / ( (a1 - a2 * Npsi) - (a2 - a1 * Npsi) * Nphi );
                sigma1 = sigma1 - lambda_s * (a1 - a2 * Npsi);
                sigma2 = sigma2 - lambda_s * a2 * (1 - Npsi);
                sigma3 = sigma3 - lambda_s * (- a1 * Npsi + a2);
                disp(['Step ', num2str(a), ': Yield in shear and tension']);
            else
                % failure only in tension
                lambda_t = yield_tension / a1;
                sigma1 = sigma1 + lambda_t * a2;
                sigma2 = sigma2 + lambda_t * a2;
                sigma3 = sigma3 + lambda_t * a1;
                sigt = 0;
                disp(['Step ', num2str(a), ': Yield only in tension']);
            end
        elseif (yield_shear < 0)
            % failure only in shear
            lambda_s = yield_shear / ( (a1 - a2 * Npsi) - (a2 - a1 * Npsi) * Nphi );
            sigma1 = sigma1 - lambda_s * (a1 - a2 * Npsi);
            sigma2 = sigma2 - lambda_s * a2 * (1 - Npsi);
            sigma3 = sigma3 - lambda_s * (- a1 * Npsi + a2);
            disp(['Step ', num2str(a), ': Yield only in shear']); 
        else
            disp(['Step ', num2str(a), ': No yield']);
        end
  
       %% ******************** UPDATE STATE ***********************
        % update index
        a = a + 1;
        
        % update stress in cartesian coordinate
        if (sigma1 == sigma3)
            cs2 = 1.0;
            si2 = 0.0;
        else
            cs2 = (sx - sy) / (s1 - s2);
            si2 = 2.0 * sxy / (s1 - s2);
        end
        
        if (mohr_flag == 1)
            dc2 = (sigma1 - sigma3) * cs2;
            dss = sigma1 + sigma3;
            Sig(1, a) = 0.5 * (dss + dc2);
            Sig(2, a) = 0.5 * (dss - dc2);
           	Sig(3, a) = sigma2;
            Sig(4, a) = 0.5 * (sigma1 - sigma3) * si2;
            Sig(5, a) = 0.0;
            Sig(6, a) = 0.0;
        elseif (mohr_flag == 2)
            dc2 = (sigma2 - sigma3) * cs2;
            dss = sigma2 + sigma3;
            Sig(1, a) = 0.5 * (dss + dc2);
            Sig(2, a) = 0.5 * (dss - dc2);
            Sig(3, a) = sigma1;
            Sig(4, a) = 0.5 * (sigma2 - sigma3) * si2;
            Sig(5, a) = 0.0;
            Sig(6, a) = 0.0;
        else
            dc2 = (sigma1 - sigma2) * cs2;
            dss = sigma1 + sigma2;
            Sig(1, a) = 0.5 * (dss + dc2);
            Sig(2, a) = 0.5 * (dss - dc2);
            Sig(3, a) = sigma3;
            Sig(4, a) = 0.5 * (sigma1 - sigma2) * si2;
            Sig(5, a) = 0.0;
            Sig(6, a) = 0.0;
        end
        
        % update strain variables
        % dEpsV = dEps(1) + dEps(2) + dEps(3);
        dEpsV = dEps1;
        EpsV(a, 1)   =  EpsV(a-1, 1) + dEpsV;      
        Eps1(a, 1)   =  Eps1(a-1, 1) + dEps1;
       
    end
    
%{   
Sig
Eps1
Sig(2,:)' ./ Eps1(:)
%}
Sig
%%*************************************************************************
% *                         PLOT RESULTS                                  *
% ************************************************************************* 
subplot(2,1,1) 
plot(-Eps1(:)*100, -Sig(2,:), '-d');
xlabel('Axial Strain [%]');
ylabel('Axial Stress [Pa]');

subplot(2,1,2)
plot(-Eps1(:)*100, -EpsV(:)*100, '-d');
xlabel('Axial Strain [%]');
ylabel('Volumetric Strain [%]');

end