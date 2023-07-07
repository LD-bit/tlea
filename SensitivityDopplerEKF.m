% AUTHOR:   DUBREIL Léa
% DATE:     09/02/2023 (created) 31/03/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Sensitivity analysis of the Doppler Extended Kalman Filter
%           regarding the uncertainties on the evolution model
% REF:      Samy Labsir, 2020, Méthodes statistiques fondées sur les 
%           groupes de Lie pour le suivi d'un amas de débris spatiaux.

clear all; close all; clc;
global s_1 s_2 s_3
%% ---------------------------- VARIABLES ZONE ----------------------------
% DISCRETISATION OF TIME
T       = 1e-1;         %[s] sampling period
mdT     = 1;            %[s] measuring period
N       = 60*8;         %[s] max duration of measurements
step    = N/T;          %[-] step
t       = 0:T:N*T-T;    %[s] time vector
nMC     = 1e2;          %[-] Monte Carlo number (keep it min 100)

% EVOLUTION MODEL
% Uncertainties in model
sigma_ck_vk = 0.05;    %[m/s] Velocity standard deviation
sigma_ck_pk = 1;      %[m] Position standard deviation
% Sensitivity analysis:
max_bound = 3; min_bound = 1e-5; n_point = 20;
e1 = linspace(min_bound,10,n_point)/T;
e2 = linspace(min_bound,0.5,n_point)/T;
                                        
% OBSERVATION MODEL
% Uncertainties in measures
sigma_nk = 1;               %[m] range
sigma_tk = 1e-5;            %[rad] theta, elevation 
sigma_lk = 1e-5;            %[rad] phi, azimuth
sigma_fk1 = (10)/(2.1e9);   %[-] normalised doppler 1
sigma_fk2 = (10)/(2.1e9);   %[-] normalised doppler 2
sigma_fk3 = (10)/(70e6);    %[-] normalised doppler 3
sigma_vvk = 1e-2;           %[m/s] measured velocities

% Initial State Vector
% the initialisation of this vector is for ops-sat state vec, on Feb 20,
% 2023 (TLE age was 0.4 days). Position in [m] and velocity in [m/s]
x0 = [-3089022.5178; %[m] x
      -3952006.3867; %[m] y
      4687365.0410;  %[m] z
      2105.6276;     %[m/s] vx
      4864.7572;     %[m/s] vy
      5469.7445];    %[m/s] vz

% RECEIVERS:
% First receiver O1: IZN-1
h_1 = 2429.662;                     %[m] altitude
lat_1 = 28.29959791;                %[°] North latitude
lon_1 = -16.51061160+360;           %[°] West longitude = -East + 360°
s_1 = spherical_to_cartesianECEF(lat_1,lon_1,h_1); %[m] Position vector 
% Second receiver O2: ESOC-1
h_2 = 125;                          %[m] altitude
lat_2 = 49.871111;                  %[°] North latitude
lon_2 = 8.622778;                   %[°] East longitude
s_2 = spherical_to_cartesianECEF(lat_2,lon_2,h_2); %[m] Position vector 
% Third receiver O3: TU Graz
h_3 = 389;                          %[m] altitude
lat_3 = 47.0678581;                 %[°] North latitude
lon_3 = 15.4500702;                 %[°] East longitude
s_3 = spherical_to_cartesianECEF(lat_3,lon_3,h_3); %[m] Position vector 

% PARAMETERS SAVED
xt = zeros(6,1); % true state of x
z = zeros(7,1); % measurements z
x = zeros(6,1); % state vector
MeasErr = zeros(N,nMC); % MSE total
MeasErrp = zeros(N,nMC); % MSE position
MeasErrv = zeros(N,nMC); % MSE velocity
RMSE_sigma = zeros(n_point,n_point); % Sensitivity study
%% ------------------------------ MAIN ZONE -------------------------------
disp('[INFO] Sensitivity analysis of the Doppler Extended Kalman Filter')
h = waitbar(0,'Please wait...');

% EXTENDED KALMAN FILTER
tic
j = 0;
for sigma_pk = e1    
    j = j+1
    i = 0;
    for sigma_vk = e2
        i = i+1;
        for loop = 1:nMC
        
            waitbar(loop/nMC,h);
        
            % INITIALISATION
            % covariance matrix for Doppler integration
            Rk = diag([sigma_nk^2; ... % range1
                        sigma_tk^2; ... % az1
                        sigma_lk^2; ... % el1
                        sigma_fk1^2; ... % fd1
                        %sigma_tk^2; ... % az2
                        %sigma_lk^2; ... % el2
                        sigma_fk2^2; ... % fd2
                        %sigma_tk^2; ... % az3
                        %sigma_lk^2; ... % el3
                        sigma_fk3^2;]); % fd3
            Qk = diag([sigma_pk^2*ones(3,1);sigma_vk^2*ones(3,1)]); % covariance matrix of evaluations
            P = [100*eye(3) zeros(3);...
                  zeros(3) eye(3)]; % covariance P[n|n]
            Ck = diag([(sigma_ck_pk)^2*ones(3,1);(sigma_ck_vk)^2*ones(3,1)]);
            for k = 1:N
                
                % INITIALISATION: for each sampling time
                % noise n_k ~ N(0,Rk)
                n_k = chol(Rk)'*randn(6,1); 
                % noise v_k ~ N(0,Qk) 
                v_k = chol(Qk)'*randn(6,1);
                c_k = chol(Ck)'*randn(6,1);
                for g = 0:(mdT/T)-1
                    % TRUE STATES
                    if k == 1
                        xt = x0; % initialisation
                    else
                        xt = xt + fk(xt)*T + c_k; % eq. 1.45, true state
                    end
                end
                % MEASUREMENT
                z = Doppler_hk(xt,xt) + n_k; % eq. 1.46, measurement z[k]
        
                % PREDICTION 
                if k == 1
                    % initialisation of x[n|n] 
                    % NB: it is != of init state (true state):
                    x = x0;% + randn(6,1).*[50 50 50 1 1 1]';
                    MeasErr(k,loop) = norm(x - xt)^2;
                    MeasErrp(k,loop) = norm(x(1:3) - xt(1:3))^2;
                    MeasErrv(k,loop) = norm(x(4:6) - xt(4:6))^2; 
                else
                    for g = 0:(mdT/T)-1
                        P = Jfk(x,T)*P*Jfk(x,T)' + Qk; % eq. 1.51, covariance estimation P[k|k-1]
                        x = x + fk(x)*T; % eq. 1.50 state estimation x[k|k-1]
                    end
                               
                    % CORRECTION
                    Sk = Doppler_Jhk(x,xt)*P*Doppler_Jhk(x,xt)' + Rk; % eq. 1.56 Update matrix
                    Kn = P*Doppler_Jhk(x,xt)'/Sk; % eq. 1.55, Kalman gain K[n]
                    x = x + Kn*(z - Doppler_hk(x,xt)); % eq. 1.53, state update x[k|k]
                    P = (eye(6) - Kn*Doppler_Jhk(x,xt))*P; % eq. 1.54, covariance update P[k|k]
                    
                    % ERRORS
                    % norm of the diff for each state between true and filtered
                    MeasErr(k,loop) = norm(x - xt)^2;
                    MeasErrp(k,loop) = norm(x(1:3) - xt(1:3))^2; 
                    MeasErrv(k,loop) = norm(x(4:6) - xt(4:6))^2;      
                end
            end % end for KF   
        end % end for of MC runs
        % CONVERGENCE
        % In order to study the convergence of the KF, the trace of Pn has to 
        % tend towards the value of the MSE for each sampling time
        MeasMSE = sqrt((1/nMC)*sum(MeasErr,2));
        MeasMSEp = sqrt((1/nMC)*sum(MeasErrp,2)); 
        MeasMSEv = sqrt((1/nMC)*sum(MeasErrv,2)); 
        RMSE_sigma(i,j) = MeasMSE(end); % sensitivity analysis
    end % end loop on sigma_vk
end % end loop on sigma_pk
toc
close(h)

%% PLOTS

figure
surf(e1*T,e2*T,10*log10(RMSE_sigma))
grid on
xlabel('Realisation of $\sigma_{pk}$ [m]','Interpreter','latex','fontsize',16)
ylabel('Realisation of $\sigma_{vk}$ [m/s]','Interpreter','latex','fontsize',16)
zlabel('$\sqrt{MSE}$ [dB]','Interpreter','latex','fontsize',16)
title('Sensitivity of uncertainties on evolution model for the Doppler EKF','Interpreter','latex','fontsize',16)
%% SAVE PLOTS FOR LaTeX
%print -deps DopplerEKF_T_4