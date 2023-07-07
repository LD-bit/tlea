% AUTHOR:   DUBREIL Léa
% DATE:     09/02/2023 (created) 06/04/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Extended Kalman Filter with Doppler observations (main
%           contribution)
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
sigma_vk = 1e-2;        %[m/s] Velocity standard deviation
sigma_pk = 1e0;         %[m] Position standard deviation
                                        
% OBSERVATION MODEL
% Uncertainties in measures
sigma_nk = 1;               %[m] range
sigma_tk = 1e-5;            %[rad] theta, elevation 
sigma_lk = 1e-5;            %[rad] phi, azimuth
% Uncertainties on the doppler are expressed via the radial velocity.
% Reason: reduce singularity in Sk (det ~0) for inversion
c = 3e8;                    %[m/s] speed of light, celerity
sigma_fk1 = 10*c/2.1e9;     %[m/s] radial velocity for 1st doppler (S-band)
sigma_fk2 = 10*c/2.1e9;     %[m/s] radial velocity for 2nd doppler (S-band)
sigma_fk3 = 10*c/70e6;      %[m/s] radial velocity for 3rd doppler (UHF)
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
%% ------------------------------ MAIN ZONE -------------------------------
disp('[INFO] Extended Kalman Filter with Doppler observations (main contrib)')
h = waitbar(0,'Please wait...');

% EXTENDED KALMAN FILTER
tic

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
    
    for k = 1:N
        
        % INITIALISATION: for each sampling time
        % noise n_k ~ N(0,Rk)
        n_k = chol(Rk)'*randn(6,1); 
        % noise v_k ~ N(0,Qk) 
        v_k = chol(Qk)'*randn(6,1);

        for g = 0:(mdT/T)-1
            % TRUE STATES
            if k == 1
                xt = x0; % initialisation
            else
                xt = xt + fk(xt)*T + v_k; % eq. 1.45, true state
            end
        end
        % MEASUREMENT
        z = Doppler_hk(xt,xt) + n_k; % eq. 1.46, measurement z[k]

        % PREDICTION 
        if k == 1
            % initialisation of x[n|n] 
            % NB: it is != of init state (true state):
            x = x0 + randn(6,1).*[50 50 50 5 5 5]';
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
            
            
            % norm of the diff for each state between true and filtered
            MeasErr(k,loop) = norm(x - xt)^2;
            MeasErrp(k,loop) = norm(x(1:3) - xt(1:3))^2; 
            MeasErrv(k,loop) = norm(x(4:6) - xt(4:6))^2;      
        end
    end % end for KF   
end % end for of MC runs
close(h)
% CONVERGENCE
MeasMSE = sqrt((1/nMC)*sum(MeasErr,2));
MeasMSEp = sqrt((1/nMC)*sum(MeasErrp,2)); 
MeasMSEv = sqrt((1/nMC)*sum(MeasErrv,2)); 
toc

%% PLOTS
figure
% MSE position
subplot(2,1,1),plot(t,MeasMSEp,'-r')
grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('RMSE [m]','Interpreter','latex','fontsize',16)
title('Root Mean Square Error')
legend('RMSE: position')
% MSE velocity
subplot(2,1,2),plot(t,MeasMSEv,'b')
grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('RMSE [m/s]','Interpreter','latex','fontsize',16)
title('Mean Square Error','Interpreter','latex','fontsize',16)
legend('RMSE: velocity','Interpreter','latex','fontsize',16)
sgtitle({['Doppler integrated EKF method with sampling time $T_s =$ ' num2str(T) ' [s] and measurements every' ' ' num2str(mdT) ' [s]']},'Interpreter','latex','fontsize',20)
%% SAVE PLOTS FOR LaTeX
%print -deps DopplerEKF_T_4