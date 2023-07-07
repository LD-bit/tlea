% AUTHOR:   DUBREIL Léa
% DATE:     16/02/2023 (created) 23/03/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Sensitivity analysis of the Classical Extended Kalman Filter
%           regarding the uncertainties on the evolution model
% REF:      Samy Labsir, 2020, Méthodes statistiques fondées sur les 
%           groupes de Lie pour le suivi d'un amas de débris spatiaux.

clear all; close all; clc;

%% ---------------------------- VARIABLES ZONE ----------------------------
% DISCRETISATION OF TIME
T       = 1e-4;       %[s] sampling period: measure every second
mdT     = 1e-4;       %[s] measuring period
N       = 60*8;       %[s] max duration of measurements
step    = N/T;        %[-] step
t       = 0:T:N*T-T;  %[s] time vector
nMC     = 2e2;        %[-] Monte Carlo number (keep it min 100)

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
sigma_nk = 1e-1;    %[m] standard deviation range
sigma_tk = 1e-6;    %[rad] theta 
sigma_lk = 1e-6;    %[rad] phi 

% Initial State Vector
% the initialisation of this vector is for ops-sat state vec, on Feb 20,
% 2023 (TLE age was 0.4 days). Position in [m] and velocity in [m/s]
x0 = [-3089022.5178; %[m] x
      -3952006.3867; %[m] y
      4687365.0410;  %[m] z
      2105.6276;     %[m/s] vx
      4864.7572;     %[m/s] vy
      5469.7445];    %[m/s] vz

% PARAMETERS SAVED
xt = zeros(6,1); % true state of x
z = zeros(3,1); % measurements z
x = zeros(6,1); % state vector
MeasErr = zeros(N,nMC); % measurement error in total
MeasErrp = zeros(N,nMC); % measurement error in position
MeasErrv = zeros(N,nMC); % measurement error in velocity
TrPn = zeros(1,N); % Trace of Pn
RMSE_sigma = zeros(n_point,n_point); % Sensitivity study
%% ------------------------------ MAIN ZONE -------------------------------
disp('[INFO] Sensitivity analysis of the Classical Extended Kalman Filter')
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
            Rk = diag([sigma_nk^2;sigma_tk^2;sigma_lk^2]); % covariance matrix of measurements
            Qk = diag([(sigma_pk*T)^2*ones(3,1);(sigma_vk*T)^2*ones(3,1)]); % covariance matrix of evaluations
            P = [100*eye(3) zeros(3);...
                  zeros(3) eye(3)]; % covariance P[n|n]
            Ck = diag([(sigma_ck_pk)^2*ones(3,1);(sigma_ck_vk)^2*ones(3,1)]);
            for k = 1:N
                
                % INITIALISATION: for each sampling time
                % noise n_k ~ N(0,Rk)
                n_k = chol(Rk)'*randn(3,1); 
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
                z = hk(xt) + n_k; % eq. 1.46, measurement z[k]
        
                % PREDICTION 
                    if k == 1
                        % initialisation of x[n|n] 
                        % NB: it is != of init state (true state):
                        x = x0 ;%+ randn(6,1).*[50 50 50 1 1 1]';
                        MeasErr(k,loop) = norm(x - xt)^2;
                        MeasErrp(k,loop) = norm(x(1:3) - xt(1:3))^2;
                        MeasErrv(k,loop) = norm(x(4:6) - xt(4:6))^2; 
                    else
                        for g = 0:(mdT/T)-1
                            P = Jfk(x,T)*P*Jfk(x,T)' + Qk; % eq. 1.51, covariance estimation P[k|k-1]
                            x = x + fk(x)*T; % eq. 1.50 state estimation x[k|k-1]
                        end
                   
                    % CORRECTION
                    Sk = Jhk(x)*P*Jhk(x)' + Rk; % eq. 1.56 Update matrix
                    Kn = P*Jhk(x)'/Sk; % eq. 1.55, Kalman gain K[n]
                    x = x + Kn*(z - hk(x)); % eq. 1.53, state update x[k|k]
                    P = (eye(6) - Kn*Jhk(x))*P; % eq. 1.54, covariance update P[k|k]
                    
                    % ERRORS
                    MeasErr(k,loop) = norm(x - xt)^2;
                    MeasErrp(k,loop) = norm(x(1:3) - xt(1:3))^2; 
                    MeasErrv(k,loop) = norm(x(4:6) - xt(4:6))^2;
                    TrPn(k) = trace(P);
              
                    end % end prediction
                    
            end % end for KF   
        end % end for of MC runs
        close(h)
        % CONVERGENCE
        % In order to study the convergence of the KF, the trace of Pn has to 
        % tend towards the value of the MSE for each sampling time
        MeasMSE = sqrt((1/nMC)*sum(MeasErr,2));
        MeasMSEp = sqrt((1/nMC)*sum(MeasErrp,2)); 
        MeasMSEv = sqrt((1/nMC)*sum(MeasErrv,2)); 
        RMSE_sigma(i,j) = MeasMSE(end); % sensitivity analysis
    end
end
toc


%% PLOTS
figure
surf(e1*T,e2*T,10*log10(RMSE_sigma),'Linewidth',2)
grid on
xlabel('Realisation of $\sigma_{pk}$ [m]','Interpreter','latex','fontsize',16)
ylabel('Realisation of $\sigma_{vk}$ [m/s]','Interpreter','latex','fontsize',16)
zlabel('$\sqrt{MSE}$ [dB]','Interpreter','latex','fontsize',16)
title('Sensitivity of uncertainties on evolution model for the Doppler EKF','Interpreter','latex','fontsize',16)

%% SAVE PLOTS FOR LaTeX
%print -deps ClassicalEKF_no_mDT