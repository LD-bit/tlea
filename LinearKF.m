% AUTHOR:   DUBREIL Léa
% DATE:     08/02/2023 (created) 06/03/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Kalman Filter of a linear observation and evolution model
% REF:      Samy Labsir, 2020, Méthodes statistiques fondées sur les 
%           groupes de Lie pour le suivi d'un amas de débris spatiaux.
clear all; close all; clc;

%% ---------------------------- VARIABLES ZONE ----------------------------
% DISCRETISATION OF TIME
T = 0.25;               %[s] sampling period: measure every second
N = 60*5;               %[s] max duration of measurements
step = N/T;             %[-] step
t = 0:T:N*T-T;          % time vector
nMC = 1e3;              % Monte Carlo number (keep it min 100)

% matrice for uncertainties, order of magnitude
% uncertainties for speed ~=0.1 m/s and for pos ~=10m
M = [eye(3) zeros(3);zeros(3) 0.01*eye(3)]; 

% EVOLUTION MODEL
% Uncertainties in model
sigma_vk = 1; % standard deviations
Fk = [eye(3) T*eye(3);zeros(3) eye(3)]; % transition matrix

% Initial State Vector
% the initialisation of this vector is for ops-sat state vec, on Feb 20,
% 2023 (TLE age was 0.4 days). Position in [m] and velocity in [m/s]
x0 = [-3089022.5178; %[m] x
      -3952006.3867; %[m] y
      4687365.0410;  %[m] z
      2105.6276;     %[m] vx
      4864.7572;     %[m] vy
      5469.7445];    %[m] vz

% OBSERVATION MODEL
% Uncertainties in measures
sigma_nk = 0.5; % standard deviation
% observation matrix
Hk = [eye(3) zeros(3);zeros(3) eye(3)]; % Observing only positions

% PARAMETERS SAVED
x = zeros(6,N,nMC); % true state of x
z = zeros(6,N,nMC); % measurements z
xe = zeros(6,N,nMC); % predicted state of xk
xf = zeros(6,N,nMC); % filtered state xk
MeasErr = zeros(N,nMC); % measurement error
EstErr = zeros(N,nMC); % estimation error
TrPn = zeros(1,N); % Trace of Pn

%% ------------------------------ MAIN ZONE -------------------------------

h = waitbar(0,'Please wait...');

% KALMAN FILTER
for loop = 1:nMC

    waitbar(loop/nMC,h);

    % INITIALISATION
    Rk = sqrt(sigma_nk)*M; % covariance matrix of measurements
    Qk = sqrt(sigma_vk)*M; % covariance matrix of evaluations
    Pn = [100*eye(3) zeros(3);
          zeros(3) 0.1*eye(3)]; % covariance P[n|n]

    for k = 1:N
        
        % INITIALISATION: noises for each sampling time
        % noise n_k ~ N(0,Rk)
        n_k = sqrt(sigma_nk)*reshape([randn(3,1) 0.1*randn(3,1)],6,1); 
        % noise v_k ~ N(0,Qk)
        v_k = sqrt(sigma_vk)*reshape([randn(3,1) 0.1*randn(3,1)],6,1); 
        
        % TRUE STATES
        if k == 1
            x(:,k,loop) = x0; % initialisation
        else
            x(:,k,loop) = Fk*x(:,k-1,loop) + v_k; % eq. 1.39, true state
        end

        % MEASUREMENT
        z(:,k,loop) = Hk*x(:,k,loop) + n_k; % eq. 1.40, measurement z[n]

        % PREDICTION 
        if k == 1
            % initialisation of x[n|n] 
            % NB: it is != of init state (true state):
            xe(:,k,loop) = x0 + reshape([randn(3,1) 0.1*randn(3,1)],6,1); 
        else 
            xe(:,k,loop) = Fk*xf(:,k-1,loop);   % eq. 1.41, state estimation x[n+1|n]
        end
        Pk = Fk*Pn*Fk' + Qk;  % eq. 1.42, covariance estimation P[n+1|n]
        
        % CORRECTION
        Sk = Rk + Hk*Pk*Hk'; % eq. p.37 Update covariance matrix S[n]
        Kn = Pk*Hk'/Sk; % eq. p.37, Kalman gain K[n]
        xf(:,k,loop) = xe(:,k,loop) + Kn*(z(:,k,loop) - Hk*xe(:,k,loop)); % eq. 1.44, state update x[n|n]
        Pn = (eye(6) - Kn*Hk)*Pk; % eq. 1.43, covariance update P[n|n]
        
        % ERRORS
        % norm of the diff for each state between true and filtered
        MeasErr(k,loop) = norm(xf(:,k,loop) - x(:,k,loop))^2; 
        % norm of the diff for each state between estimated and filtered
        EstErr(k,loop) = norm(xe(:,k,loop) - xf(:,k,loop))^2;
        TrPn(k) = trace(Pn);

    end % end for KF   
end % end for of MC runs
close(h)

% CONVERGENCE
% In order to study the convergence of the KF, the trace of Pn has to 
% tend towards the value of the MSE for each sampling time
MeasMSE = (1/nMC)*sum(MeasErr,2); % Mean square error on measurements
EstMSE = (1/nMC)*sum(EstErr,2); % Mean square error on estimation

% PLOTS
figure
plot(t,MeasMSE,'-r',t,TrPn','b')
xlabel('Time [s]'), ylabel('MSE')
title('Study of Linear KF convergence')
legend('Measurement Mean Square Error','Trace of covariance matrix')

% figure
% plot(t,EstMSE,'-g')
% xlabel('Time [s]'), ylabel('MSE')
% title('Study of Linear KF convergence')
% legend('Estimation Mean Square Error')