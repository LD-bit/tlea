% AUTHOR:   DUBREIL Léa
% DATE:     16/02/2023 (created) 07/03/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Baby Extended Kalman Filter of a non linear observation and 
%           non linear evolution model
% REF:      Samy Labsir, 2020, Méthodes statistiques fondées sur les 
%           groupes de Lie pour le suivi d'un amas de débris spatiaux.

clear all; close all; clc;

%% ---------------------------- VARIABLES ZONE ----------------------------
% DISCRETISATION OF TIME
T = 1e-4;               %[s] sampling period: measure every second
N = 60*5;               %[s] max duration of measurements
step = N/T;             %[-] step
t = 0:T:N*T-T;          % time vector
nMC = 1e2;              % Monte Carlo number (keep it min 100)

% matrice for uncertainties, order of magnitude
% uncertainties for speed ~=0.1 m/s and for pos ~=10m
M = [eye(3) zeros(3);zeros(3) 0.01*eye(3)]; 

% EVOLUTION MODEL
% Uncertainties in model
sigma_vk = 1e-3; % standard deviations vel
sigma_pk = 1e-4; % pos
% Initial State Vector
% the initialisation of this vector is for ops-sat state vec, on Feb 20,
% 2023 (TLE age was 0.4 days). Position in [m] and velocity in [m/s]
x0 = [-3089022.5178; %[m] x
      -3952006.3867; %[m] y
      4687365.0410;  %[m] z
      2105.6276;     %[m/s] vx
      4864.7572;     %[m/s] vy
      5469.7445];    %[m/s] vz
%x0 = [200 0 0 0 0 0]';
% OBSERVATION MODEL
% Uncertainties in measures
sigma_nk = 1; %[m] standard deviation range
sigma_tk = 1e-5; % theta [rad]
sigma_lk = 1e-5; % phi [rad]
% VARIABLES 
% fs  = 70e6;             %[Hz] satellite frequency (emitted)
% c   = 299792458;        %[m/s] speed of light
% Om  = 7.292115e-5;      %[rad/s] rotation rate of Earth
% G   = 6.6743015e-11;    %[m^3/(kg.s^2)] gravitational constant
% M   = 5.9722e24;        %[kg] Earth mass

% PARAMETERS SAVED
xt = zeros(6,N,nMC); % true state of x
z = zeros(3,N,nMC); % measurements z
x = zeros(6,1); % state vector
MeasErr = zeros(N,nMC); % measurement error
EstErr = zeros(N,nMC); % estimation error
TrPn = zeros(1,N); % Trace of Pn

%% ------------------------------ MAIN ZONE -------------------------------

h = waitbar(0,'Please wait...');

% EXTENDED KALMAN FILTER
for loop = 1:nMC

    waitbar(loop/nMC,h);

    % INITIALISATION
    Rk = diag([sigma_nk^2;sigma_tk^2;sigma_lk^2]); % covariance matrix of measurements
    Qk = diag([sigma_pk^2*ones(3,1);sigma_vk^2*ones(3,1)]); % covariance matrix of evaluations
    P = [100*eye(3) zeros(3);
          zeros(3) eye(3)]; % covariance P[n|n]

    for k = 1:N
        
        % INITIALISATION: for each sampling time
        % noise n_k ~ N(0,Rk)
        n_k = chol(Rk)'*randn(3,1); 
        % noise v_k ~ N(0,Qk)
        %v_k = sigma_vk^2*reshape(100*[randn(3,1) 0.1*randn(3,1)],6,1); 
        v_k = chol(Qk)'*randn(6,1);
        
        % TRUE STATES
        if k == 1
            xt(:,k,loop) = x0; % initialisation
        else
            xt(:,k,loop) = xt(:,k-1,loop) + fk(xt(:,k-1,loop))*T + v_k; % eq. 1.45, true state
        end

        % MEASUREMENT
        z(:,k,loop) = hk(xt(:,k,loop)') + n_k; % eq. 1.46, measurement z[k]

        % PREDICTION 
        if k == 1
            % initialisation of x[n|n] 
            % NB: it is != of init state (true state):
            x = x0 + randn(6,1).*[10 10 10 1e-3 1e-3 1e-3]';
            MeasErr(k,loop) = norm(x - xt(:,k,loop));
        else
            P = Jfk(x,T)*P*Jfk(x,T)' + Qk; 
            x = x + fk(x)*T; % eq. 1.50 state estimation x[k|k-1]
        
           % eq. 1.51, covariance estimation P[k|k-1]
            
            % CORRECTION
            Sk = Jhk(x)*P*Jhk(x)' + Rk; % eq. 1.56 Update matrix
            Kn = P*Jhk(x)'/Sk; % eq. 1.55, Kalman gain K[n]
            x = x + Kn*(z(:,k,loop) - hk(x)); % eq. 1.53, state update x[k|k]
            P = (eye(6) - Kn*Jhk(x))*P; % eq. 1.54, covariance update P[k|k]
            
            % ERRORS
            % norm of the diff for each state between true and filtered
            MeasErr(k,loop) = norm(x - xt(:,k,loop)); 
            % norm of the diff for each state between estimated and filtered
            % EstErr(k,loop) = norm(x- xf(:,k,loop))^2;
            TrPn(k) = trace(P);
        end
        
    end % end for KF   
end % end for of MC runs

close(h)

% CONVERGENCE
% In order to study the convergence of the KF, the trace of Pn has to 
% tend towards the value of the MSE for each sampling time
MeasMSE = (1/nMC)*sum(MeasErr,2); % Mean square error on measurements
%EstMSE = (1/nMC)*sum(EstErr,2); % Mean square error on estimation

%% PLOTS
figure
plot(t,MeasMSE,'-r')
xlabel('Time [s]'), ylabel('MSE')
title('Study of Extended KF convergence')
legend('Mean Square Error')