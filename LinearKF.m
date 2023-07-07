% AUTHOR:   DUBREIL Léa
% DATE:     08/02/2023 (created) 04/04/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Kalman Filter of a linear observation and evolution model
% REF:      Samy Labsir, 2020, Méthodes statistiques fondées sur les 
%           groupes de Lie pour le suivi d'un amas de débris spatiaux.
clear all; close all; clc;

%% ---------------------------- VARIABLES ZONE ----------------------------
% DISCRETISATION OF TIME
T = 0.25;               %[s] sampling period: measure every second
N = 60*8;               %[s] max duration of measurements
step = N/T;             %[-] step
t = 0:T:N*T-T;          %[s] time vector
nMC = 1e3;              %[-] Monte Carlo number (keep it min 100)

% EVOLUTION MODEL
% Uncertainties in model
sigma_vk = 1e1;         %[m/s] Velocity standard deviation
sigma_pk = 1e-2;        %[m] Position standard deviation
Fk = [eye(3) T*eye(3);zeros(3) eye(3)]; % transition matrix

% Initial State Vector
% the initialisation of this vector is for ops-sat state vec, on Feb 20,
% 2023 (TLE age was 0.4 days). Position in [m] and velocity in [m/s]
x0 = [-3089022.5178; %[m] x
      -3952006.3867; %[m] y
      4687365.0410;  %[m] z
      2105.6276;     %[m/s] vx
      4864.7572;     %[m/s] vy
      5469.7445];    %[m/s] vz

% OBSERVATION MODEL
% Uncertainties in measures
sigma_nk    = 1e0;    %[m] standard deviation position
sigma_nnk   = 1e-1;   %[m/s] standard deviation velocity
Hk = [eye(3) zeros(3);zeros(3) eye(3)]; % observation matrix

% PARAMETERS SAVED
x = zeros(6,1); % true state of x
z = zeros(6,N); % measurements z
MeasErr = zeros(N,nMC); % measurement error
TrPn = zeros(1,N); % Trace of Pn

%% ------------------------------ MAIN ZONE -------------------------------
disp('[INFO] Kalman Filter: linear case')
h = waitbar(0,'Please wait...');

% KALMAN FILTER
for loop = 1:nMC

    waitbar(loop/nMC,h);

    % INITIALISATION
    Rk = diag([sigma_nk^2*ones(3,1);sigma_nnk^2*ones(3,1)]); % covariance matrix of measurements
    Qk = diag([sigma_pk^2*ones(3,1);sigma_vk^2*ones(3,1)]); % covariance matrix of evaluations
    Pn = [1*eye(3) zeros(3);
          zeros(3) eye(3)]; % covariance P[n|n]
    
    for k = 1:N
        
        % INITIALISATION: noises for each sampling time
        % noise n_k ~ N(0,Rk)
        n_k = chol(Rk)'*randn(6,1); 
        % noise v_k ~ N(0,Qk)
        v_k = chol(Qk)'*randn(6,1);
        
        % TRUE STATES
        if k == 1
            xt = x0; % initialisation
        else
            xt = Fk*xt + v_k; % eq. 1.39, true state
        end

        % MEASUREMENT
        z = Hk*xt + n_k; % eq. 1.40, measurement z[n]

        % PREDICTION 
        if k == 1
            % initialisation of x[n|n] 
            % NB: it is != of init state (true state):
            x = x0 + reshape([randn(3,1) 0.1*randn(3,1)],6,1); 
            % ERRORS
            % norm of the diff for each state between true and filtered
            MeasErr(k,loop) = norm(x - xt)^2; 
            TrPn(k) = trace(Pn);
        else 
            x = Fk*x;   % eq. 1.41, state estimation x[n+1|n]
        
            Pk = Fk*Pn*Fk' + Qk;  % eq. 1.42, covariance estimation P[n+1|n]
            
            % CORRECTION
            Sk = Rk + Hk*Pk*Hk'; % eq. p.37 Update covariance matrix S[n]
            Kn = Pk*Hk'/Sk; % eq. p.37, Kalman gain K[n]
            x = x + Kn*(z - Hk*x); % eq. 1.44, state update x[n|n]
            Pn = (eye(6) - Kn*Hk)*Pk; % eq. 1.43, covariance update P[n|n]
            
            % ERRORS
            % norm of the diff for each state between true and filtered
            MeasErr(k,loop) = norm(x - xt)^2; 
            TrPn(k) = trace(Pn);
        end
    end % end for KF   
end % end for of MC runs
close(h)

% CONVERGENCE
% In order to study the convergence of the KF, the trace of Pn has to 
% tend towards the value of the MSE for each sampling time
MeasMSE = (1/nMC)*sum(MeasErr,2); % Mean square error on measurements

%% PLOTS
figure
plot(t,MeasMSE,'-sr','Linewidth',2), hold on
plot(t,TrPn','-.ob','Linewidth',2)
grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('[-]','Interpreter','latex','fontsize',16)
title('Convergence of Linear KF','Interpreter','latex','fontsize',16)
legend('Total $MSE$','$Tr(P_{n})$','Interpreter','latex','fontsize',16)