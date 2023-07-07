% AUTHOR:   DUBREIL Léa
% DATE:     09/03/2023 (created) 06/04/2023 (modified)
% PROJECT:  Master Thesis
% NAME:     Comparison of Classical EKF and Doppler EKF (main contribution)
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
sigma_pk = 1e-0;        %[m] Position standard deviation

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
z = zeros(4,N,nMC); % measurements z
x = zeros(6,1); % state vector
MeasErrC = zeros(N,nMC); MeasErrD = MeasErrC; % MSE total
MeasErrCp = zeros(N,nMC); MeasErrDp = MeasErrC; % MSE position
MeasErrCv = zeros(N,nMC); MeasErrDv = MeasErrC; % MSE velocity
MeasErrCvx = zeros(N,nMC); MeasErrCvy = MeasErrCvx; MeasErrCvz = MeasErrCvx;
MeasErrDvx = zeros(N,nMC); MeasErrDvy = MeasErrDvx; MeasErrDvz = MeasErrCvx;
%% ------------------------------ DOPPLER EKF -------------------------------
disp('[INFO] Comparison of Classical EKF and Doppler EKF (main contrib)')
h = waitbar(0,'Please wait...');

% EXTENDED KALMAN FILTER
for loop = 1:nMC

    waitbar(loop/nMC,h);

    % INITIALISATION
    % covariance matrix of measurements for classical method
    Rkc = diag([sigma_nk^2;sigma_tk^2;sigma_lk^2]);
    % covariance matrix for Doppler integration
    Rkd = diag([sigma_nk^2; ... % range1
                sigma_tk^2; ... % az1
                sigma_lk^2; ... % el1
                sigma_fk1^2; ... % fd1
                %sigma_tk^2; ... % az2
                %sigma_lk^2; ... % el2
                sigma_fk2^2; ... % fd2
                %sigma_vvk^2;...
                sigma_fk3^2]); % fd3
    % covariance matrix of evaluations
    Qk = diag([sigma_pk^2*ones(3,1);sigma_vk^2*ones(3,1)]); 
    % covariance P[n|n]
    Pc = [100*eye(3) zeros(3);
          zeros(3) eye(3)]; 
    Pd = Pc;

    for k = 1:N
        
        % INITIALISATION: for each sampling time
        % noise n_k ~ N(0,Rk)
        n_kc = chol(Rkc)'*randn(3,1); 
        n_kd = chol(Rkd)'*randn(6,1);
        % noise v_k ~ N(0,Qk)
        v_k = chol(Qk)'*randn(6,1);
        
        for g = 0:(mdT/T)-1 % loop to implement the delta of measurement
            % TRUE STATES
            if k == 1
                xt = x0; % initialisation for k=1
            else
                xt = xt + fk(xt)*T + v_k; % eq. 1.45, true state
            end
        end
        % MEASUREMENT
        zc = hk(xt) + n_kc; % eq. 1.46, measurement z[k]
        zd = Doppler_hk(xt,xt) + n_kd; % eq. 1.46, measurement z[k]
        
        % PREDICTION 
        if k == 1
            % initialisation of x[n|n] 
            % NB: it is != of init state (true state):
            xc = x0 + randn(6,1).*[50 50 50 5 5 5]/2';
            xd = xc;%x0 + randn(6,1).*[50 50 50 5 5 1]';
            % ERRORS
            % Classic
            MeasErrC(k,loop) = norm(xc - xt)^2;
            MeasErrCp(k,loop) = norm(xc(1:3) - xt(1:3))^2;
            MeasErrCvx(k,loop) = norm(xc(4) - xt(4))^2;
            MeasErrCvy(k,loop) = norm(xc(5) - xt(5))^2;
            MeasErrCvz(k,loop) = norm(xc(6) - xt(6))^2;
            % Doppler
            MeasErrD(k,loop) = norm(xd - xt)^2;
            MeasErrDp(k,loop) = norm(xd(1:3) - xt(1:3))^2;
            MeasErrDvx(k,loop) = norm(xd(4) - xt(4))^2;
            MeasErrDvy(k,loop) = norm(xd(5) - xt(5))^2;
            MeasErrDvz(k,loop) = norm(xd(6) - xt(6))^2;
        else
            for g = 0:(mdT/T)-1
            % CLASSIC
                Pc = Jfk(xc,T)*Pc*Jfk(xc,T)' + Qk; % eq. 1.51, covariance estimation P[k|k-1]
                xc = xc + fk(xc)*T; % eq. 1.50 state estimation x[k|k-1]
            % DOPPLER
                Pd = Jfk(xd,T)*Pd*Jfk(xd,T)' + Qk; % eq. 1.51, covariance estimation P[k|k-1]
                xd = xd + fk(xd)*T; % eq. 1.50 state estimation x[k|k-1]
            end
            
            % CORRECTION CLASSIC
            Skc = Jhk(xc)*Pc*Jhk(xc)' + Rkc; % eq. 1.56 Update matrix
            Knc = Pc*Jhk(xc)'/Skc; % eq. 1.55, Kalman gain K[n]
            xc = xc + Knc*(zc - hk(xc)); % eq. 1.53, state update x[k|k]
            Pc = (eye(6) - Knc*Jhk(xc))*Pc; % eq. 1.54, covariance update P[k|k]
            % CORRECTION DOPPLER
            Skd = Doppler_Jhk(xd,xt)*Pd*Doppler_Jhk(xd,xt)' + Rkd; % eq. 1.56 Update matrix
            Knd = Pd*Doppler_Jhk(xd,xt)'/Skd; % eq. 1.55, Kalman gain K[n]
            xd = xd + Knd*(zd - Doppler_hk(xd,xt)); % eq. 1.53, state update x[k|k]
            Pd = (eye(6) - Knd*Doppler_Jhk(xd,xt))*Pd; % eq. 1.54, covariance update P[k|k]
            
            % ERRORS
            % Classic
            MeasErrC(k,loop) = norm(xc - xt)^2;
            MeasErrCp(k,loop) = norm(xc(1:3) - xt(1:3))^2;
            MeasErrCvx(k,loop) = norm(xc(4) - xt(4))^2;
            MeasErrCvy(k,loop) = norm(xc(5) - xt(5))^2;
            MeasErrCvz(k,loop) = norm(xc(6) - xt(6))^2;
            % Doppler
            MeasErrD(k,loop) = norm(xd - xt)^2;
            MeasErrDp(k,loop) = norm(xd(1:3) - xt(1:3))^2;
            MeasErrDvx(k,loop) = norm(xd(4) - xt(4))^2;
            MeasErrDvy(k,loop) = norm(xd(5) - xt(5))^2;
            MeasErrDvz(k,loop) = norm(xd(6) - xt(6))^2;
        end
        
    end % end for KF   
end % end for of MC runs

close(h)

% RMSE
% CLASSIC
MeasMSEC = sqrt((1/nMC)*sum(MeasErrC,2)); 
MeasMSECp = sqrt((1/nMC)*sum(MeasErrCp,2)); 
MeasMSECvx = sqrt((1/nMC)*sum(MeasErrCvx,2));
MeasMSECvy = sqrt((1/nMC)*sum(MeasErrCvy,2));
MeasMSECvz = sqrt((1/nMC)*sum(MeasErrCvz,2));
% DOPPLER
MeasMSED = sqrt((1/nMC)*sum(MeasErrD,2));
MeasMSEDp = sqrt((1/nMC)*sum(MeasErrDp,2)); 
MeasMSEDvx = sqrt((1/nMC)*sum(MeasErrDvx,2));
MeasMSEDvy = sqrt((1/nMC)*sum(MeasErrDvy,2));
MeasMSEDvz = sqrt((1/nMC)*sum(MeasErrDvz,2));

%% PLOTS
% TOTAL RMSE and POSITION RMSE
figure
% MSE total
subplot(2,1,1),plot(t,[MeasMSEC,MeasMSED],'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ [dB]','Interpreter','latex','fontsize',16)
title('RMSE Total','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
% MSE position
subplot(2,1,2),plot(t,10*log10([MeasMSECp,MeasMSEDp]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ in m [dB]','Interpreter','latex','fontsize',16)
title('RMSE in position','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
sgtitle({['Comparison between Classical and Doppler method with sampling time $T_s =$ ' num2str(T) ' [s] and measurements every' ' ' num2str(mdT) ' [s]']},'Interpreter','latex','fontsize',20)
% VELOCITY 3-AXIS RMSE (dB)
figure
% x-axis
subplot(3,1,1),plot(t,10*log10([MeasMSECvx,MeasMSEDvx]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ in m/s [dB]','Interpreter','latex','fontsize',16)
title('On x-axis (ECEF)','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
% y-axis
subplot(3,1,2),plot(t,10*log10([MeasMSECvy,MeasMSEDvy]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ in m/s [dB]','Interpreter','latex','fontsize',16)
title('On y-axis (ECEF)','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
% z-axis
subplot(3,1,3),plot(t,10*log10([MeasMSECvz,MeasMSEDvz]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ in m/s [dB]','Interpreter','latex','fontsize',16)
title('On z-axis (ECEF)','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
sgtitle({['Velocity: Comparison with sampling time $T_s =$ ' num2str(T) ' [s] and measurements every' ' ' num2str(mdT) ' [s]']},'Interpreter','latex','fontsize',20) 

%% % VELOCITY 3-AXIS RMSE (no dB)
figure
% x-axis
subplot(3,1,1),plot(t,([MeasMSECvx,MeasMSEDvx]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ [m/s]','Interpreter','latex','fontsize',16)
title('On x-axis (ECEF)','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
% y-axis
subplot(3,1,2),plot(t,([MeasMSECvy,MeasMSEDvy]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ [m/s]','Interpreter','latex','fontsize',16)
title('On y-axis (ECEF)','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
% z-axis
subplot(3,1,3),plot(t,([MeasMSECvz,MeasMSEDvz]),'Linewidth',2), grid on
xlabel('Time [s]','Interpreter','latex','fontsize',16)
ylabel('$\sqrt{MSE}$ [m/s]','Interpreter','latex','fontsize',16)
title('On z-axis (ECEF)','Interpreter','latex','fontsize',16)
legend('Classical method','Doppler','Interpreter','latex','fontsize',16)
sgtitle({['Velocity: Comparison with sampling time $T_s =$ ' num2str(T) ' [s] and measurements every' ' ' num2str(mdT) ' [s]']},'Interpreter','latex','fontsize',20) 

%% SAVE PLOTS FOR LaTeX
%print -deps ComparisonEKF_velocity_3axis
%% COMPARISON WITH TLE
%tle = readtle('ops-sat_tle.txt', 'OPS-SAT');