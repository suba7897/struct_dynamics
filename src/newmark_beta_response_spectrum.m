function newmark_beta_response_spectrum()

clc; close all;

disp('-----------------------------------------------------------');
disp('  SDOF RESPONSE SPECTRUM — NEWMARK β METHOD (AVERAGED)     ');
disp('-----------------------------------------------------------');

%% ====================== MASS UNIT MENU ==========================

mass_unit = menu('Select MASS unit:', ...
    'kg (mass)', ...
    'kip (weight)', ...
    'kN (weight)' );

switch mass_unit
    case 1
        m_input = input('Enter mass (kg): ');
    case 2
        W_kip = input('Enter weight (kip): ');
        W_N = W_kip * 4448.22;     
        m_input = W_N / 9.80665; 
    case 3
        W_kN = input('Enter weight (kN): ');
        W_N  = W_kN * 1000;       
        m_input = W_N / 9.80665;
end

%% ====================== DAMPING RATIO ===========================

zeta = input('Enter damping ratio ζ (e.g., 0.05): ');

%% ====================== TIME STEP ===============================

dt = input('Enter time step Δt (seconds): ');

%% ====================== LOAD INPUT ===============================

disp('Select Excel file containing: time | load');
[file,path] = uigetfile({'*.xlsx;*.xls'});
raw = readmatrix(fullfile(path,file));

time = raw(:,1);
load_raw = raw(:,2);
N = length(time);

%% ====================== LOAD UNIT MENU ==========================

load_unit = menu('Select LOAD type:', ...
    'Force: N', ...
    'Force: kN', ...
    'Force: kip', ...
    'Ground Acceleration: m/s^2', ...
    'Ground Acceleration: g');

neg_flag = menu('Apply negative sign for F = - m * a ?', ...
    'YES (-m a)', ...
    'NO (+m a)');

sign_factor = -1 + 2*(neg_flag==2);

%% =============== CONVERT LOAD / ACCELERATION ====================

switch load_unit
    case 1
        p_raw = load_raw;
    case 2
        p_raw = load_raw * 1000;
    case 3
        p_raw = load_raw * 4448.22;
    case 4
        p_raw = sign_factor * m_input * load_raw;
    case 5
        a = load_raw * 9.80665;
        p_raw = sign_factor * m_input * a;
end

%% =================================================================
%          ** USER INPUT: ARRAY OF PERIODS FOR RESPONSE SPECTRUM **
%% =================================================================

T_array = input('Enter period array (e.g., 0.05:0.05:5): ');

% Storage for spectrum
Sd = zeros(length(T_array),1);
Sv = zeros(length(T_array),1);
Sa = zeros(length(T_array),1);

%% =================================================================
%                   ** LOOP OVER EACH PERIOD **
%% =================================================================

for j = 1:length(T_array)

    T = T_array(j);
    omega = 2*pi/T;

    % ------------------- SDOF parameters -------------------
    k = m_input * omega^2;
    c = 2*zeta*sqrt(k*m_input);

    % ------------------- Arrays -----------------------------
    v  = zeros(N,1);
    vd = zeros(N,1);
    vdd = zeros(N,1);

    vdd(1) = (p_raw(1) - c*vd(1) - k*v(1)) / m_input;

    % ------------------- Newmark constants ------------------
    beta  = 1/2;
    gamma = 1/2;

    a0 = 1/(beta*dt^2);
    a1 = gamma/(beta*dt);
    a2 = 1/(beta*dt);
    a3 = 1/(2*beta) - 1;
    a4 = gamma/beta - 1;
    a5 = dt*(gamma/(2*beta) - 1);

    k_eff = k + a1*c + a0*m_input;

    % ------------------ Time stepping -----------------------
    for i = 1:N-1
        
        p_eff = p_raw(i+1) ...
                + m_input*( a0*v(i) + a2*vd(i) + a3*vdd(i) ) ...
                + c*( a1*v(i) + a4*vd(i) + a5*vdd(i) );

        v(i+1) = p_eff / k_eff;

        vdd(i+1) = a0*(v(i+1)-v(i)) - a2*vd(i) - a3*vdd(i);
        vd(i+1)  = vd(i) + dt*((1-gamma)*vdd(i) + gamma*vdd(i+1));
    end

    % ================= STORE PEAK RESPONSES ===================
    Sd(j) = max(abs(v));
    Sv(j) = max(abs(vd));
    Sa(j) = max(abs(vdd));

end

%% =================================================================
%                  ** PLOT RESPONSE SPECTRA **
%% =================================================================

figure;

subplot(3,1,1);
plot(T_array, Sd, 'LineWidth', 2);
xlabel('Period T (s)');
ylabel('Spectral Displacement Sd (m)');
title('Response Spectrum — Displacement');
grid on;

subplot(3,1,2);
plot(T_array, Sv, 'LineWidth', 2);
xlabel('Period T (s)');
ylabel('Spectral Velocity Sv (m/s)');
title('Response Spectrum — Velocity');
grid on;

subplot(3,1,3);
plot(T_array, Sa, 'LineWidth', 2);
xlabel('Period T (s)');
ylabel('Spectral Acceleration Sa (m/s^2)');
title('Response Spectrum — Acceleration');
grid on;

sgtitle('Newmark-β Response Spectrum');

end
