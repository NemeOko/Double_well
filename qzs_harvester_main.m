% QZS_HARVESTER_MAIN  Frequency response of a bimorph piezoelectric cantilever
% energy harvester, with or without a magnetic quasi-zero-stiffness (QZS) element.
%
%   Model: Euler-Bernoulli beam, single-mode expansion, piezoelectric coupling,
%   resistive load, point-dipole magnetic force. Solved with ode45.
%   Requires: MATLAB with the Symbolic Math Toolbox. Run from this folder.
%
%   Prompts: magnetic force on/off, frequency range and step, load resistances.
%   For a reverse sweep, enter a minimum above the maximum with a negative step.
%
%   Author: Neme Okoroafor, University of Manchester (MEng dissertation, 2025)
clear;

%% ============================ USER INPUT ============================
Mag_force = input('Magnetic force included? (1 = Yes): ');
calc_dis = 1;

min_freq = input('Minimum frequency to be assessed: ');
max_freq = input('Maximum frequency to be assessed: ');
int_freq = input('Frequency interval (1 if single value): ');

Frequency_range = min_freq:int_freq:max_freq;
w_range = Frequency_range * 2 * pi;

Resistance = input('Resistances in ohms to be assessed (vector): ');

clear min_freq max_freq;

%% ========================= BEAM & MATERIAL PROPERTIES ===============
% Beam dimensions (m)
L  = 60e-3;
hp = 0.267e-3;
hs = 0.3e-3;
b  = 25e-3;

% Material properties
Ys   = 72e9;        % Young's modulus (substrate)
Yp   = 62e9;        % Young's modulus (piezoelectric)
rhop = 7800;        % Density (piezoelectric)
rhos = 2700;        % Density (substrate)
d31  = -320e-12;    % Piezoelectric constant
ep33t = 3800 * 8.854187817e-12;     % Permittivity at constant stress
ep33s = ep33t - (d31^2) * Yp;       % Permittivity at constant strain
xi    = 1.09e-2;    % Damping ratio of identical beam with no tip mass
Omega = -(d31 * Yp * b * (hp + hs)) / 2;
m     = b * (2 * hp * rhop + hs * rhos);   % Mass per unit length
mt    = 4.3e-3;     % Tip magnet mass (kg)

% Derived beam constants
YI = b * (((Ys * hs^3) / 12) + Yp * (((hp * hs^2) / 2) + (hp^2 * hs) + ...
    (2/3) * hp^3));                 % Flexural stiffness of the whole beam
Cp = (ep33s * b * L) / hp;          % Capacitance
zc = (hp + hs) / 2;                 % Mid-point depth of the beam

%% ============================ MODAL ANALYSIS ========================
lambdar0 = 1.87510407;  % Clamped-free beam eigenvalue with no tip mass

% Eigenvalue with tip mass
lambdar = fzero(@(lam) find_lambda_r(lam, mt, m*L), lambdar0);
sigmar = (sinh(lambdar) - sin(lambdar) + lambdar * (mt / (m * L)) * ...
         (cosh(lambdar) - cos(lambdar))) / ...
         (cosh(lambdar) + cos(lambdar) + lambdar * (mt / (m * L)) * ...
         (sinh(lambdar) - sin(lambdar)));

% Mode shape and integrals (symbolic)
syms x;
psir_x = cosh((lambdar/L)*x) - cos((lambdar/L)*x) ...
        - sigmar * (sinh((lambdar/L)*x) - sin((lambdar/L)*x));
psir_x_sqr = psir_x^2;

int_psi_r_sqr = double(int(psir_x_sqr, x, 0, L));
psir_x_sqr_L = (cosh(lambdar) - cos(lambdar) ...
               - sigmar * (sinh(lambdar) - sin(lambdar)))^2;

Mr = m * int_psi_r_sqr + mt * psir_x_sqr_L;
phir_x = sqrt(1 / Mr) * psir_x;     % Mass-normalised mode shape

phir_L = sqrt(1 / Mr) * (cosh(lambdar) - cos(lambdar) ...
          - sigmar * (sinh(lambdar) - sin(lambdar)));

int_phi_r      = double(int(phir_x, x, 0, L));
diff_phi_r_L   = double(subs(diff(phir_x, x), x, L));
int_diff_phir_L  = double(int(diff(phir_x, x), x, 0, L));
int_2diff_phir_L = double(int((diff(diff(phir_x, x), x))^2, x, 0, L));

%% ============================ MODEL CONSTANTS =======================
wr  = lambdar^2  * sqrt(YI / (m * L^4));    % Natural frequency with tip mass
wr0 = lambdar0^2 * sqrt(YI / (m * L^4));    % Natural frequency, no tip mass
xi_mass = wr * (xi / wr0);                  % Damping term with tip mass

Xr     = Omega * diff_phi_r_L;                      % Piezoelectric coupling term
alph_r = -Yp * d31 * b * zc * diff_phi_r_L;

% Constants for the equation of motion
A1 = -2 * xi_mass * wr;
A2 = -wr^2;
A3 = -m * int_phi_r - mt * phir_L;
A4 = phir_L;
B1 = 2 / Cp;

%% ========== QZS OPTIMISATION VIA MAGNETIC FORCE TUNING ==============
% Magnet properties (point-dipole approximation)
Sa = 200e-6;            % Surface area (m^2)
Sb = Sa;
Ma = 0.96e6;            % Magnetisation (A/m)
Mb = Ma;
u0 = 4 * pi * 1e-7;
a  = (2.75 / 2) * 1e-3; % Half magnet length (m)

if Mag_force == 1
    wa = linspace(-0.02, 0.02, 200);    % Tip deflection range (m)
    best_d = NaN;
    min_derivative = Inf;

    % Scan the magnet gap for the smallest net stiffness at zero deflection
    for dk = 5:0.05:15
        d = dk * 1e-3;

        Fmagy = zeros(size(wa));    % Magnetic tip force
        Fbeam = zeros(size(wa));    % Elastic restoring tip force

        for c = 1:length(wa)
            w  = wa(c);
            du = diff_phi_r_L * w / phir_L;

            Fmag_1 = (Ma * Sa * (-Mb * Sb * w / ((2*a + d)^2 + w^2)^(3/2) ...
                     + Mb * Sb * w / ((4*a + d)^2 + w^2)^(3/2)));
            Fmag_2 = (Ma * Sa * (-Mb * Sb * (w + 2*a*du) / (d^2 + (w + 2*a*du)^2)^(3/2) ...
                     + Mb * Sb * (w + 2*a*du) / ((2*a + d)^2 + (w + 2*a*du)^2)^(3/2)));

            Fmagy(c) = (u0 / (4*pi)) * (Fmag_1 - Fmag_2);
            Fbeam(c) = 3 * (YI / L^3) * w;
        end

        Ftotal = -Fmagy + Fbeam;
        % Potential energy from numerical integration of the magnetic force
        U_mag = -cumtrapz(wa, Fmagy);
        dw = wa(2) - wa(1);
        derivative_values = gradient(Ftotal, dw);

        [~, idx_zero] = min(abs(wa - 0));
        derivative_at_zero = derivative_values(idx_zero);

        if abs(derivative_at_zero) < abs(min_derivative)
            min_derivative = derivative_at_zero;
            best_d = d;
            best_Fmagy = Fmagy;
            best_U_mag = U_mag;
        end
    end

    % Static elastic potential energy of the beam
    U_beam  = ((wa * (1/phir_L)).^2) * 0.5 * YI * int_2diff_phir_L;
    U_total = U_beam + best_U_mag;
    Ftotal  = -best_Fmagy + Fbeam;

    fprintf('Optimal magnet gap (d): %.4f m\n', best_d);
    fprintf('Derivative at wa=0: %.4e\n', min_derivative);

    figure;
    plot(wa, Fbeam, 'k', 'LineWidth', 0.5)
    hold on
    plot(wa, best_Fmagy, 'k--', wa, Ftotal, 'k', 'LineWidth', 1.5);
    hold off
    xlabel('Tip Deflection (m)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    ylabel('F_{restoring} (N)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    legend('Beam', 'Magnet', 'Total', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Cambria')
    set(gca, 'FontSize', 12, 'FontName', 'Cambria');

    figure;
    plot(wa, U_beam, 'k', 'LineWidth', 0.5)
    hold on
    plot(wa, best_U_mag, 'k--', wa, U_total, 'k', 'LineWidth', 1.5);
    hold off
    legend('Beam', 'Magnet', 'Total', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Cambria');
    xlabel('Tip Deflection (m)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    ylabel('Potential Energy (J)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    set(gca, 'FontSize', 12, 'FontName', 'Cambria');
else
    best_d = 0;
end

%% ========================== RUNNING THE MODEL =======================
volt = zeros(length(Frequency_range), length(Resistance));
if calc_dis == 1
    tip_dis = zeros(length(Frequency_range), length(Resistance));
end

% Indices of the two frequencies (4 Hz and 70 Hz) used for time histories
[~, idx_zero]  = min(abs(w_range - (4*2*pi)));
[~, idx_zero2] = min(abs(w_range - (70*2*pi)));

y0 = [0, 0, 0];     % Initial state [eta, etadot, v]; carried over between steps
for iR = 1:length(Resistance)
    R = Resistance(iR);
    for c = 1:length(w_range)
        w = w_range(c);
        [t, sol] = ode45(@(t, q) piezo_beam_resistor(t, q, A1, A2, A3, B1, R, w, ...
            alph_r, Xr, A4, diff_phi_r_L, Mag_force, best_d), [0 4], y0);
        y0 = sol(end, :);

        % Amplitude measured over the final part of the run (t > 3.8 s)
        idxStart = sum(t <= 3.8);
        volt(c, iR) = (max(sol(idxStart:end, 3)) - min(sol(idxStart:end, 3))) / 2;

        if calc_dis == 1
            tip_dis(c, iR) = (max(sol(idxStart:end, 1)) - min(sol(idxStart:end, 1))) / 2;
            disp(w / (2 * pi))
        end

        if ismember(c, [idx_zero, idx_zero2])
            xdot  = sol(:, 2) * phir_L;
            x_dis = sol(:, 1) * phir_L;
            figure;
            plot(t, x_dis * 1e3);
            xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
            ylabel('Tip Displacement (mm)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
            set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria')

            figure;
            plot(x_dis * 1e3, xdot);
            xlabel('Tip Displacement (mm)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
            ylabel('Velocity (m/s)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria')
            set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria')

            if Mag_force == 1
                % Dynamic potential energy along the trajectory (magnet case)
                eta    = sol(:, 1);
                etadot = sol(:, 2);
                vol    = sol(:, 3);
                wa_t   = eta * phir_L;
                du_t   = diff_phi_r_L * eta;
                Fmag_y = zeros(size(wa_t));

                for k = 1:length(wa_t)
                    wk = wa_t(k);
                    Fmag1 = (Ma * Sa * (-Mb * Sb * wk / ((2*a + best_d)^2 + wk^2)^(3/2) ...
                             + Mb * Sb * wk / ((4*a + best_d)^2 + wk^2)^(3/2)));
                    Fmag2 = (Ma * Sa * (-Mb * Sb * (wk + 2*a*du_t(k)) / (best_d^2 + (wk + 2*a*du_t(k))^2)^(3/2) ...
                             + Mb * Sb * (wk + 2*a*du_t(k)) / ((2*a + best_d)^2 + (wk + 2*a*du_t(k))^2)^(3/2)));
                    Fmag_y(k) = (u0 / (4*pi)) * (Fmag1 - Fmag2);
                end
                Umag = -cumtrapz(wa_t, Fmag_y);
                vdot = gradient(vol, t);
                Ubeam = (eta.^2) * 0.5 * YI * int_2diff_phir_L ...
                    + (eta .* (alph_r*etadot - ((Cp/2) * gradient(vol, t)))) * 2*Omega*R*int_diff_phir_L ...
                    - ((Cp/2) * (R * (alph_r*etadot - ((Cp/2) * vdot))).^2);
                Utotal = Umag + Ubeam;

                figure
                plot(wa_t * 1e3, Utotal, 'b', 'LineWidth', 1.5)
                hold on
                plot(wa_t * 1e3, Umag, 'r--', 'LineWidth', 1.5);
                plot(wa_t * 1e3, Ubeam, 'g--', 'LineWidth', 1.5);
                hold off
                legend('Total', 'Magnet', 'Beam', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Cambria');
                xlabel('Tip Deflection (mm)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
                ylabel('Dynamic Potential Energy (J)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
                set(gca, 'FontSize', 12, 'FontName', 'Cambria');
            end
        end
    end
    disp(['R = ', num2str(R), ' Ohms']);
end

if calc_dis == 1
    tip_dis = tip_dis * phir_L;     % Convert modal coordinate to tip displacement
end

%% ========================= PLOTTING RESULTS =========================
% Voltage frequency response
if input('Plot Voltage FRF (1 = Yes): ') == 1
    figure;
    for iR = 1:length(Resistance)
        semilogy(Frequency_range, volt(:, iR), '-');
        hold on;
    end
    xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    ylabel('Voltage (V)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance, 'UniformOutput', false), ...
        'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Cambria');
    set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria');
end

% Tip displacement frequency response
if input('Plot Tip Displacement FRF (1 = Yes): ') == 1
    figure;
    for iR = 1:length(Resistance)
        plot(Frequency_range, tip_dis(:, iR));
        hold on;
    end
    xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    ylabel('Tip Displacement (m)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Cambria');
    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance, 'UniformOutput', false), ...
        'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria');
    set(gca, 'FontSize', 12, 'FontName', 'Cambria');
end
