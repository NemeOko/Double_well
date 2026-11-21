function qprime = piezo_beam_resistor(t, q, A1, A2, A3, B1, R, w, alph_r, Xr, ...
    phir_L, diff_phi_r_L, Mag_force, best_d)
% PIEZO_BEAM_RESISTOR  State derivatives for a piezoelectric beam harvester
% with a resistive load, base excitation and optional magnetic tip force.
%
%   State vector:      q(1) = eta (modal displacement)
%                      q(2) = etadot
%                      q(3) = v (voltage across the resistor)
%   Returned:          qprime = [etadot; etaddot; vdot]

if Mag_force == 1
    Sa = 200e-6;            % Magnet surface area (m^2)
    Sb = Sa;
    Ma = 0.96e6;            % Magnetisation (A/m)
    Mb = Ma;
    u0 = 4 * pi * 1e-7;     % Permeability of free space (T m/A)
    a  = (2.75 / 2) * 1e-3; % Half magnet length (m)

    wa = q(1) * phir_L;         % Tip displacement
    d  = best_d;                % Magnet gap
    du = diff_phi_r_L * q(1);   % Tip slope

    Fmag_1 = Ma * Sa * (-Mb * Sb * wa / ((2*a + d)^2 + wa^2)^(3/2) ...
        + Mb * Sb * wa / ((4*a + d)^2 + wa^2)^(3/2));
    Fmag_2 = Ma * Sa * (-Mb * Sb * (wa + 2*a*du) / (d^2 + (wa + 2*a*du)^2)^(3/2) ...
        + Mb * Sb * (wa + 2*a*du) / ((2*a + d)^2 + (wa + 2*a*du)^2)^(3/2));
    Fmagy = (u0 / (4*pi)) * (Fmag_1 - Fmag_2);
else
    Fmagy = 0;
end

qprime(1) = q(2);
ub2dot = 1 * 9.81 * sin(w * t);     % Base acceleration amplitude: 1 g (edit for 0.1 g, 2 g)
qprime(2) = A1*q(2) + A2*q(1) - Xr*q(3) + ub2dot*A3 + phir_L*Fmagy;
qprime(3) = B1 * (alph_r*q(2) - (q(3) / R));
qprime = qprime(:);
end
