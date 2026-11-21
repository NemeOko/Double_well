function f = find_lambda_r(lambdar, mt, mb)
% FIND_LAMBDA_R  Characteristic equation of a clamped-free beam with a tip mass.
%   Roots give the eigenvalues lambda_r. mt = tip mass, mb = beam mass (m*L).
f = 1 + cos(lambdar)*cosh(lambdar) + (mt/mb)*lambdar * ...
    (cos(lambdar)*sinh(lambdar) - sin(lambdar)*cosh(lambdar));
end
