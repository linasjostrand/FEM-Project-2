function [KinE, IntE] = plan3gEd(Me,vel,ep,eff)
% PLAN3gEd Calculates kinetic and internal energy for a triangular 3-node
% element
% Me element mass matrix
% vel velocity
% ep=[E v] material properties 
% eff deformation gradient

% Kinetic energy
KinE = 0.5*vel'*Me*vel;

% Strain energy
F = [eff(1:2)' 0; eff(3:4)' 0; 0 0 1];
E = ep(1); v = ep(2);
K = E/(3*(1-2*v)); G = E/(2*(1+v));
J = det(F);
C = F'*F;

IntE = 0.5*K*(0.5*(J^2-1) - log(J)) + 0.5*G*(J^(-2/3)*trace(C) - 3);

end