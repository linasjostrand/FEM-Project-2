clc; clearvars; close all;

% Load geometry
load input_data_2017.mat

% Modulus of elasticity
E = 10e9; %Pa/m^2
% Poisson's ratio
v = 0.3;
ep = [E v];
k = 1e3; % bar stiffness?

% Thickness, plane strain conditions assumed
t = 10e-3; %m

load_dof = 248;
free_nodes = (1:ndof)';
free_nodes(bc(:,1)) = [];

% Extract coordinates for bars
exb = zeros(nelmb,2);
eyb = zeros(nelmb,2);
ezb = zeros(nelmb,2);
for i = 1:nelmb
    
    exb(i,:) = ecb{i}(1,:);
    eyb(i,:) = ecb{i}(2,:);
    ezb(i,:) = ecb{i}(3,:);
    
end

%% Draw geometry
figure
eldraw2(ex,ey,[1 4 1])
hold on
viscircles([xc1 yc1; xc2 yc2],[r;r]);
xlim([-0.08 0.31])
title('Geometry','Fontsize',13);
xlabel('{\it x }/ m');
ylabel('{\it y }/ m');

%% Draw with bar elements

eldraw2(exb,eyb,[1 2 1])


