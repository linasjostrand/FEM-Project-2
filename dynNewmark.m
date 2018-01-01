
% Load geometry and material parameters
run data.m

n_max = 2000;
i_max = 50;
beta = 1/3;
gamma = 1/2;
h = 0.00001; % time increment 1 ms?

c1 = 1/(beta*h^2); c2 = 1/(beta*h); c3 = (1-2*beta)/(2*beta);
c4 = gamma/(beta*h); c5 = (gamma-beta)/beta; c6 = h*(gamma-2*beta)/(2*beta);

load def_data.mat
u_def = u(:,end);

% Initialize displacements
u = zeros(ndof,n_max + 1);
u(:,1) = u_def;
% Initialize velocities
vel = zeros(ndof,n_max + 1); % 0?
% Initialize force
f_ext = zeros(ndof,n_max + 1);
% f_ext(load_dof,1) = -10e3;
f_ext(load_dof,1:101) = -10e3:10e1:0;
% Initialize accelerations
acc = zeros(ndof,n_max + 1);

% Convergence tolerance
tol = 1e-6;

for n=2:n_max+1
    
    % Star?
    acc_star = c1*u(:,n-1) + c2*vel(:,n-1) + c3*acc(:,n-1);
    vel_star = c4*u(:,n-1) + c5*vel(:,n-1) + c6*acc(:,n-1);
    
    % Predictor (acceleration same)
    acc(:,n) = acc(:,n-1);
    vel(:,n) = vel(:,n-1) + h*acc(:,n-1);
    u(:,n) = u(:,n-1) + h*vel(:,n-1) + h^2/2*acc(:,n-1);
    
    % Newton-Raphson algorithm
    for i=1:i_max
        
        % Stiffness matrix and out of balance forces for 3-node elements
        K = zeros(ndof,ndof);
        M = zeros(ndof,ndof);
        G = -f_ext(:,n);
        for el = 1:nelm
            ec = [ex(el,:); ey(el,:)];
            ed = u(edof(el,2:7),n);
            [ee,eff] = plan3gs(ec,ed);
            es = stresscal(ep,eff); % stress             
            D = mstiff(ep,eff); % material stiffness
            % Assemble stiffness matrix
            Ke = plan3ge(ec,t,D,ed,es);
            Me = plan3gm(ec,t,rho);
            M(edof(el,2:7),edof(el,2:7)) = M(edof(el,2:7),edof(el,2:7)) + Me;
            K(edof(el,2:7),edof(el,2:7)) = K(edof(el,2:7),edof(el,2:7)) + Ke;
            % Out of balance forces
            ef = plan3gf(ec,t,ed,es); 
            G(edof(el,2:7)) = G(edof(el,2:7)) + ef + c1*Me*u(edof(el,2:7),n) - Me*acc_star(edof(el,2:7));%- Me*acc(edof(el,2:7),n);%+ 
        end
        % Bars
        for el = 1:nelmb
            ec = ecb{el};
            ed = u(edofb(el,2:7),n)';
            [~,ee] = bar3gs(ec,ep,ed);
            es = norfb(ec,ee,k,r);
            D_bar = bstiff(ec,ee,k,r);
            % Assemble stiffness matrix
            Ke = bar3ge(ec,D_bar,ed,es);
            K(edofb(el,2:7),edofb(el,2:7)) = K(edofb(el,2:7),edofb(el,2:7)) + Ke;
            % Out of balance forces
            ef = bar3gf(ec,ed,es); 
            G(edofb(el,2:7)) = G(edofb(el,2:7)) + ef;
        end
        K_star = K + 1/(beta*h^2)*M;

        % Solve for du
        du = solveq(K_star,-G,bc);
        
        % Update displacement u
        u(:,n) = u(:,n) + du;
        vel(:,n) = vel(:,n) + c4*du;
        acc(:,n) = acc(:,n) + c1*du;
        
%         fprintf('norm G = %g\n',norm(G(free_nodes)))

        % Check convergence
        if norm(G(free_nodes)) < tol
            fprintf('it %g conv in i = %g steps\n',n,i)
            break
        end
        

    end
    
end

%% Plot deformed geometry
eldisp2(ex,ey,extract(edof,u(:,end)),[1 4 1],1);
hold on
eldisp2(exb,eyb,extract(edofb,u(:,end)),[1 2 1],1);
viscircles([xc1 yc1; xc2 yc2],[r;r]);
% ylim([-0.2 yc2+r]);
title(['Deformed geometry, k = ' num2str(k,'%10.1e')],'Fontsize',13);
xlabel('{\it x }/ m');
ylabel('{\it y }/ m');

%% Plot deformation of the geometry

figure

for n = 1:50:n_max + 1
    
    clf
    eldisp2(ex,ey,extract(edof,u(:,n)),[1 4 1],1);
    hold on
    eldisp2(exb,eyb,extract(edofb,u(:,n)),[1 2 1],1);
    viscircles([xc1 yc1; xc2 yc2],[r;r]);
    grid on
    ylim([-0.2 yc2+r]);
    pause

end

%% Plot the force vs displacement 

figure
plot(-u(load_dof,:),-f_ext(load_dof,:))
hold on
xlabel('-u','Fontsize',13); ylabel('-f','Fontsize',13)
