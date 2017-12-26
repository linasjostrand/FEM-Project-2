
% Load geometry and material parameters
run data.m

load_step = 100;
n_max = 10e3/load_step;
i_max = 50;

u = zeros(ndof,n_max+1);
df = zeros(ndof,1);
df(load_dof) = -load_step; % just dof D
f_ext = zeros(ndof,n_max+1);

% Convergence tolerance
tol = 1e-5;

for n=2:n_max+1
    f_ext(:,n) = f_ext(:,n-1) + df; 
    u(:,n) = u(:,n-1); 
    
    % Newton-Raphson algorithm
    for i=1:i_max
        
        % Stiffness matrix and out of balance forces for 3-node elements
        K = zeros(ndof,ndof);
        G = -f_ext(:,n);
        for el = 1:nelm
            ec = [ex(el,:); ey(el,:)];
            ed = u(edof(el,2:7),n);
            [ee,eff] = plan3gs(ec,ed);
            es = stresscal(ep,eff); % stress             
            D = mstiff(ep,eff); % material stiffness
            % Assemble stiffness matrix
            Ke = plan3ge(ec,t,D,ed,es);
            K(edof(el,2:7),edof(el,2:7)) = K(edof(el,2:7),edof(el,2:7)) + Ke;
            % Out of balance forces
            ef = plan3gf(ec,t,ed,es); 
            G(edof(el,2:7)) = G(edof(el,2:7)) + ef;
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
     
        % Solve for du
        du = solveq(K,-G,bc);
        % Update displacement u
        u(:,n) = u(:,n) + du;
        % Check convergence
        if norm(G(free_nodes)) < tol*norm(f_ext(:,n))
            fprintf('it %g conv in i = %g steps\n',n,i)
            break
        end
        

    end
    
end

%% Save deformations

save('def_data.mat','u')

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

for n = 1:5:n_max + 1
    
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
