function rayleigh_modes()
% =========================================================================
%  MULTI-MODE RAYLEIGH–STODOLA TOOL WITH MASS ORTHOGONALIZATION
%
%  - Computes first N modes
%  - Orthogonalizes each higher mode against all previous modes
%  - Normalizes top DOF = 1
%  - FINAL PLOT: shows all converged mode shapes side-by-side
% =========================================================================

clc; close all;

disp('----------------------------------------------------------');
disp('   MULTI-MODE RAYLEIGH–STODOLA MODE SHAPE COMPUTATION');
disp('----------------------------------------------------------');

%% ===================== USER INPUT =======================

m = input('Enter story masses m (kg), row vector (base included):\nm = ');
k = input('\nEnter story stiffness k (N/m), row vector (base included):\nk = ');
num_modes = input('\nHow many modes do you want to extract?  num_modes = ');

m = m(:); 
k = k(:);
N = length(m);

if num_modes > N-1
    error('Number of modes cannot exceed N-1 DOF.');
end

%% ==================== STORAGE ============================
phi = zeros(N, num_modes);      % final mode shapes
omega = zeros(num_modes,1);     % natural frequencies

%% ==================== MODE LOOP ==========================

for mode = 1:num_modes

    fprintf('\n-------------------------------------------\n');
    fprintf('   COMPUTING MODE %d\n', mode);
    fprintf('-------------------------------------------\n');

    % Initial shape guess
    shape = linspace(0,1,N)';  

    % Initial estimate
    omega_old = angular_frequency(m, shape, k);

    iter = 1;
    tol = 1e-6;

    while true

        % Internal forces
        F_I = m .* (omega_old^2) .* shape;

        % Shear forces
        SF = flip(cumsum(flip(F_I)));

        % Story drifts
        story_dis = SF ./ k;
        story_dis(isinf(story_dis)) = 0;

        % Floor displacements
        tot_dis = cumsum(story_dis);

        % Updated (unnormalized) shape
        shape_new = tot_dis;

        % ------------------------------------------------------
        %        MASS ORTHOGONALIZATION
        % ------------------------------------------------------
        if mode > 1
            for j = 1:mode-1
                proj = (phi(:,j)' * diag(m) * shape_new) / ...
                       (phi(:,j)' * diag(m) * phi(:,j));
                shape_new = shape_new - proj * phi(:,j);
            end
        end

        % Normalize (top DOF = 1)
        shape_new = shape_new / shape_new(end);

        % New frequency
        omega_new = angular_frequency(m, shape_new, k);

        fprintf('Mode %d Iter %d:  ω = %.6f\n', mode, iter, omega_new);

        % Convergence check
        if abs(omega_new - omega_old) < tol
            break;
        end
        
        shape = shape_new;
        omega_old = omega_new;
        iter = iter + 1;
    end

    % Save final mode
    phi(:,mode) = shape_new;
    omega(mode) = omega_new;

end

%% ======================== OUTPUT =========================

fprintf('\nFINAL RESULTS\n');
for i = 1:num_modes
    fprintf('Mode %d:  ω = %.6f rad/s,   T = %.6f s\n', ...
        i, omega(i), 2*pi/omega(i));
end

%% ======================== FINAL MODE PLOT ==========================

figure('Position',[200 50 1200 700]);
floors = (0:N-1)';

for mode = 1:num_modes
    subplot(1, num_modes, mode);
    plot(phi(:,mode), floors, 'r-o','LineWidth',1.8);
    set(gca,'YDir','normal');
    xlabel('Shape Value');
    ylabel('Floor Level');
    title(sprintf('Mode %d Shape', mode));
    grid on;
end

sgtitle('Final Converged Mode Shapes (Mass-Orthogonalized)');

end


% ===========================================================
%   SUBFUNCTION: ANGULAR FREQUENCY FROM RAYLEIGH QUOTIENT
% ===========================================================
function omega = angular_frequency(m, shape_f, k)

    F_I = m .* shape_f;
    SF = flip(cumsum(flip(F_I)));
    story_dis = SF ./ k;
    story_dis(isinf(story_dis)) = 0;

    num = sum(k .* story_dis.^2);
    den = sum(m .* shape_f.^2);

    omega = sqrt(num / den);
end
