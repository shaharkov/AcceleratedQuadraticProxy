%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script visualizes the results of the examples provided in this
%%% code package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% init stuff
solver_colors = {'b','m','c','r','k','g','y'};
switch size(F,2)
    case 3
        visF = F;
    case 4
        visF = boundaryFaces(F);
end


%% print summary
for ii = 1:n_solvers
    fprintf('%s - %g sec\n', logs{ii}.tag, sum(logs{ii}.t_iter));
end


%% visualize init state
figure(100); clf;
clear hP;
for ii = 1:n_solvers
    hP{ii} = patch('vertices', optimProblem.x0, 'faces', visF, 'FaceColor', solver_colors{ii}, 'FaceAlpha', 0.4);
end
axis equal;
legend(cellfun(@(x) x.tag, logs, 'UniformOutput', false));
if optimProblem.dim==3
    cameratoolbar;
    cameratoolbar('SetCoordSys','y');
end


%% go through iterations
figure(100);
if visualizeIterations
    iterVec = 1:max(cellfun(@(x) length(x.iter), logs));
else
    iterVec = num_iter;
end
for iter = iterVec
    % update solvers
    for ii = 1:n_solvers
        try
            hP{ii}.Vertices = logs{ii}.X(:,:,min(iter, end));
        catch
            warning('solver %d failed to draw', ii);
        end
    end
    title(iter);
    drawnow;
    pause(0.001);
end


%% plot some other information
figure(101); clf;
for ii = 1:n_solvers
    semilogy(logs{ii}.f, solver_colors{ii});
    hold on;
end
ylabel('functional value');
xlabel('iteration #');
legend(cellfun(@(x) x.tag, logs, 'UniformOutput', false));

figure(102); clf;
for ii = 1:n_solvers
    semilogy(logs{ii}.t, solver_colors{ii});
    hold on;
end
ylabel('step size (t)');
xlabel('iteration #');
legend(cellfun(@(x) x.tag, logs, 'UniformOutput', false));

figure(103); clf;
for ii = 1:n_solvers
    semilogy(cumsum([0; colStack(logs{ii}.t_iter(2:end))]), logs{ii}.f(1:end), solver_colors{ii});
    hold on;
end
ylabel('functional value');
xlabel('time [sec]');
legend(cellfun(@(x) x.tag, logs, 'UniformOutput', false));

figure(104); clf;
for ii = 1:n_solvers
    plot( abs(diff(logs{ii}.f))./(1+abs(logs{ii}.f(1:end-1))), solver_colors{ii});
    hold on;
end
ylabel('relative functional error')
xlabel('iteration #');
legend(cellfun(@(x) x.tag, logs, 'UniformOutput', false));


%% show per element energies (2d target domains only)
if optimProblem.dim==2
    % compute per element energies
    clear e;
    for ii = 1:n_solvers
        try
            e{ii} = optimProblem.evaluatePerElementEnergy(solver{ii}.x);
        catch
            warning('solver %d failed to draw', ii);
        end
    end
    c_min = min(cellfun(@(x) min(x), e));
    c_max = min(cellfun(@(x) max(x), e));
    
    % compute boundary
    temp.TR = triangulation(optimProblem.F, optimProblem.x0);
    temp.indB = temp.TR.freeBoundary();
    temp.indB = [temp.indB(:,1); temp.indB(1,1)];
    if ~exist('boundary_line_width','var')
        boundary_line_width = 1;
    end
    
    for ii = 1:n_solvers
        try
            figure(400+ii); clf;
            patch('faces',optimProblem.F,'vertices',logs{ii}.X(:,:,end),'FaceVertexCData',e{ii},'facecolor','flat','EdgeAlpha',0.03);
            title(sprintf('%s\nmean energy: %g\nmedian energy: %g\nmax energy: %g',logs{ii}.tag, mean(e{ii}), median(e{ii}), max(e{ii})));
            axis off;
            axis equal;
            colorbar;
            colormap(getDistortionColormap());
            caxis([c_min c_max]);
            % boundary
            line(logs{ii}.X(temp.indB,1,end), logs{ii}.X(temp.indB,2,end),...
                'color', solver_colors{ii}, 'linewidth', boundary_line_width);
        catch
            warning('solver %d failed to draw', ii);
        end
    end
end