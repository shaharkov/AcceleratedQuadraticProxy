%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initialize

% reset timer;
tic; dummy = toc;

% set paths
if isempty(whos('global','path_def')),
    fprintf('- Adding toolbox paths\n');
    
    % common path
    addpath(genpath('toolbox'));
    addpath(genpath('code'));
    addpath(genpath('mex'));
    
    if isunix,
        opengl neverselect % disable opengl
        set(0, 'DefaultFigureRenderer', 'zbuffer');
    end
    
    global path_def
end;
