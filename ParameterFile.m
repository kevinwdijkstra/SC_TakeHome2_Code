D2pList = 2:7;          % do not go above 7 for MATLAB implementation or 5 for user implementation
D3pList = 2:5;          % do not go above 5 for MATLAB implementation or 3 for user implementation

epsilon = 1e-10;        % convergence criteria IC BIM

use_MATLAB  = true;     % when true use matlab own implementation for code speed up
plot_figure = true;     % plot figures
use_symrcm  = true;     % use matrix reordering


M = 1e4;                % max number of iterations to plot convergence for preallocation of memory

do_DIRECT   = false;     % execute direct solvers
do_IC_BIM   = true;     % execute IC BIM solvers
do_ICCG     = false;     % execute ICCG solvers

%% some setup stuff
D2nList = 2.^D2pList;
D3nList = 2.^D3pList;



%% 
solve_options.use_MATLAB    = use_MATLAB;
solve_options.plot_figure   = plot_figure;
solve_options.use_symrcm    = use_symrcm;
solve_options.epsilon       = epsilon;
solve_options.M             = M;
