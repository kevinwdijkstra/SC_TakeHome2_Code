D2pList = 2:7; % do not go above 7 for MATLAB implementation of 5 for user implementation
D3pList = 2:5; % do not go above 5 for MATLAB implementation of 3 for user implementation

epsilon = 1e-10; % convergence criteria IC BIM

use_MATLAB = true; % when true use matlab own implementation for code speed up

M = 2000; % max number of iterations to plot convergence



%% some setup stuff
D2nList = 2.^D2pList;
D3nList = 2.^D3pList;