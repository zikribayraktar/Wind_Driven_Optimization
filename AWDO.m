function AWDO_v03()
%-------------------------------------------------------------------------
% Sample Matlab Code for the Adaptive Wind Driven Optimization.
% The coefficients of the WDO are optimized by the CMEAS algorithm. 
% No additional cost function calls are needed.
%
% CHANGES:
% v02: Removed SORTING in WDO, because it messes up CMAES indexing.
% v03: Included more comments.
%
% Contact Information:
% Zikri Bayraktar, Ph.D.
% Email: thewdoalgorithm@gmail.com
%
% DISCLAIMER: This sample code is provided for educational purposes only. 
% Please also read the referred papers to better understand the variables
% as well as the inner-workings of the algorithms.
% Use at own your risk! There is no guarantee that the code is bug free.
%-------------------------------------------------------------------------
%
% WDO REFERENCE PAPERS:
% Please refer to the following TWO articles in your publications if you use this algorithm:
%
% 1.) Z. Bayraktar, and M. Komurcu, "Adaptive Wind Driven Optimization," 
% Proceedings of the 9th EAI International Conference on Bio-inspired Information and 
% Communications Technologies (Formerly BIONETICS) on 9th EAI International Conference 
% on Bio-inspired Information and Communications Technologies (Formerly BIONETICS), 
% New York City, NY, Dec. 3-5, 2015.
% http://dl.acm.org/citation.cfm?id=2954811&CFID=806983354&CFTOKEN=67304345
%
% 2.) Z. Bayraktar, M. Komurcu, J. A. Bossard and D. H. Werner, "The Wind 
% Driven Optimization Technique and its Application in Electromagnetics," 
% IEEE Transactions on Antennas and Propagation, Volume 61, Issue 5, 
% pages 2745 - 2757, May 2013.
% http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6407788&newsearch=true&queryText=wind%20driven%20optimization
%
%-------------------------------------------------------------------------

tic; clear; close all; clc;  
format long g;
delete('WDOoutput.txt');  
delete('WDOpressure.txt');  
delete('WDOposition.txt');
fid=fopen('WDOoutput.txt','a');
%--------------------------------------------------------------

% User defined parameters:
param.popsize = 20;     % population size.
param.npar = 10;		% Dimension of the problem.
param.maxit = 1000;		% Maximum number of iterations.
maxV = 0.3;             % maximum allowed speed.
dimMin = -100;			% Lower dimension boundary.
dimMax = 100;			% Upper dimension boundary.

% Unlike the classical WDO, Adaptive WDO does not need the coefficients to
% be predetermined. Coefficients; alpha, RT, g, and c, will be selected by
% the CMAES algorithm:
rec.arx = rand(4,param.popsize);   %consistent with the CMAES indexing
%---------------------------------------------------------------

% Initialize WDO population, position and velocity:
% Randomize population in the range of [-1, 1]:
% Please note that WDO/AWDO always works in the range of [-1, 1], but is mapped
% to the upper/lower bounds of the problem before calling the pressure (cost,fitness) function.

pos = 2*(rand(param.popsize,param.npar)-0.5);
% Randomize velocity:
vel = maxV * 2 * (rand(param.popsize,param.npar)-0.5);  
%---------------------------------------------------------------

% Call the pressure (cost,fitness) function
% Evaluate initial population one member at a time: (simple Sphere Function is utilized here)
for K=1:param.popsize,
	% map the AWDO 'pos' vector to problem itself using the upper/lower boundaries of each dimension.
	x = (dimMax - dimMin) * ((pos(K,:)+1)./2) + dimMin;
	% call the pressure (cost,fitness) function:
	pres(K,:) = fsphere(x);   %here insert your own cost function and make sure the mapping above carried out properly.
end
%----------------------------------------------------------------

% Finding best air parcel in the initial population :
[globalpres,indx] = min(pres);
globalpos = pos(indx,:);
minpres(1) = min(pres);			% minimum pressure
%-----------------------------------------------------------------

% Rank the air parcels:
[sorted_pres rank_ind] = sort(pres);
% Do not sort the air parcels! Sorting mixes up CMAES indexing.
% pos = pos(rank_ind,:);
keepglob(1) = globalpres;
%-----------------------------------------------------------------

% Start iterations :
iter = 1;   % iteration counter
for ij = 2:param.maxit,
    	% Update the velocity:
    	for i=1:param.popsize
		% choose random dimensions:
		a = randperm(param.npar);        			
		% choose velocity based on random dimension:
    		velot(i,:) = vel(i,a);				
        	vel(i,:) = (1-rec.arx(1,i))*vel(i,:)-(rec.arx(2,i)*pos(i,:))+ ...
				    abs(1-1/rank_ind(i))*((globalpos-pos(i,:)).*rec.arx(3,i))+ ...
				    (rec.arx(4,i)*velot(i,:)/rank_ind(i));
        end
    
        % Check velocity:
        vel = min(vel, maxV);
        vel = max(vel, -maxV);
        % Update air parcel positions:
    	pos = pos + vel;
        pos = min(pos, 1.0);
        pos = max(pos, -1.0); 
			
		% Call the pressure (cost,fitness) function
		% Evaluate initial population one member at a time: (simple Sphere Function is utilized here)
		for K=1:param.popsize,	
			% map the AWDO 'pos' vector to problem itself using the upper/lower boundaries of each dimension.
			x = (dimMax - dimMin) * ((pos(K,:)+1)./2) + dimMin;
			% call the pressure (cost,fitness) function:
			pres(K,:) = fsphere(x);   %here insert your own cost function and make sure the mapping above carried out properly.
		end		

        % call CMAES with the coefficients and corresponding pressure, so
        % that CMAES can return the new set of coefficient values for next
        % iteration:
        [rec] = purecmaes_wdo(ij,rec,param.popsize,pres);

    	%----------------------------------------------------
    	% Finding best particle in population
    	[minpres,indx] = min(pres);
    	minpos = pos(indx,:);          % min location for this iteration
    	%----------------------------------------------------
    	% Rank the air parcels:
    	[sorted_pres rank_ind] = sort(pres);
    	% Do not sort the air parcels position, velocity and pressure!
    	% Instead use "rank_ind" if needed.
    	% pos = pos(rank_ind,:);
    	% vel = vel(rank_ind,:);
    	% pres = sorted_pres;  
    
    	% Updating the global best:
    	better = minpres < globalpres;
    	if better
        		globalpres = minpres		% global minimum pressure (cost,fitness) value
        		globalpos = minpos;			% global minimum position vector, note that it is in the range of [-1, 1].
		end
		% Keep a record of the progress:
    	keepglob(ij) = globalpres;
%     	save WDOposition.txt pos -ascii -tabs;
end
		%Save values to the final file.
    	pressure = transpose(keepglob);
        filenamestr = ['WDOpressure.txt'];
    	save(filenamestr, 'pressure', '-ascii' , '-tabs');
		
		% note that the 'globalpos' is the best solution vector found.
		% 'globalpos' is in the range of [-1 1] and needs to be scaled with upper/lower bounds of the specific problem that you are trying to solve.

end
% end-of-WDO.
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

function [rec] = purecmaes_wdo(counteval,rec,npop,pressure)
%   coeff -- CMAES optimized coefficients returned to WDO.
%   counteval -- Iteration counter.
%   rec -- Record of prior values used in CMAES.
%   npop -- number of population members from WDO, each member gets their
%          own set of coefficients determined by the CMAES.
%   pressure -- pressure(cost function) computed by WDO for the set of
%               coefficients that CMEAS picked last iteration
%---------------------------------------------------------
%   This script is modified version of the publicly available purecmaes.m
%   Original "purecmaes.m" can be accessed from:
%   URL: http://www.lri.fr/~hansen/purecmaes.m
%
%   This function takes in pressure values from the WDO and computes new
%   set of coefficients for the WDO and returns them. In simpler terms, it
%   optimizes the WDO coefficients of "alpha, g, RT, c"
%
%   This function is not optimized, only proof-of-concept code to
%   illustrate that the WDO coefficients can be optimized by CMAES, 
%   which in turn creates and adaptive WDO algorithm.
%---------------------------------------------------------

if counteval==2    %initialization only happens when the CMAES called for the first time.
    rec.N = 4;  % problem dimension is fixed to 4: alp, g, RT, c 
    rec.xmean = rand(rec.N,1);  % objective variables initial point
    rec.sigma = 0.5;        % coordinate wise standard deviation (step size)
  
    % Strategy parameter setting: Selection  
%     lambda = 4+floor(3*log(N));  % population size, offspring number
    rec.lambda = npop; %this is defined by the wdo population size.
    rec.mu = rec.lambda/2;               % number of parents/points for recombination
    rec.weights = log(rec.mu+1/2)-log(1:rec.mu)'; % muXone array for weighted recombination
    rec.mu = floor(rec.mu);        
    rec.weights = rec.weights/sum(rec.weights);     % normalize recombination weights array
    rec.mueff=sum(rec.weights)^2/sum(rec.weights.^2); % variance-effectiveness of sum w_i x_i

    % Strategy parameter setting: Adaptation
    rec.cc = (4 + rec.mueff/rec.N) / (rec.N+4 + 2*rec.mueff/rec.N); % time constant for cumulation for C
    rec.cs = (rec.mueff+2) / (rec.N+rec.mueff+5);  % t-const for cumulation for sigma control
    rec.c1 = 2 / ((rec.N+1.3)^2+rec.mueff);    % learning rate for rank-one update of C
    rec.cmu = min(1-rec.c1, 2 * (rec.mueff-2+1/rec.mueff) / ((rec.N+2)^2+rec.mueff));  % and for rank-mu update
    rec.damps = 1 + 2*max(0, sqrt((rec.mueff-1)/(rec.N+1))-1) + rec.cs; % damping for sigma 
                                                      % usually close to 1
    % Initialize dynamic (internal) strategy parameters and constants
    rec.pc = zeros(rec.N,1); 
    rec.ps = zeros(rec.N,1);   % pc, ps: evolution paths for C and sigma
    rec.B = eye(rec.N,rec.N);                       % B defines the coordinate system
    rec.D = ones(rec.N,1);                      % diagonal D defines the scaling
    rec.C = rec.B * diag(rec.D.^2) * rec.B';            % covariance matrix C
    rec.invsqrtC = rec.B * diag(rec.D.^-1) * rec.B';    % C^-1/2 
    rec.eigeneval = 0;                      % track update of B and D
    rec.chiN=rec.N^0.5*(1-1/(4*rec.N)+1/(21*rec.N^2));  % expectation of  ||N(0,I)|| == norm(randn(N,1))
end

    % get the fitness from WDO pressure:
    rec.arfitness = pressure';
    
    % Sort by fitness and compute weighted mean into xmean
    [rec.arfitness, rec.arindex] = sort(rec.arfitness);  % minimization
    rec.xold = rec.xmean;
    rec.xmean = rec.arx(:,rec.arindex(1:rec.mu)) * rec.weights;  % recombination, new mean value
    
    % Cumulation: Update evolution paths
    rec.ps = (1-rec.cs) * rec.ps ... 
          + sqrt(rec.cs*(2-rec.cs)*rec.mueff) * rec.invsqrtC * (rec.xmean-rec.xold) / rec.sigma; 
    rec.hsig = sum(rec.ps.^2)/(1-(1-rec.cs)^(2*counteval/rec.lambda))/rec.N < 2 + 4/(rec.N+1);
    rec.pc = (1-rec.cc) * rec.pc ...
          + rec.hsig * sqrt(rec.cc*(2-rec.cc)*rec.mueff) * (rec.xmean-rec.xold) / rec.sigma; 

    % Adapt covariance matrix C
    rec.artmp = (1/rec.sigma) * (rec.arx(:,rec.arindex(1:rec.mu)) - repmat(rec.xold,1,rec.mu));  % mu difference vectors
    rec.C = (1-rec.c1-rec.cmu) * rec.C ...                   % regard old matrix  
         + rec.c1 * (rec.pc * rec.pc' ...                % plus rank one update
                 + (1-rec.hsig) * rec.cc*(2-rec.cc) * rec.C) ... % minor correction if hsig==0
         + rec.cmu * rec.artmp * diag(rec.weights) * rec.artmp'; % plus rank mu update 

    % Adapt step size sigma
    rec.sigma = rec.sigma * exp((rec.cs/rec.damps)*(norm(rec.ps)/rec.chiN - 1)); 
    
    % Update B and D from C
    if counteval - rec.eigeneval > rec.lambda/(rec.c1+rec.cmu)/rec.N/10  % to achieve O(N^2)
      rec.eigeneval = counteval;
      rec.C = triu(rec.C) + triu(rec.C,1)'; % enforce symmetry
      [rec.B,rec.D] = eig(rec.C);           % eigen decomposition, B==normalized eigenvectors
      rec.D = sqrt(diag(rec.D));        % D contains standard deviations now
      rec.invsqrtC = rec.B * diag(rec.D.^-1) * rec.B';
    end
    
    % Generate and evaluate lambda offspring
    for k=1:rec.lambda,
        rec.arx(:,k) = rec.xmean + rec.sigma * rec.B * (rec.D .* randn(rec.N,1)); % m + sig * Normal(0,C) 
    end

end
% end-of-CMAES.
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------



%%% TEST FUNCTIONS:
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function f=fsphere(x)
  f=sum(x.^2);
end
% -----------------------------


%----------------------------------------------------------------------
%----------------------------------------------------------------------
%-----------------------||| CMAES REFERENCES |||-----------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%
% Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% Strategy on Multimodal Test Functions.  Eighth International
% Conference on Parallel Problem Solving from Nature PPSN VIII,
% Proceedings, pp. 282-291, Berlin: Springer. 
% (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
% 
% Further references:
% Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% 9(2), pp. 159-195.
% (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
%
% Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% Time Complexity of the Derandomized Evolution Strategy with
% Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
%-----------------------------------------------------
% end-of-file