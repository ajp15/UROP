%function that controls the parameter optimisations for the fungal macrophage model
function ControlJY(initialParameters,optimisation)
 
    %time how long the parallelised code takes to run
    tStart=tic;
 
    if(nargin<2)
        %default to the MCSA optimisation
        optimisation = 1;
    end
 
    if(nargin<1)
        error 'input parameters!'
    end
 
    if(length(initialParameters)<19)
        error 'incorrect form of parameters'
    end
 
 
    %%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%
    %define the name of the model
    model='FungalMacrophageModel';
 
    %example parameter sets
 
 %   initialParameters(12)=0;%as there is no killing due to the inflammation in any of the experiments

    %calculate the constants needed for this model
    constants = CalculateConstantsJY();
    %%%%%%%%%%%%
 
 
    if(optimisation==1)
        %optimise using simulated Annealing
        OptimisationOutput = OptimiseModelJY(model,initialParameters,constants)
    else
        %optimise using MonteCarlo steps (faster)
        OptimisationOutput.parameters = MonteCarloMinimisationJY(model,initialParameters)
    end
 
    %OptimisationOutput.mincost%DEBUG: for comparison with below
    constants.cost = CostFunctionJY(model,OptimisationOutput.parameters,constants,1)
    tElapsed = toc(tStart)
 
end
