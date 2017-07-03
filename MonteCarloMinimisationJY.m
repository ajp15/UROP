%Want to optimise a set of parameters in a small region around this one. NOTE this optimisation is not very good for a spiky landscape.
function minParameters = MonteCarloMinimisationJY(model,initialParameters)
 
 
%need to run this command for matlab to work in parallel.
%matlabpool
 
%time how long the parallelised code takes to run
tStart=tic;
 
%parfor i=1:4
 
    %%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%
    %model='FungalMacrophageModel';
 
    %optimum parameter set
    %initialParameters = []
 
    %calculate the constants needed for this model
    constants = CalculateConstantsJY();
    setting='10pcVar';
    delta = .1; %change this to alter how fair it 'jumps', ori=0.1
    %%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%
 
 
%     fname = sprintf('MCminimised_Parameters_%s_%s.txt', model,setting);
%     fid = fopen(fname, 'a');
    
    mincost = CostFunctionJY(model,initialParameters,constants);
    minParameters = initialParameters;
 
    progressbar(0)
    for k = 1:constants.numRanParameters
 
        parameters=minParameters
 
        for j = (constants.ignoreFirstNPar+1):length(minParameters)
 
            %pick a random parameter point variation within a hypercube
            parameters(j) = parameters(j) + minParameters(j)*delta*(rand-.5);
 
            while(parameters(j)<0)
                parameters(j) = parameters(j) + minParameters(j)*delta*2*(rand-.5);
            end
 
        end
 
        for j = constants.ignorePar
            parameters(j) = minParameters(j);
        end
 
        cost = CostFunctionJY(model,parameters,constants);
 
        if ((cost<mincost)&(imag(cost)==0))
            minParameters = parameters
            mincost=cost;
        end
        
        progressbar(k/constants.numRanParameters);
 
    end
 
        %Write the data to a file
%         fprintf(fid, '\n%f', mincost);
%         for i=1:length(minParameters)
%             fprintf(fid, ',%f', minParameters(i));
%             showminParam(i)=minParameters(i);
%         end
%    fclose(fid);
 
%end
 
 
 
tElapsed=toc(tStart)
 
end
 
 
