%Optimisation function - optimises parameters for a given model
 
%Applies a customised Monte Carlo moves algorithm with simulated annealing
%so a set of parameters and an objective/cost function
 
function OptimisationOutput = OptimiseModelJY(model,initialParameters,constants)    
 
    cost = CostFunctionJY(model,initialParameters,constants);
    cost_old = cost;
    mincost = cost;
    parameters = initialParameters;
    minParameters = parameters;
    lengthParameters=length(parameters);
    startk=constants.ignoreFirstNPar + 1;
 
%     fname = sprintf('minimised_Parameters_%s.txt', model);
%     fid = fopen(fname, 'a');
%     
    %simulated Annealing
    progressbar(0)
    for j=1:constants.annealingSteps
          
          try
          %T=0.298-(0.298/m)*j;
          T=1/j;
 
          %standard Monte Carlo moves algorithm
          for i=1:constants.optSteps
            %maximum displacement is restricted by Delta_max.
            Delta_max=0.137*sqrt(T);
            for k=1:lengthParameters
              oldParameters(k)=parameters(k);
            end            
            
            %create a new point - ignore first startk-1 paramters
            for k=startk:lengthParameters
              parameters(k)=parameters(k)+(rand-0.5)*Delta_max;
            end
            %Assuming the first startk-1 parameters are not necessarily
            %positive
            %Ensure others are bigger than zero
            for k=startk:lengthParameters
              while parameters(k)<0
                parameters(k)=oldParameters(k)+(rand-0.5)*Delta_max;
              end
            end
           
            %ensure we do not change specific parameters
            for k=constants.ignorePar
                parameters(k)=oldParameters(k);
            end
 
            a = (parameters(4)*parameters(7)*parameters(5))/(parameters(2)*parameters(6));
            b = 1+ (parameters(4)/parameters(2));
            c = - 1;
            NFs0   = ( -b + sqrt( b^2 - (4*a*c)))*(1/(2*a));
            RCAN10 = NFs0*parameters(5)/parameters(6); 
            parameters(19)=NFs0*(1.001 + (1.4-1.001)*rand);
            parameters(14)=RCAN10*(1.001 + (1.4-1.001)*rand);
            %calculate the cost function for that point
            cost = CostFunctionJY(model,parameters,constants);
            cost_new=cost;
            
            % Monte Carlo criteria for accepting the new point
            if (cost_new < cost_old)
              %accept move
              disp('accept move');
              disp(parameters);
              cost_old = cost_new;
     
            elseif (rand < exp(-(cost_new-cost_old)/T))
              %accept move
              cost_old = cost_new;
              disp('accept move');
              disp(parameters);
            else
              %reject move
              for k=1:lengthParameters
                parameters(k)=oldParameters(k);
              end
            end
          end
          catch
              continue
          end
          
          progressbar(j/constants.annealingSteps);
    end
 
    cost = CostFunctionJY(model,parameters,constants);
    mincost=cost;
    minParameters=parameters
    minParameters
    
    
    %What is the absolute minimum from all runs? Return this
    OptimisationOutput.mincost = mincost
    OptimisationOutput.parameters = minParameters
    
%     %Write the data to a file
%         fprintf(fid, '\n%f', mincost);
%         for i=1:length(minParameters)
%             fprintf(fid, ',%f', minParameters(i));
%             showminParam(i)=minParameters(i);
%         end
%  
%     fclose(fid);
 
end
