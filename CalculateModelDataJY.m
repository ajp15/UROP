%Calculates given a set of parameters, the model results
 
function modelData = CalculateModelDataJY(model,parameters,constants,timeWindow,state)
 
    tStart=tic; 
    
    if nargin < 4
        timeWindow=[0,6];
    end
 
    if nargin < 5
        state = constants.ICs;
    end
 
    timeStart=timeWindow(1);
    timeEnd=timeWindow(2);
 
    if(isempty(state)==1)
        error 'state is empty'
    end
 
    tolerances = 1e-1;%change to -3
    options = odeset('Jacobian',@StateJac,'RelTol',tolerances,'AbsTol',[tolerances]);%add jacobian, 
    %instead of above, use the below. but adjust accordingly. need to
    %derive our own StateJac
    %options_state = odeset('Jacobian',@(t,x)StateJac(t,x,par,HybridP),'RelTol',10^-2);
    
    %StateJac function looks like this:
%   function jac = StateJac(t,x,par,HybridP)
%         jac = [ -(par(5) + HybridP(1)*par(4)) , -(par(1)*par(2)*par(3))/(par(3)*x(2) + 1)^2 
%         -par(8)*HybridP(2)*x(2)       , -(par(6)/(HybridP(1)*par(7) + 1) - par(8)*(HybridP(3) - HybridP(2)*x(1)))];
% 
%   end

%   have to write 4 Jacobian matrices for 4 different switch combination?
%   R on N on, R off N off, R on N off, R off N on

tspan=[timeStart;timeEnd];

    if(strfind(model,'FungalMacrophageModel')>0)
        [T,Y] = ode15s(@FMMOdesJY,tspan,state,options,parameters); %NEED TO CHANGE THIS BACK TO O.01 &15s!!
        %disp('ODE-ing..');
    else
        error 'incorrectly  defined model'
    end
 
    modelData=[T,Y];
   
    tElapsed_CalcModelData=toc(tStart)
end

