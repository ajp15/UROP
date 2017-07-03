function [modelData] = HybridSystem_RecursiveCodeJY(parameters,constants,timeWindow,state)
%state here is the initial conditions
    tStart=tic; 

    if nargin < 3
        timeWindow=[0,6];
    end
 
    if nargin < 4
        state = constants.ICs;
    end

    timeStart=timeWindow(1);
    timeEnd=timeWindow(2);
    
    if(isempty(state)==1)
        error 'state is empty' %check if initial condition is 0
    end
    
    %base cases
    time=timeStart;
    IC_ODE=state;
    n=0; %to switch between different cases
    counter=0;
    
    tolerances = 1e-8;
    
    
    %different option accounts for different event when switch is on or off
    %can change to different jacobian for different switch combination?
    options_ODE_Low = odeset('Events',@event_function_Hybrid_Low,'Jacobian',@StateJacRoff,'RelTol',tolerances,'AbsTol',[tolerances]);
    % ^when R switch OFF, want the simulation to stop when hits threshold
    % with +ve gradient
    options_ODE_High = odeset('Events',@event_function_Hybrid_High,'Jacobian',@StateJac,'RelTol',tolerances,'AbsTol',[tolerances]);
       
    
    while time ~= timeEnd
        switch n
            case 0 %LOW = Initial case
                solutions_ODE_Low = ode15s(@FMMOdesLow,[time;timeEnd],IC_ODE,options_ODE_Low,parameters);
                %solver returns column vector
                tint_ODE = linspace(time,solutions_ODE_Low.x(end),10000); %row vector of size 1 x 10k
                Solutions_Low = deval(solutions_ODE_Low,tint_ODE); %returns row vector (6 rows, 10k cols)
               
                %make new matrix containing the solution
                if counter==0
                    Sol_Assembly=[tint_ODE; Solutions_Low];
                else
                    Sol_Assembly=[Sol_Assembly [tint_ODE; Solutions_Low]]; %size 7x10k
                end 
                
                time=tint_ODE(end);
                IC_ODE = Solutions_Low(:,end);
                
                if time ~= timeEnd
                    n = 1; %move on to the next case, R switch ON
                    counter=counter+1;
                    %disp('n->1');
                end
                
             case 1
                solutions_ODE_High = ode15s(@FMMOdesHigh,[time;timeEnd],IC_ODE,options_ODE_High,parameters);
                tint_ODE_high = linspace(time,solutions_ODE_High.x(end),10000);
                Solutions_High = deval(solutions_ODE_High,tint_ODE_high);
              
                if Solutions_High(3,end-1)<parameters(14)
                    n=0;
                    disp('n->0')
                else
                    Sol_Assembly=[Sol_Assembly [tint_ODE_high;Solutions_High]];
                    size_solAssembly=size(Sol_Assembly);
                    time=tint_ODE_high(end);


                    if time ~= timeEnd
                        IC_ODE = Solutions_High(:,end);
                        IC_ODE_case1=size(IC_ODE); 
                        n = 0; %move on to the next case, R switch Off
                        counter=counter+1;
                        %disp('n->0')
                    end
                end
        end
                
    end
    
    modelData=transpose(Sol_Assembly); %because in printfittingres, every variable is col vector
    initialmodelData=modelData(1,:)
    
    tElapsed_HybridSys=toc(tStart)

end

