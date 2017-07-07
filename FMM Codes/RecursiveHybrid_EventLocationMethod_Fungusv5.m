%Solve model ODEs using the initial conditions passed into the function
function [modelData,Switch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5(t_ODE,IC_States,IC_Switches,parameters)
    % Switch_state = 2xN
    tStart=tic; 
    time = t_ODE(1); % base case
    timeEnd=t_ODE(2);
    SwitchStatus=IC_Switches;
    tau=parameters(15);

    if(isempty(IC_States)==1)
        error 'state is empty' %check if initial condition is 0
    end

    IC_ODE=IC_States;
    Rc=parameters(14);
    Nc=parameters(19);

    T = [];
    Y = [];
    Switch_state = [];
    SwitchIDhistory=[];
    timeThreshold=timeEnd;
    WhichEvent =[];


    while time ~= timeEnd
        EventFcn=@(t,y)HybridSwitchingEvent(t,y,Rc,Nc,timeThreshold,tau,SwitchStatus);
        %Tips: sometimes, you may adjust RelTol or AbsTol to get different
        %results. To smaller the value the higher the resolution
        options = odeset('Events',EventFcn,'Jacobian',@(t,y)StateJac(t,y,parameters,SwitchStatus),'RelTol',1e-8,'AbsTol',1e-8);
        solution  = ode15s(@(t,y)CompleteHybridSystem(t,y,parameters,SwitchStatus),t_ODE,IC_ODE,options);
        t_ODE_extended = linspace(time,solution.x(end),1000); % NOTE: 1000 is arbitrary increase or decrease depending on requirements
        Solution = deval(solution,t_ODE_extended);     
        T = [T,t_ODE_extended]; 
        Y = [Y,Solution]; 
        Switch_state = [Switch_state,repmat(SwitchStatus',1,numel(t_ODE_extended))];
        % Update recursive variable
        time = solution.x(end);
        if time ~= timeEnd
            WhichEvent=solution.ie(end);
            if WhichEvent==1  %hits threshold
                timeThreshold=solution.xe(end);
                t_ODE = [solution.xe(end),timeEnd];
                IC_ODE = Solution(:,end);
            
            elseif WhichEvent==3 || WhichEvent==2 %after time delay or N hits threshold, change the switch
                if WhichEvent==3 
                    timeThreshold = timeEnd;
                end
                SwitchStatus = ChangeSwitchState(WhichEvent,SwitchStatus);
                t_ODE = [solution.xe(end),timeEnd];
                IC_ODE = Solution(:,end);
                %display('........');
            end     
        end  
    end
    
    %Switchsize=size(Switch_state) %2 rows n columns
    modelData=[T;Y];
    modelData=transpose(modelData); %because in printfittingres, every variable is col vector
  
    tElapsed_HybridSys=toc(tStart);
end

function dydt = CompleteHybridSystem(t,y,parameters,SwitchStatus)
    a           = parameters(1);%fungal growth rate
    kN         = parameters(2);%NFAT activation rate
    yNK         = parameters(3);%FK506 inhibition rate (of NFAT activation)
    deltaN      = parameters(4);%NFAT deactivation rate
    kRN          = parameters(5);%RCAN1 production rate
    deltaR      = parameters(6);%RCAN1 degredation rate
    yNR         = parameters(7);%RCAN1 inhibition rate (of NFAT activation)
    kR         = parameters(8);%RCAN1 production rate
    betaNR         = parameters(9);%Enhanced NFAT activation rate (by RCAN1-p)
    deltaRP     = parameters(10);%RCAN1-p degredation rate
    kT          = parameters(11);%TNFalpha production
    deltaTF          = parameters(12);%Fungal killing rate
    deltaT      = parameters(13);%TNFalpha degredation rate
    Rc           = parameters(14);%RCAN1 switch threshold
    x           = parameters(15);%RCAN1 switch hill coefficient
    betaNF         = parameters(16);%NFAT activation enhancement by F
    betaRF         = parameters(17);%RCAN1-p production enhancement by F
    betaRK         = parameters(18);%RCAN1-p production enhancement by F
    Nc         = parameters(19);%Threshold activated NFAT causes inflammation

    %variable
    dydt = zeros(6,1);
 
    F       = y(1);
    NFs     = y(2);
    R       = y(3);
    Rp      = y(4);
    T       = y(5);
    FK506   = y(6);
 
   
    dydt(1) = a*F - deltaTF*T*F;
    dydt(2) = (kN*(1+betaNF*F)*(1+betaNR*Rp)*(1-NFs))/((1+yNK*FK506)*(1+yNR*R))- (deltaN*NFs);
    dydt(3) = kRN*NFs-SwitchStatus(1)*kR*(1+betaRF*F)*(1+betaRK*FK506)*R-deltaR*R;
    dydt(4) = SwitchStatus(1)*kR*R*(1+betaRF*F)*(1+betaRK*FK506) - deltaRP*Rp;
    dydt(5) = SwitchStatus(2)*kT*NFs - deltaT*T;
    dydt(6) = 0;%FK506 levels remains constant

    if(y(1)<1e-6)
        dydt(1)=0;
    end
    
    
end

function [value,isterminal,direction] = HybridSwitchingEvent(t,y,Rc,Nc,TimeAtThreshold,timedelay,SwitchStatus)
    
    ActivationTime = TimeAtThreshold + timedelay; 
    
    if SwitchStatus(1) == 0 && SwitchStatus(2) == 0 %R and N both OFF
        value = [y(3)- Rc; y(2)- Nc; t-ActivationTime;];
        isterminal = [1;1;1]; % 0 = don't stop
        direction = [+1;+1;0];  % 0 = both negative and positive direction
    elseif SwitchStatus(1) == 1 && SwitchStatus(2) == 0 %R on, N off
        value = [y(3) - Rc; y(2) - Nc; t-ActivationTime;];
        isterminal = [1;1;1]; % 0 = don't stop
        direction = [-1;+1;0];  % 0 = both negative and positive direction 
    elseif SwitchStatus(1) == 0 && SwitchStatus(2) == 1 %R off N on
        value = [y(3) - Rc; y(2) - Nc; t-ActivationTime;];
        isterminal = [1;1;1]; % 0 = don't stop
        direction = [+1;-1;0];  % 0 = both negative and positive direction 
    elseif SwitchStatus(1) == 1 && SwitchStatus(2) == 1 %R and N both ON
        value = [y(3) - Rc; y(2) - Nc; t-ActivationTime;];
        isterminal = [1;1;1]; % 0 = don't stop
        direction = [-1;-1;0];  % 0 = both negative and positive direction 
    end

     
end

function SwitchStatus = ChangeSwitchState(EventID_history,SwitchStatus)
        
        if EventID_history == 1
            disp('Event=1,error!');
        end
    % Reversible switch (R)
        if EventID_history == 3 %that means time delay has reached, change R switch
            if SwitchStatus(1) == 0 % if initially switched off then switch on
                SwitchStatus(1) = 1;
            elseif SwitchStatus(1) == 1 % if initially switched on then switch off
                SwitchStatus(1) = 0;
            end
        % Reversible switch (N)
        elseif EventID_history == 2
            if SwitchStatus(2) == 0 % switch on Gata3 switch
                SwitchStatus(2) = 1;
            elseif SwitchStatus(2) == 1;
                SwitchStatus(2) = 0;
            end
        end 
        
        
    end
