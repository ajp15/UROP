function [SwitchState] = GetSwitchState(parameters, FK506)
    
    states = CalculateSteadyStatesJY(parameters,FK506);
    
    if states.switchOff(3)>parameters(14)
            Rswitch=1;
    else
            Rswitch=0;
    end
        
    if states.switchOff(2)>parameters(19)
            Nswitch=1;
    else
            Nswitch=0;
    end
        
    SwitchState=[Rswitch,Nswitch];
 end