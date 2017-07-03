%the fungal macrophage model 
function dy = FMMOdesJY(t,y,parameters) %where do we get the y? from calc model data when we call this fn
%parameter is an array with the specific sequence 
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
    dy = zeros(6,1);
 
    F       = y(1);
    NFs     = y(2);
    R       = y(3);
    Rp      = y(4);
    T       = y(5);
    FK506   = y(6);
 
    %sw = (R^x)/(k^x + R^x); DECIDE IF THE SWITCH IS ON OR OFF
    
 
 
    dy(1) = a*F - deltaTF*T*F;
    dy(2) = (kN*(1+betaNF*F)*(1+betaNR*Rp)*(1-NFs))/((1+yNK*FK506)*(1+yNR*R))- (deltaN*NFs);
    if(y(3)>Rc) %switch on
        dy(3) = kRN*NFs-(kR*(1+betaRF*F)*(1+betaRK*FK506)+deltaR)*R;
        dy(4) = kR*R*(1+betaRF*F)*(1+betaRK*FK506) - deltaRP*Rp;
    else
        dy(3) = kRN*NFs - deltaR*R;
        dy(4) = -1*deltaRP*Rp;
    end
     
    if(y(2)>Nc)  %If N* is > threshold then T switch ON, then dT equation becomes;
        dy(5) = kT*NFs - deltaT*T;
    else    %if T switch OFF
        dy(5) = - deltaT*T;
    end
    dy(6) = 0;%FK506 levels remains constant
 
    if(y(1)<1e-6)
        dy(1)=0;
    end
    
    %disp('FMMODE-ing..');
 
end
