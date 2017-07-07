function jac=StateJac(t,y,parameters,IC_Switches)

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

    F       = y(1);
    NFs     = y(2);
    R       = y(3);
    Rp      = y(4);
    T       = y(5);
    FK506   = y(6); 
    
    jac=zeros(6,6);
    Rsw=IC_Switches(1);
    Nsw=IC_Switches(2);
    
jac=[a-T*deltaTF,0,0,0,-F*deltaTF,0;
     (betaNF*kN*(Rp*betaNR + 1)*(1-NFs))/((FK506*yNK + 1)*(R*yNR + 1)), -deltaN-(kN*(F*betaNF+1)*(Rp*betaNR+1))/((FK506*yNK+1)*(R*yNR+1)), (kN*yNR*(F*betaNF+1)*(Rp*betaNR+1)*(NFs-1))/((FK506*yNK+1)*(R*yNR+1)^2), (betaNR*kN*(F*betaNF+1)*(1-NFs))/((FK506*yNK + 1)*(R*yNR + 1)),0, (kN*yNK*(F*betaNF + 1)*(Rp*betaNR + 1)*(NFs - 1))/((FK506*yNK + 1)^2*(R*yNR + 1));
     -R*Rsw*betaRF*kR*(FK506*betaRK + 1), kRN, - deltaR - Rsw*kR*(F*betaRF + 1)*(FK506*betaRK + 1),0,0,-R*Rsw*betaRK*kR*(F*betaRF + 1);
     R*Rsw*betaRF*kR*(FK506*betaRK + 1), 0, Rsw*kR*(F*betaRF + 1)*(FK506*betaRK + 1), -deltaRP,0,R*Rsw*betaRK*kR*(F*betaRF + 1);
     0,Nsw*kT,0,0,-deltaT,0;
     0,0,0,0,0,0];
 
end
