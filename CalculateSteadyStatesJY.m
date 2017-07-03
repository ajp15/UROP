%the fungal macrophage model, steady state when R=off, N=off
%will be used as initial condition
function states = CalculateSteadyStatesJY(parameters,FK506)
    
    if nargin < 2
        FK506=0;
    end
    
    alpha           = parameters(1);%fungal growth rate
    kN         = parameters(2);%NFAT activation rate
    yNK         = parameters(3);%FK506 inhibition rate (of NFAT activation)
    deltaN          = parameters(4);%NFAT deactivation rate
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
    betaNF        = parameters(16);%NFAT activation enhancement by F
    betaRF         = parameters(17);%RCAN1-p production enhancement by F
    betaRK         = parameters(18);%RCAN1-p production enhancement by F
    Nc         = parameters(19);%Threshold activated NFAT causes inflammation
 
    if FK506 == 0
        %cdt when R switch is OFF, N off and Fss=0 (fungal is cleared)
        a = deltaN*yNR*kRN/(kN*deltaR);
        b = 1+(deltaN/kN);
        c = -1;
        Foff    = 0;
        NFsoff  = (1/(2*a))*(-b+sqrt(b^2-4*a*c)); %value of N when R switch is off
        Roff    = kRN*NFsoff/deltaR; %using NFoff to calc Roff
        Rpoff   = 0;
        Toff    = 0; 
        FK506off= 0;

        states.switchOff = [Foff,NFsoff,Roff,Rpoff,Toff,FK506off];
        states.Roff_Non_Fclear = [Foff,NFsoff,Roff,Rpoff,kT*NFsoff/deltaT,FK506off];


        %When R On N On (Fungus cleared) 
        a=(1/(kR+deltaR))*(betaNR*kR*kRN/deltaRP+deltaN*yNR*kRN/kN);
        b=1+deltaN/kN-betaNR*kR*kRN/(deltaRP*(kR+deltaR));
        c=-1;
        Fon    = 0;
        NFson  = (1/(2*a))*(-b+sqrt(b^2-4*a*c)); %calc the value of N when R is ON
        Ron    = (kRN*NFson)/(kR+deltaR); %value of R when R switch ON
        Rpon   = kR*Ron/deltaRP; %what we get when R switch is ON
        Ton    = kT*NFson/deltaT; %value when N is ON
        FK506on= 0;

        states.switchOn = [Fon,NFson,Ron,Rpon,Ton,FK506on];
        %When R on N off F cleared
        states.switchRonNoff_Fclear=[Fon,NFson,Ron,Rpon,0,FK506on];

        %When R Off N On (Fungus not cleared)
        Rp_RoffNon=0;
        Tss_RoffNon=alpha/deltaTF;
        NFs_RoffNon=Tss_RoffNon*deltaT/kT;
        R_RoffNon=kRN*NFs_RoffNon/deltaR;
        F_RoffNon= (1/betaNF)*(deltaN*NFs_RoffNon*(1+yNR*R_RoffNon)/(kN*(1-NFs_RoffNon))-1);
        FK506_RoffNon=0;

        states.Roff_Non_Fnotclc=[F_RoffNon,NFs_RoffNon,R_RoffNon,Rp_RoffNon,Tss_RoffNon,FK506_RoffNon];

        %When R on N on F not cleared

        T_RonNonFnotclc=alpha/deltaTF;
        NFs_RonNonFnotclc=deltaT*T_RonNonFnotclc/kT;
        a=(1-NFs_RonNonFnotclc)*kN*betaNF*betaRF*(kR+betaNR*kR*kRN*NFs_RonNonFnotclc/deltaRP);
        b=(1-NFs_RonNonFnotclc)*kN*(betaRF*(kR+betaNR*kR*kRN*NFs_RonNonFnotclc/deltaRP)+betaNF*(kR+deltaR+betaNR*kR*kRN*NFs_RonNonFnotclc/deltaRP))-deltaN*NFs_RonNonFnotclc*kR*betaRF;
        c=(1-NFs_RonNonFnotclc)*kN*(kR+deltaR+betaNR*kR*kRN*NFs_RonNonFnotclc/deltaRP)-deltaN*NFs_RonNonFnotclc*(kR+deltaR+yNR*kRN*NFs_RonNonFnotclc);
        F_RonNonFnotclc=(1/(2*a))*(-b+sqrt(b^2-4*a*c));
        R_RonNonFnotclc=kRN*NFs_RonNonFnotclc/(kR*(1+betaRF*F_RonNonFnotclc)+deltaR);
        Rp_RonNonFnotclc=kR*(1+betaRF*F_RonNonFnotclc)*R_RonNonFnotclc/deltaRP;

        states.Ron_Non_Fnotclc=[F_RonNonFnotclc,NFs_RonNonFnotclc,R_RonNonFnotclc,Rp_RonNonFnotclc,T_RonNonFnotclc,0];

     elseif FK506 == 1
%             
        % F cleared
        kR2=kR*(1+betaRK*FK506);
        kN2=kN/(1+yNK*FK506);
        a=(1/(kR2+deltaR))*(betaNR*kR2*kRN/deltaRP+deltaN*yNR*kRN/kN2);
        b=1+deltaN/kN2-betaNR*kR2*kRN/(deltaRP*(kR2+deltaR));
        c=-1;
        NFs_FK506_PI=(1/(2*a))*(-b+sqrt(b^2-4*a*c));
        R_FK506_PI=kRN*NFs_FK506_PI/(kR2+deltaR);
        F_FK506_PI=0;
        Rp_FK506_PI=kR*R_FK506_PI/deltaRP;
        T_FK506_PI=kT*NFs_FK506_PI/deltaT;
        
        states.switchOn=[F_FK506_PI,NFs_FK506_PI,R_FK506_PI,Rp_FK506_PI,T_FK506_PI,FK506];

        
        a = deltaN*yNR*kRN/(kN2*deltaR);
        b = 1+(deltaN/kN2);
        c = -1;
        F_Roff_Noff_FK506    = 0;
        NFs_Roff_Noff_FK506  = (1/(2*a))*(-b+sqrt(b^2-4*a*c)); %value of N when R switch is off
        R_Roff_Noff_FK506   = kRN*NFs_Roff_Noff_FK506/deltaR; %using NFoff to calc Roff
        Rp_Roff_Noff_FK506    = 0;
        T_Roff_Noff_FK506     = 0; 
        FK506_Roff_Noff_FK506 = 1;

        states.Roff_Noff_FK506 = [F_Roff_Noff_FK506 ,NFs_Roff_Noff_FK506 ,R_Roff_Noff_FK506 ,Rp_Roff_Noff_FK506 ,T_Roff_Noff_FK506 ,FK506_Roff_Noff_FK506 ];
        states.switchOff = [F_Roff_Noff_FK506 NFs_Roff_Noff_FK506 R_Roff_Noff_FK506 0 0 1];
 
%FK506: kNF--> kNF/(1+yNK*FK506);
%FK506: kRP--> kRP*(1+cKR*FK506);
 
    end
end
