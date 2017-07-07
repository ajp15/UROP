%Plot results
function PrintResults(modelData_F,modelData_NFs,states,parameters,Switch_state,Switch_state_NFs,plotSwitch)
    
    constants = CalculateConstantsJY();
    
    if plotSwitch==1
        Switch_state(1,:)=30.*Switch_state(1,:);
        Switch_state(2,:)=1.*Switch_state(2,:);

        Switch_state_NFs(1,:)=90.*Switch_state_NFs(1,:);
        Switch_state_NFs(2,:)=(1/2).*Switch_state_NFs(2,:);
    end
    
    %Fungal Growth data*
    times = constants.data.Fgrowth.t;
    data = constants.data.Fgrowth.F_MOI3_withM;
    figure();
    plot(modelData_F(:,1),modelData_F(:,2), 'b', times, data, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('Fungus (F)','FontSize', constants.fs);
    title('Fungal Growth data','FontSize', constants.fs)

% TNFalpha data  
    t = constants.data.TNFalpha.t;
    TNFalpha = constants.data.TNFalpha.UT;
    
    figure();
    if plotSwitch==1
        plot(modelData_F(:,1),100.*Switch_state(2,:),'.y','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    end
    hold on;
    plot(modelData_F(:,1),modelData_F(:,6), 'b', t, TNFalpha, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
    hold on;
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('TNFalpha (T)','FontSize', constants.fs);
    title('TNFalpha (T) time course','FontSize', constants.fs)
 
    %NFAT: has to redefine IC
    t = constants.data.actNFAT.t;
    actNFAT = constants.data.actNFAT.UT;
    
    figure();
    if plotSwitch==1
        plot(modelData_NFs(:,1),Switch_state_NFs(2,:),'.y','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    end
    hold on;
    plot(modelData_NFs(:,1),modelData_NFs(:,3), 'b', t, actNFAT, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
    hold on;
    plot([0,6],[parameters(19),parameters(19)], ':k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2) 
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('active NFAT (N*)','FontSize', constants.fs);
    title('active NFAT (N*) time course','FontSize', constants.fs)

    figure,
    if plotSwitch==1
        plot(modelData_F(:,1),Switch_state(2,:),'.y','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    end
    hold on;
    plot(modelData_F(:,1),modelData_F(:,3), 'b', 'LineWidth',constants.lw)
    hold on;
    plot([0,6],[parameters(19),parameters(19)], ':k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2) 
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('active NFAT (N*)','FontSize', constants.fs);
    title('active NFAT (N*) time course MOI=3','FontSize', constants.fs)
   
    %RCAN1 data: plot model against experimental
    t = constants.data.RCAN1.t;
    RCAN1levels = constants.data.RCAN1.RCAN1levels;
    pRCAN1levels = constants.data.RCAN1.pRCAN1levels;
    totalRCAN1levels = RCAN1levels + pRCAN1levels;
    
    figure();
    if plotSwitch==1
        plot(modelData_F(:,1),Switch_state(1,:), '.y', 'MarkerSize', constants.ms,'LineWidth',constants.lw/10)
    end
    hold on;
    plot(modelData_F(:,1),modelData_F(:,4)+modelData_F(:,5), 'b','MarkerSize', constants.ms,'LineWidth',constants.lw)
    hold on
    plot(modelData_F(:,1),modelData_F(:,4), 'c',modelData_F(:,1),modelData_F(:,5), 'r' ,'MarkerSize', constants.ms,'LineWidth',constants.lw)
    hold on;
    plot([0,6],[parameters(14),parameters(14)], ':k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    hold on;
    plot(t, totalRCAN1levels, 'bo', t, RCAN1levels, 'co', t, pRCAN1levels, 'ro', 'MarkerSize', constants.ms/2,'LineWidth',constants.lw/2)
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('RCAN1levels','FontSize', constants.fs);
    title('RCAN1levels time course','FontSize', constants.fs)

end