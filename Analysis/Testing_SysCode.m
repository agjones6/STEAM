%% Testing System Code for the BWRX300 Type Reactor 
%  STEAM: Scaled Time-dependent Esbwr Analysis Model
%  Author:  Andy Jones
%  Updated: February 19, 2019 

% ------------------------- IMPORTING DATA --------------------------------
%--> Declaring filenames

    % **** THESE NEED TO BE CHANGED IN ORDER TO USE ON LOCAL MACHINE ****
    % The file path leading to GitHub should be changed to match user 
    filepath = "/Users/AndyJones/Documents/SPRING 2019/GitHub/SysCoX300/";
    my_file1 = filepath + "velocityMat.csv";
    my_file2 = filepath + "uMat.csv";
    my_file3 = filepath + "densMat.csv";
    my_file4 = filepath + "voidMat.csv"; 
    % my_file5 = filepath + "convMat.csv"; 
    
  
%--> Creating Matrices 
    % Raw inputs 
    raw.vel    = csvread(my_file1); %,0,0,[0,0,25,13]);
    raw.u      = csvread(my_file2);
    raw.dens   = csvread(my_file3); 
    raw.void   = csvread(my_file4); 
    % raw.conv = csvread(my_file5);
    
    

    % Velocity matrix 
    t     = raw.vel(:,1); 
    Psys  = raw.vel(:,2);
    dt    = raw.vel(:,3);
    vel   = raw.vel(:,4:end);

    % Internal Energy matrix
    Tn     = raw.u(:,4:end);
    TFD    = raw.u(:,2);
    q0     = raw.u(:,3);
    ag11   = raw.u(:,1);
    relPwr = q0/900e6;

    % Density Matrix 
    rhon   = raw.dens(:,4:end);
    tcvPos = raw.dens(:,2);
    stDem  = raw.dens(:,3);
    t2     = raw.dens(:,1);
%     kloop = raw.dens(:,2); 

    % Void Fraction Matrix
    agn    = raw.void(:,4:end);
    msteam = raw.void(:,1);
    mfeed  = raw.void(:,2);
%     react  = raw.void(:,3);

    % Convergence Matrix (my_file5)
    % conv   = raw.conv(:,4:end);

    % Relative velocity Matrix 
    my_file6 = filepath + "relVMat.csv";
    raw.relV   = csvread(my_file6); 
    N_I    = raw.relV(:,1);
    CPR    = raw.relV(:,2);
    N_Xe   = raw.relV(:,3);
    Xj     = raw.relV(:,4:end);
    
    % All Reactivity Inputs 
    %rho_net,rho_ex,rho_void,rho_fuel,rho_CR,rho_Xe
%     my_file7 = filepath + "reactivity.csv";
%     raw.areact = csvread(my_file7);
%     r_net  = raw.areact(1:end,1);
%     r_ex   = raw.areact(1:end,2);
%     r_void = raw.areact(1:end,3);
%     r_fuel = raw.areact(1:end,4);
%     r_CR   = raw.areact(1:end,5);
%     r_Xe   = raw.areact(1:end,6);
%     t3     = raw.areact(1:end,7);

    % Random Variables
    my_file8 = filepath + "random.csv";
    raw.rand = csvread(my_file8);
    q_p_avg = raw.rand(:,1);
    test = raw.rand(:,2);
    
%--> Error Checks 
    % Making sure the arrays are of equal size
%     if size(Psys,1) ~= size(Tn,1)
%         t     = [t(:,1); t(end,1)];
%         Psys  = [Psys(:,1); Psys(end,1)];
%         dt    = [dt(:,1); dt(end,1)];
%         vel   = [vel(:,:); vel(end,:)];
%     end 

    % Preallocating Size
    [y,j] = size(vel); 
    [y2,n] = size(Tn);
    str.j = cell(j,1);
    str.n = cell(j,1);

%--> Ease of Plotting 
    % Making cell arrays for strings used in legends
    for i = 1:j
        str.j{i} = "j" + i;
        if i ~= j
            str.n{i} = "n" + i;
        else
            str.n{i} = "Sat";
        end 
    end

    % Creating a saturation Line for internal energies
%     satl.u = zeros(y2,1);
%     satl.rho = zeros(y2,1);
%     for i = 1:y2
%         satl.u(i,1)   = XSteamUS('uL_P',Psys(i)); 
%         satl.rho(i,1) = XSteamUS('rhoL_P',Psys(i)); 
%     end

    % Font Size 
    fz = 5.5;
    
    % Nominal Steam Flow 
    nomSteam = 3906386.77595359;
    
    % Nominal Poisons 
    nomXe = 78794112268.469009;
    nomI  = 54589134913.952164;
    
    str.names = ["inlet plenum","core 1","core 2","core 3","core 4","core 5",...
                 "core 6","chimney","dryers/seperators","downcomer","steam dome"];
    
%% ------------------------- PLOTTING DATA --------------------------------

% ANS 1                                 DOPPLER AND NEUTRON POWER
figure
yyaxis left
plot(t3,r_fuel.*1e5,'LineWidth',1.2)
ylabel('reactivity (pcm)')
yyaxis right
plot(t2,relPwr,'LineWidth',1.2)
ylabel('(relative to nominal)') 
xlabel('time (hr)')
legend('Doppler Feedback','Neutron Power') 
grid on

%% ANS 2                                CR REACTIVITY
figure 
plot(t3,r_CR.*1e5,'LineWidth',1.2,'Color',[0.8500, 0.3250, 0.0980])
ylabel('reactivity (pcm)')
xlabel('time (hr)')
grid on

%% ANS                                  STEAM DEMAND  
figure
plot(t2,stDem,'LineWidth',1.2)
ylabel('(relative to nominal)')
xlabel('time (hr)')
ylim([0.4,1.1])
grid on


%% ANS 3                                PRESSURE AND TCV
figure 
yyaxis left
plot(t,Psys,'LineWidth',1.2)
ylabel('(psia)')
% ylim([1049.9,1050.1])
yyaxis right
plot(t2,tcvPos.*100,'LineWidth',1.2)
ylabel('valve position (% full)')
xlabel('time (hr)')
legend('System Pressure','Turbine Control Valve','Location','SouthWest')
grid on

%% ANS 4                                VOID REACTIVITY AND M_STEAM
figure 
yyaxis left
plot(t3,r_void.*1e5,'LineWidth',1.2)
ylabel('reactivity (pcm)')
yyaxis right
plot(t2,msteam./nomSteam,'LineWidth',1.2)
ylabel('steam flow rate (relative to nominal)')
xlabel('time (hr)')
legend('Void Reactivity','Steam Flow Rate','Location','Northeast')
grid on

%% ANS 5                                STEAM DEMAND VS STEAM OUTPUT
figure 
plot(t2,stDem,'LineWidth',1.1)
% title('Steam Output and Steam Demand')
xlabel('time (hr)')
hold on 
plot(t2,msteam./nomSteam,'LineWidth',1.1)
ylabel('(relative to nominal)')
legend('Steam Demand','Thermal Power','FontSize',fz+2)
grid on 


%% ANS 6                                ALL REACTIVITIES 
figure
plot(t3, r_void*1e5,'LineWidth',1.2)
hold on
plot(t3, r_fuel*1e5,'LineWidth',1.2)
plot(t3, r_CR*1e5,'LineWidth',1.2)
plot(t3, r_Xe*1e5,'LineWidth',1.2)
plot(t3, r_ex.*1e5,'LineWidth',1.2)
plot(t3, r_net.*1e5,'LineWidth',1.2)

legend('void','fuel','control rods','xenon','excess','net','Location','West')
xlabel('time (hr)')
ylabel('reactivity (pcm)')
grid on

%% ANS 7 (backup)                       ALL VELOCITIES 
figure
plot(t,vel(:,1:end-1)./3600,'LineWidth',1.2)
legend(str.names,'FontSize',fz+4,'Location','best')
xlabel('time (hr)')
ylabel('velocity (ft/s)')
% title('Junction Velocity')
grid on

%% ANS 8 (backup)                       TEMEPERATURES
figure
plot(t2,Tn(:,1:end),'LineWidth',1.2)
legend(str.n{1:end-1},'FontSize',fz+1)
xlabel('time (hr)')
ylabel('temperature (\circF)')
% title('Junction Velocity')
grid on

%% ANS 9 (backup)                       CPR 
figure 
plot(t, CPR,'LineWidth', 1.2)
title('Critical Power Ratio')
ylabel('ratio')
xlabel('time (hr)')
grid on

%% ANS 10                               HEAT LAG 
figure
plot(t2, msteam./nomSteam,'LineWidth',1.3)
% title('Heat Lag in the System')
xlabel('time (hr)')
hold on 
plot(t2,relPwr,'LineWidth',1.3)
plot(t2,stDem,'--k','LineWidth',0.9)
ylabel('(fraction of nominal)')
legend('Thermal Power','Neutron Power','Steam Demand','FontSize',fz+2)
grid on 

%% STARTUP 1                            Fission Poisons
figure 
plot(t2,N_Xe./nomXe,'LineWidth',1.2)
hold on
plot(t2,N_I./nomI,'LineWidth',1.2)
legend('Xenon','Iodine','Location','best')
xlabel('time (hr)')
ylabel('fraction to nominal')
% title('Fission Product Poison')
grid on

%% STARTUP 2                            VOID FRACTION 
figure
plot(t2,agn(:,:),'Linewidth',1.4)
title('Node Void')
xlabel('Time (hr)')
ylabel('')
legend(str.names,'FontSize',fz+2,'Location','best')
grid on

%% Simpler 4 Plot of Values 
figure 
set(gcf, 'Position',  [20, 50, 1250, 700])

%--> Plot of Velocities 
subplot(2,2,1) 
plot(t,vel(:,1:end-1)./3600,'LineWidth',1.2)
legend(str.j{1:end-1},'FontSize',fz)
xlabel('Time (hr)')
ylabel('Velocity (ft/s)')
title('Junction Velocity')
grid on

%--> Control Rod Reactivity, Relative Power, and Xenon Concentration 
subplot(2,2,2)
title('Reactivity')
plot(t3,1+(r_CR./max(abs(r_CR))),'LineWidth',1.1) 
ylabel('Reactivity (pcm)')
hold on
plot(t2,relPwr,'LineWidth',1.1)
plot(t2,N_Xe./nomXe,'LineWidth',1.1)
maxCRrho = string(round(max(abs(r_CR)) * 1e5));
legend("CR Reactivity [" + maxCRrho + "pcm]",'Power','Xenon Conc','Location','Best')
title('Control Rod Reactivity and Power Level')
xlabel('Time (hr)')
ylabel('Relative Value')
grid on 
% plot(t2,agn(:,:),'LineWidth',1.2)
% xlabel('Time (hr)')
% ylabel('Void')
% title('Nodal Vapor Void Fraction')
% legend(str.n{1:end-1},'FontSize',fz)
% grid on

%--> Plot of Steam and Feed Flow
subplot(2,2,3)
plot(t2,mfeed./nomSteam,'LineWidth',1.1)
hold on
plot(t2,msteam./nomSteam,'LineWidth',1.1)
xlabel('Time (hr)')
ylabel('mass flow rate (rel)')
legend('Feed','Steam','Location','Best')
title('Feed and Steam Mass Flow Rates')
grid on 

%--> Steam Output vs Steam Demand 
subplot(2,2,4)
plot(t2,stDem,'LineWidth',1.1)
title('Steam Output and Steam Demand')
xlabel('Time (hr)')
hold on 
plot(t2,msteam./nomSteam,'LineWidth',1.1)
ylabel('Fraction of Full Power Value')
legend('Steam Demand','Rx Power','FontSize',fz+2)
grid on 

%% Plot the Fission Poisons  
figure 
plot(t2,N_Xe./nomXe,'LineWidth',1.2)
hold on
plot(t2,N_I./nomI,'LineWidth',1.2)
legend('Xenon','Iodine','Location','best')
xlabel('time (hr)')
ylabel('fraction to nominal')
% title('Fission Product Poison')
grid on

%% 
    figure 
    plot(t3,r_CR.*1e5,'LineWidth',1.2)
    hold on
    plot(t3,r_Xe.*1e5,'LineWidth',1.2)
    ylabel('reactivity (pcm)')
    xlabel('time (hr)')
    legend('control rods','xenon','Location','best')
    grid on
%% Plotting All of the Variables of Interest in 6 Plots 
figure 
set(gcf, 'Position',  [20, 50, 1250, 700]) 

%--> Velocities 
subplot(2,3,1)
plot(t,vel(:,1:11)./3600)
xlabel('Time (hr)')
ylabel('Velocity (ft/s)') 
title('Velocity')
legend(str.j{1:11},'FontSize',fz,'Position',[0.3243,0.681963470319635,0.0388,0.142313546423135]);

%--> Temperature 
subplot(2,3,2)
plot(t2,Tn)
title('Node Temperature')
xlabel('Time (hr)')
ylabel('Temperature (\circF)')
hold on 
legend(str.n{1:end-1},'FontSize',fz,'Position',[0.605,0.683713850837138,0.0404,0.142313546423135])
hold off

%--> Density 
% subplot(2,3,3)
% plot(t,rhon(:,:))
% title('Node Density')
% xlabel('Time (hr)')
% ylabel('Density (lbm/ft^3)')
% hold on
% legend(str.n{1:end-1},'FontSize',fz,'Position',[0.8858,0.683713706326368,0.0404,0.142313546423135])
% hold off

%---> Flow in Mlbm/hr
subplot(2,3,3)
plot(t2,mfeed./1e6)
hold on
plot(t2,msteam./1e6)
xlabel('Time (hr)')
ylabel('mass flow rate (Mlbm/hr)')
legend('Feed','Steam','Location','Best')
title('Feed and Steam Mass Flow Rates')

%--> Void Fraction
subplot(2,3,4)
plot(t2,agn(:,:))
title('Node Void')
xlabel('Time (hr)')
ylabel('')
legend(str.n{1:end-1},'FontSize',fz,'Position',[0.3243,0.230209992566635,0.0404,0.13089802130898])
% ylim([0,max(max(agn))+0.01])

%--> System Pressure
subplot(2,3,5)
yyaxis left
plot(t,Psys,'LineWidth',1.1)
hold on
% plot([t(1),t(end)],[1000,1000])
title('Pressure')
xlabel('Time (hr)')
ylabel('P (psia)')
% ylim([999.5,1000.5])
hold off
yyaxis right
plot(t2,tcvPos,'LineWidth',1.1)
ylabel('Turbine Control Valve Position')
legend("Actual","TCV Position",'FontSize',fz+5,'Position',[0.415,0.45250684931507,0.0532,0.039573820395738])

%--> Power and Steam Demand
subplot(2,3,6)
plot(t2,stDem,'LineWidth',1.1)
title('Rx Power and Steam Demand')
xlabel('Time (hr)')
hold on 
plot(t2,relPwr,'LineWidth',1.1)
ylabel('Fraction of Full Power Value')
legend('Steam Demand','Rx Power','FontSize',fz+2)

% ramp   = ((relPwr(1) - relPwr(3296))/(t2(1) - t2(3296)))/60


%% Plot of Demand, Steam Output, and Rx Power 
figure
plot(t2, msteam./nomSteam,'LineWidth',1.3)
% title('Heat Lag in the System')
xlabel('time (hr)')
hold on 
plot(t2,relPwr,'LineWidth',1.3)
plot(t2,stDem,'--k','LineWidth',0.9)
ylabel('Fraction of Full Power Value')
legend('Thermal Power','Neutron Power','Steam Demand','FontSize',fz+2)
grid on 


%% Plot of Feed Temperature over time 
figure
plot(t2, TFD)

%% Plot of CPR over time 
figure 
plot(t2, CPR,'LineWidth', 1.2)
title('Critical Power Ratio')

ylabel('CPR')
xlabel('Time (hr)')

minCPR = min(CPR);

%% Plot of the Control Rod Reactivity Over Time 
figure 
title('Reactivity')
yyaxis left
plot(t3,r_CR*1e5,'LineWidth',1.1) 
ylabel('Reactivity (pcm)')
yyaxis right
plot(t2,relPwr,'LineWidth',1.1)
% legend('CR Reactivity','Relative Power')
title('Control Rod Reactivity and Power Level')
ylabel('Relative Power')
xlabel('Time (hr)')
grid on 

%% Plotting Mass Flow Rate
%---> Core Constants ****NOT PULLED FROM ACTUAL CODE****
D_o      = 0.041083333333333;
s        = 0.053333333333333;  
D_core   = 10.7222;
n_frods  = 22080;
n_rods   = 24000;
Ax       = (pi/4)*D_core^2 - ((pi/4)*D_o^2)*n_frods;
G_core   = rhon(:,2:7) .* vel(:,2:7);
m_core   = G_core .* Ax;
m_corekg = m_core/2.2046226;

 
figure
%---> Plotting Mass flow in [Mlbm/hr]
plot(t,m_core(:,1)./10^6)
xlabel('Time (hr)')
ylabel('Mlbm/hr')

%---> Plotting Mass flow in [kg/s]
% plot(t,m_corekg(:,1)./(60^2))
% xlabel('Time (s)')
% ylabel('kg/s')

title('Core Mass Flow Rate')

%% Plot of Feed and Steam Flow 
mfeed_kg  =  mfeed./2.2046226;
msteam_kg = msteam./2.2046226;

figure
%---> Flow in [kg/s]
% plot(t,mfeed_kg./60^2)
% hold on
% plot(t,msteam_kg./60^2)
% xlabel('Time (s)')
% ylabel('mass flow rate (kg/s)')

%---> Flow in Mlbm/hr
plot(t2,mfeed./1e6)
hold on
plot(t2,msteam./1e6)
xlabel('Time (hr)')
ylabel('mass flow rate (Mlbm/hr)')

legend('Feed','Steam','Location','NorthWest')
title('Feed and Steam Mass Flow Rates')

%% Plot of Average Void in the Core
avg_void = mean(agn(:,2:7),2);

figure
plot(t2,avg_void)
title('Average Core Void (Nodes 2-7)')

%% plot of dt over time 
figure
title('Change in Time (dt)')
plot(t,dt)

%% Solo Plot of Junction Velocity 
figure
title('Junction Velocity')
plot(t,vel(:,1:11)./3600)

%% Plot of Junction Quality 
figure
title('Junction Quality')
plot(t,Xj)
% ylim([0,1])

%% Plotting Pressure and Turbine Control Valve Seperately 
figure
plot(t,Psys,'.')
hold on
plot([t(1),t(end)],[1000,1000])
legend("Actual","Reference",'FontSize',fz,'Location','NorthWest')
title('System Pressure')
xlabel('Time (hr)')
ylabel('P (psia)')
grid on
hold off

figure
title('TCV Position')
plot(t2,tcvPos,'.')
ylabel('Turbine Control Valve Position')

%% Plot of Convergence for Debbugging 
% THE FORTRAN CODE DOES NOT OUTPUT A CSV FOR CONVERGENCE DEBUGGING 
% If this is desired, uncomment the my_file5 code at the top of the script
% after adding the functionality in fortran 

% %--> Defining Variables to make viewing easier 
% % Normalizing the convergence matrix 
% convL  = raw.conv(:,8).*0.0001;
% % Setting 0's equal to NaN's so they dont show up on the plot 
% convL(convL==0) = NaN;
% % Variables used to define bounds that are desired for viewing  
% a = numel(conv(:,13)); % End Bound
% b = 600;               % Range from end bound (ie end - b)
% 
% %--> Plotting convergence 
% figure 
% plot(conv(a-b:end,13),'.')
% hold on
% plot(convL(a-b:end),'.')
% grid on


% Real feed flow controller
% level - level_ref (0,1) = E1
% (mFD - mST)/mST_nom     = E2
%
% delta_mST / m_nom

% Change in reactor power based on steam demand
% mST_demand is controlled
%
% delta_RxPwr_rel = ( mST_demand / mSTnom ) - ( mST/mSTnom ) 

% instantaneous CPR calculation 
% 
%% Calculating TCV Full Open Loss Coefficient 
gc       = 32.17405*3600^2; 
Pcond    = 1 * 144;
Psys_c   = 1000 * 144;
rhog     = XSteamUS('rhoV_p',Psys_c/144);
Vst      = max(vel(:,12)); % ft/hr
Ksys     = 40; 
TCV_full = 2*gc*((Psys_c) - Pcond)/(rhog*Vst^2) - Ksys; 

%% Calculating Rod Insertion Limits 
%--> Knowns
dz_dt   = 28 * 0.00328084 * 60^2; % [mm/s] --> [ft/hr] "Max rate of change of the control rod position (DCD)"
h_CR    = 12.5; % [ft] "Height of a control rods"
nCR     = 74; % number of control rods

%--> Estimates
rho_1CR  = 1000/nCR; % [pcm] "total reactivity from one control Rod"
mvm_rods = 1/4; % "fraction of rods availible for control" 

%--> Calculations 
rho_tot = rho_1CR * nCR; % [pcm]  "total reactivity from all control rods"
drho_dz = mvm_rods*rho_tot/h_CR;  % [pcm/ft] "reactivtiy per foot of control rod"
drho_dt = drho_dz * dz_dt ; % [pcm / hr] "Change in reactivity over time"
fprintf("The Max Insertion limit is: %.3f (pcm/hr)\n",drho_dt)

%% 
% Notes

% Fradial * Pavg = Phighest rod 
% [Phighest rod]/H = q0_hot_avg 
% q0_hot_avg -> q0_hot_hot 

% HEAT TRANSFER MODEL 
% Can do a steady state conduction model across the fuel rod to get a temp
    % distribution. Can get the heat flux at the rod surface. Need a UA
    % value. Assume Nucleate Boiling (ie Tom Correlation) you can get a
    % distribution 
% q_pp -> q_ppp -> Tf(r) 
% q_pp * n * pi * Do * H = UA_rx * (T_rx - T_sat)
% UA = Q_rx / (T_rx - T_sat)
% T_rx = (1/(pi*R^2)) * int(T(r)*2*pi*r, dr, 0, R)

% M is the mass of the fuel in the reactor 
% The heat coming off of the fuel is used in the CPR calc 
% 

% Xe worth 
% In the notes online, Chapter 15 of Brunswick FSAR there is a list of
%   accidents. 
% Trips 
    % High Power 
    % 
    
% Design Objectives 
    % CPR ? 
    % 

% 806004976.00839758        38.746833013442775











