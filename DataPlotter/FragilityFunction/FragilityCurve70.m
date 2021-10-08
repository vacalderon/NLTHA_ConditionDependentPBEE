%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORTH CAROLINA STATE UNIVERSITY%
%Designed by: Eng. Diego Martinez Pineda (ingdrmp@gmail.com)
%Checked by: Dr. Mervyn J. Kowalsky
%DETERMINATION OF FRAGILITY FUNCTIONS. Baker, J. W. (2015). “Efficient analytical fragility function fitting 
% using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.
clc, clear all, close all
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%I. INPUT DATA%%%
%1. Names Files
%NameOut='Fragility Bridges Kong';                      %Name Output File
InputFile='DataFragilityCurves';                        %Name of Excel File with Data for Fragility Curves
NumFile=70;                                             %Number associated to the input file
Bridges=[1;2;3;4;5;6;7;8;9;10];                         %Identification of Bridges to be analyzed
NEq=100;                                                %Number of Earthquakes
NRot=2;                                                 %Number of Rotations Earthquake

%2. Damping Models
IDR=3;                                                  %Type of Damping Model [o, 2, 6...]
NDM=4;                                                  %Number of Rayleigh Damping Models
DMFC=[1;2;4;6;7;9;12];                                  %Damping Model for Fragility Curves <IDR*NDM 
NameDM={'DM1. M+K0-1&2-5';'DM2. MODAL';'DM3. M+KS-1&2-5';'DM3. M+K0-1&10-5%';'DM5';'DM4. M+KS-1&10-5%'; 'DM5. K0-1&2-5%';'DM8';...
        'DM6. KS-1&2-5%'; 'DM10'; 'DM11.'; 'DM7. ZERO VISCOUS'};
%2= Modal Damping, 1=Initial Rayleigh 1&2, 4=Initial Rayleigh 1&10, 6=Secant Rayleigh 1&10,
%7= Initial Stiffness 1&2, 9=Secant Stiffness 1&2, 12=Secant Stiffness1&10, 15=Zero Viscous Damping
DMFIG=[2;6;9;12];                                       %Damping Models to see variation against DM1, in total 4

%3. Analysis Options
Ana=2;                                                  %1=Fragility Function for Yield Limit State
                                                        %2=Seismic Demand Model
Est=1;                                                  %1=Positive and Negative Displacements, 2=Abs Max Displacements                                              
ValProb=0.5;                                            %Value of Probability to obtain SdT1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Data From Excel file with information of Sd and Displacements for the bridges
namexcel=([InputFile,num2str(NumFile), '.xlsx']);       %Name excel file with input for Ruaumoko3D
[Sd1, TEXT1, MC1] = xlsread(namexcel, 'Sd-All');        %Cell array with numbers and strings Spectral Displacements (Sd)
[Input, TEXT, MC] = xlsread(namexcel, 'Inputs');        %Cell array with numbers and strings Demands in Z direction (Dz)
nr1=length(DMFC);                                       %Number of Damping Models for Fragility Analysis
nr2=length(DMFIG);                                      %Number of Damping Models for Variability Analysis
NSt=length(Bridges);                                    %Number of Structures to analyze
[dx,dy]=size(Input);

%Data For Fragility Functions
for i=1: dx
    hcol(i,1:5)=Input(i,2:6);                                   %Information of Height of Columns per Bridge
    Lspan(i,1:6)=Input(i,7:12);                                 %Information of Length of Spans per Bridge
    Ncol(i,1)=Input(i,13);                                      %Information of Column to Conduct Fragility Analysis
    Cap(i,1:3)=Input(i,14:16);                                  %Information of Capacities of column to Yield, Serviceability and Damage Control
    SdFig(i,1)=Input(i,17);                                     %Information of Capacities of column to Yield, Serviceability and Damage Control
    Ncol(i,:)=length(hcol(i,:)) - sum(isnan(hcol(i,:)));        %Number of Columns per Bridge
    Nspan(i,:) = length(Lspan(i,:)) - sum(isnan(Lspan(i,:)));   %Number of Spans per Bridge
end

%1. Figure 1: Displacement Ratios vs Fundamental Period
Ncols=IDR*NDM*NEq*NRot; 
figure
for i=1:nr2
    xdm=DMFIG(i);
    pdm1=(xdm-1)*NRot*NEq+1;                                    %Initial Position of DMi
    pdm2=(xdm)*NRot*NEq;                                        %Final Position of DMi
         
    for n=1:NSt
    xpb=Bridges(n);
    T(1,1:NRot*NEq)=Sd1(xpb+2,2);                               %Period Bridge
    NodeFC=Input(xpb,13);
    [Dz1, TEXT, MC] = xlsread(namexcel, ['Dem',num2str(xpb)]);  %Cell array with numbers and strings Demands in Z direction (Dz)
        
        for j=1:1*Ncols
        Dz(1,j)=max(Dz1(NodeFC+6,j+1), abs(Dz1(NodeFC+6,IDR*NDM*NEq*NRot+j+1)));    
        end
        DIS(1,1:NRot*NEq)=Dz(1,1:NRot*NEq);                     %Displacements for DM1
        Disp(1,1:NRot*NEq)=Dz(1,pdm1:pdm2);
    
        for k=1:NRot*NEq
        Dratio(1,k)=Disp(1,k)/DIS(1,k);
        end
        
        %3. Figure
        hold on
        subplot(2,2,i)
        scatter(T,Dratio,'filled')
        %title('FF-SERVICEABILITY','fontweight','bold','FontSize',22); 
        grid on
        NameBr(n,1)={['Bridge ',num2str(xpb)]};
        nam=(['Disp. Ratio',NameDM(xdm),'/DM1']);
        hx = xlabel('Fundamental Period (s)', 'Fontsize', 14);
        hy = ylabel(nam, 'Fontsize', 14);
    end   
    legend(NameBr,'Location','best','fontsize', 12)
    axis([2 4 0 5])
end
hold off
toc

%2. Figure 2: Displacement Ratios vs Ductility
figure
for i=1:nr2
    xdm=DMFIG(i);
    pdm1=(xdm-1)*NRot*NEq+1;                                    %Initial Position of DMi
    pdm2=(xdm)*NRot*NEq;                                        %Final Position of DMi
         
    for n=1:NSt
    xpb=Bridges(n);                                             %identification Bridge
    NodeFC=Input(xpb,13);                                       %Node of Critical Column
    Dyield=Input(xpb,14);                                       %Yield Displacement
    [Dz1, TEXT, MC] = xlsread(namexcel, ['Dem',num2str(xpb)]);  %Cell array with numbers and strings Demands in Z direction (Dz)
        
        for j=1:1*Ncols
        Dz(1,j)=max(Dz1(NodeFC+6,j+1), abs(Dz1(NodeFC+6,IDR*NDM*NEq*NRot+j+1)));    
        end
        DIS(1,1:NRot*NEq)=Dz(1,1:NRot*NEq);                     %Displacements for DM1
        Disp(1,1:NRot*NEq)=Dz(1,pdm1:pdm2);                     %Displacements for rest of DM
    
        for k=1:NRot*NEq
        Dratio(1,k)=Disp(1,k)/DIS(1,k);
        Du(1,k)=Disp(1,k)/Dyield;
        end
        
        %3. Figure
        hold on
        subplot(2,2,i)
        scatter(Du,Dratio,'filled')
        %title('FF-SERVICEABILITY','fontweight','bold','FontSize',22); 
        grid on
        NameBr(n,1)={['Bridge ',num2str(xpb)]};
        nam=(['Disp. Ratio',NameDM(xdm),'/DM1']);
        hx = xlabel('Displacement Ductility', 'Fontsize', 14);
        hy = ylabel(nam, 'Fontsize', 14);
    end   
    legend(NameBr,'Location','best','fontsize', 12)
    %axis([2 4 0 5])
end
hold off
toc

%3. Fragility Functions
for n=1:NSt
    xpb=Bridges(n);
%2. Spectral Displacement for Structure   
for j=1:NEq*NRot
    Sd(1,j)=Sd1(xpb+2,j+2);
    if Est==1
    Sd(1,NEq*NRot+j)=Sd1(xpb+2,j+2);
    end
end

%3. Determination of Collapse of Structure due to Earthquake
if Est==1
Ncols=2*IDR*NDM*NEq*NRot;                                       %Number of Columns of Demands per Structure
Num_gms=ones(1,2*NEq*NRot);                                     %Number of Ground Motions, they are 2, because it is included positive and negative values
else
Ncols=IDR*NDM*NEq*NRot;   
Num_gms=ones(1,NEq*NRot);                                       %Number of Ground Motions, they are 2, because it is included positive and negative values
end

[Dz1, TEXT, MC] = xlsread(namexcel, ['Dem',num2str(xpb)]);      %Cell array with numbers and strings Demands in Z direction (Dz)
CapYield=Input(xpb,14);
CapSer=Input(xpb,15);
CapDC=Input(xpb,16);
NodeFC=Input(xpb,13);  
Sdmax=Input(xpb,17);

    for j=1:1*Ncols
    if Est==1    
    Dz(1,j)=Dz1(NodeFC+6,j+1);
    else
    Dz(1,j)=max(Dz1(NodeFC+6,j+1), abs(Dz1(NodeFC+6,IDR*NDM*NEq*NRot+j+1)));    
    end
    %1. Yield
        if abs(Dz(1,j))>CapYield
        CollapseYield(1,j)=1;
        else
        CollapseYield(1,j)=0;
        end      
    
    %2. Serviceability
        if abs(Dz(1,j))>CapSer
        CollapseServ(1,j)=1;
        else
        CollapseServ(1,j)=0;
        end
    
    %3. Damage Control
        if abs(Dz(1,j))>CapDC
        CollapseDC(1,j)=1;
        else
        CollapseDC(1,j)=0;
        end  
    end

%3.1 Collapse for Damping Model
for i=1:nr1
    if Est==1
    %Positive Displacements
    z1=DMFC(i);
    y1=(z1-1)*NEq*NRot+1;
    y2=z1*NEq*NRot;
    CollapseYield1(1,1:NEq*NRot)=CollapseYield(1,y1:y2);
    CollapseServ1(1,1:NEq*NRot)=CollapseServ(1,y1:y2);
    CollapseDC1(1,1:NEq*NRot)=CollapseDC(1,y1:y2);
    
    %Negative Displacements
    x1=NEq*NRot+1;
    x2=y1+0.5*Ncols;
    x3=x2+NEq*NRot-1;
    CollapseYield1(1,x1:2*NEq*NRot)=CollapseYield(1,x2:x3);
    CollapseServ1(1,x1:2*NEq*NRot)=CollapseServ(1,x2:x3);
    CollapseDC1(1,x1:2*NEq*NRot)=CollapseDC(1,x2:x3);
    else
    z1=DMFC(i);
    y1=(z1-1)*NEq*NRot+1;
    y2=z1*NEq*NRot;    
    CollapseYield1(1,1:NEq*NRot)=CollapseYield(1,y1:y2);
    CollapseServ1(1,1:NEq*NRot)=CollapseServ(1,y1:y2);
    CollapseDC1(1,1:NEq*NRot)=CollapseDC(1,y1:y2); 
    end
    
    %Fragility Function
    %Theta: Median of fragility function, Beta: Lognormal standard deviation of fragility function
    [theta_hat_yield, beta_hat_yield] = fn_mle_pc(Sd(1,:), Num_gms, CollapseYield1(1,:));
    [theta_hat_serv, beta_hat_serv] = fn_mle_pc(Sd(1,:), Num_gms, CollapseServ1(1,:));
    [theta_hat_DC, beta_hat_DC] = fn_mle_pc(Sd(1,:), Num_gms, CollapseDC1(1,:));
    ThetaYield(1,i)=theta_hat_yield;    BetaYield(1,i)=beta_hat_yield;
    ThetaServ(1,i)=theta_hat_serv;  BetaServ(1,i)=beta_hat_serv;
    ThetaDC(1,i)=theta_hat_DC;    BetaDC(1,i)=beta_hat_DC;
end

%FIGURES
%Compute fragility functions using estimated parameters
x_IM = 0.000001:0.000001:Sdmax;                                                               %IM levels to plot fragility function at
nIM=length(x_IM);
%NameLeg(1,1)={'Observed fractions of collapse'};

for i=1:nr1
    z1=(i-1)*nIM+1;
    z2=i*nIM;
    x=DMFC(i);
    p_collapse_Yield(1,z1:z2) = normcdf((log(x_IM/ThetaYield(1,i)))/BetaYield(1,i));    %Compute fragility function using equation 1-Yield
    p_collapse_Serv(1,z1:z2) = normcdf((log(x_IM/ThetaServ(1,i)))/BetaServ(1,i));       %Compute fragility function using equation 1-Serviceability
    p_collapse_DC(1,z1:z2) = normcdf((log(x_IM/ThetaDC(1,i)))/BetaDC(1,i));             %Compute fragility function using equation 1-Damage Control
    NameLeg(i,1)=NameDM(x);                                                             %Legends
end

%Geometry Bridge
Ncolt=3*Ncol(xpb,1)+2;
Ybridge=zeros(1,Ncolt);

    for i=1:Ncol(xpb,1)+2
        if i==1; Coor(xpb,i)=0; else; Coor(xpb,i)=Lspan(xpb,i-1)+Coor(xpb,i-1);
        end
    end

    for i=1:Ncolt
        x=ceil((i-1)/3);
        if i==1; Xbridge(1,i)=0; elseif i==Ncolt; Xbridge(1,i)=max(Coor(xpb,:)); 
        else; Xbridge(1,i)=Coor(xpb,x+1);
        end
    end
    
    for i=1:Ncol(xpb,1)
        x=3*i;
        Ybridge(1,x)=-hcol(xpb,i);
    end
    Ymin=min(Ybridge)-2; Ymax=max(Ybridge)+2; tm=max(Xbridge);

%Figures
    figure
    %1. Geometry Bridge   
    subplot (2,2,1);  
    plot (Xbridge,Ybridge,'Linewidth',6,'color',[0 0 0]); title('LONGITUDINAL VIEW OF BRIDGE','fontweight','bold','FontSize',20); 
    xlabel('Length (m)','fontweight','bold','FontSize',18);
    ylabel('Height (m)','fontweight','bold','FontSize',18); grid on;
    set(gca, 'XLim', [0 tm], 'YLim', [Ymin Ymax]); 
    
    if Ana==1
    %2.1. Plot Resulting Fragility Functions-Yield    
    for i=1:nr1
        z1=(i-1)*nIM+1;
        z2=i*nIM;    
        hold on
        subplot(2,2,2)
        plot(x_IM,p_collapse_Yield(1,z1:z2), 'linewidth', 1.5)
        title('FF-YIELD','fontweight','bold','FontSize',24); 
        legend(NameLeg,'Location','best','fontsize', 10)
        grid on
        hx = xlabel('Intensity Measure: Sd(m)', 'Fontsize', 18);
        hy = ylabel('Prob. [P>LS-Yield]', 'Fontsize', 18);
        axis([0 Sdmax 0 1])
    end
    hold off

    else
    %2.2. Plot Relationship Enginering Demand Parameter vs Intensity Measure
    %x_text=1; 
    %y_text=1;
    for i=1:nr1
    %Positive Displacements
    z1=DMFC(i);
    y1=(1-1)*Ncols+(z1-1)*NEq*NRot+1;
    y2=(1-1)*Ncols+z1*NEq*NRot;
    Dzz1(1,1:NEq*NRot)=Dz(1,y1:y2);
    if Est==1
    %Negative Displacements
        x1=NEq*NRot+1;
        x2=y1+0.5*Ncols;
        x3=x2+NEq*NRot-1;
        Dzz1(1,x1:2*NEq*NRot)=abs(Dz(1,x2:x3));
        Sdf=Sd;
        Dzzf=Dzz1;
    elseif Est==2     
        Sdf=Sd;
        Dzzf=Dzz1;
    end
        hold on
        subplot(2,2,2)
        scatter(Sdf,Dzzf,'filled')
        r= corrcoef(Sdf, Dzzf);
        rfig(n,i)=r(1,2)^2;
        name1=(['DM',num2str(DMFC(i)),' r^2=',num2str(round(rfig(n,i),2))]);
        NameLeg2(i,1)={name1};  
        title('SEISMIC DEMAND MODELS','fontweight','bold','FontSize',22); 
        grid on
        hx = xlabel('Spectral Displacement, SdT1(m)', 'Fontsize', 16);
        hy = ylabel('Trans. Displ. (m)', 'Fontsize', 16);
    end
        hold off    
        legend(NameLeg2,'Location','best','fontsize', 10)
    end
        
    %3. Plot Resulting Fragility Functions-Serviceability
    for i=1:nr1
        z1=(i-1)*nIM+1;
        z2=i*nIM;    
        hold on
        subplot(2,2,3)
        prob1=round(p_collapse_Serv(1,z1:z2),2);
        B=(prob1==ValProb);
        xm=x_IM(B); 
        ProbServ(n,i)=xm(1,1);
        plot(x_IM,p_collapse_Serv(1,z1:z2), 'linewidth', 1.5)
        title('FF-SERVICEABILITY','fontweight','bold','FontSize',22); 
        legend(NameLeg,'Location','best','fontsize', 10)
        grid on
        hx = xlabel('Intensity Measure: Sd(m)', 'Fontsize', 16);
        hy = ylabel('Prob. [P>LS-Service]', 'Fontsize', 16);
        axis([0 Sdmax 0 1])
    end
    hold off

    %4. Plot Resulting Fragility Functions-Damage Control
    for i=1:nr1
        z1=(i-1)*nIM+1;
        z2=i*nIM;    
        hold on
        subplot(2,2,4)
        prob2=round(p_collapse_DC(1,z1:z2),2);
        C=(prob2==ValProb);
        xm2=x_IM(C); 
        ProbDC(n,i)=xm2(1,1);
        plot(x_IM,p_collapse_DC(1,z1:z2), 'linewidth', 1.5)
        title('FF-DAMAGE CONTROL','fontweight','bold','FontSize',22); 
        legend(NameLeg,'Location','best','fontsize', 10)
        grid on
        hx = xlabel('Intensity Measure: Sd(m)', 'Fontsize', 16);
        hy = ylabel('Prob. [P>LS-D.Control]', 'Fontsize', 16);
        axis([0 Sdmax 0 1])
    end
    hold off
    clearvars B C Dz1 Dz Xbridge Ybridge x_IM p_collapse_Yield p_collapse_Serv p_collapse_Serv 
end

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTIONS%
%Objective function to be optimized
function [theta, beta] = fn_mle_pc(Sd, num_gms, CollapseServ1)
% Initial guess for the fragility function parameters theta and beta. 
% These initial choices should not need revision in most cases, but they could be altered if needed.
x0 = [0.8 0.4];
%Run optimization
options = optimset('MaxFunEvals',1000, 'GradObj', 'off'); %Maximum 1000 iterations, gradient of the function not provided
x = fminsearch(@mlefit, x0, options, num_gms, CollapseServ1, Sd) ;
theta = x(1);
beta = x(2);
end

%Function to calculate likelihood
function [loglik] = mlefit(theta, num_gms, CollapseServ1, Sd)
if theta(1)<0 % don't let median of fragility function go below zero
    theta(1)=0;
end
% estimated probabilities of collapse, given the current fragility functionparameter estimates
p = normcdf(log(Sd), log(theta(1)), theta(2)); 

% likelihood of observing num_collapse(i) collapses, given num_gms
% observations, using the current fragility function parameter estimates
likelihood = binopdf(CollapseServ1', num_gms', p'); % 

% sum negative log likelihood (we take the nevative value because we want
% the maximum log likelihood, and the function is searching for a minimum)
loglik = -sum(log(likelihood));
end



