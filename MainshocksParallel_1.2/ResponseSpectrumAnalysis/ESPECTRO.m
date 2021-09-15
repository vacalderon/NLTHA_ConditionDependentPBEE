clc; clear all; close all
EQ_List={'RSN95_MANAGUA_A-ESO090.AT2.g3','RSN210_LIVERMOR_A-A3E146.AT2.g3','RSN212_LIVERMOR_A-DVD156.AT2.g3','RSN214_LIVERMOR_A-KOD180.AT2.g3','RSN215_LIVERMOR_A-SRM070.AT2.g3','RSN230_MAMMOTH.I_I-CVK090.AT2.g3','RSN231_MAMMOTH.I_I-LUL000.AT2.g3','RSN232_MAMMOTH.I_I-MLS254.AT2.g3','RSN285_ITALY_A-BAG000.AT2.g3','RSN286_ITALY_A-BIS000.AT2.g3','RSN289_ITALY_A-CTR000.AT2.g3','RSN290_ITALY_A-MER000.AT2.g3','RSN291_ITALY_A-VLT000.AT2.g3','RSN292_ITALY_A-STU000.AT2.g3','RSN368_COALINGA.H_H-PVY045.AT2.g3','RSN543_CHALFANT.B_B-BEN270.AT2.g3','RSN544_CHALFANT.B_B-LAD180.AT2.g3','RSN545_CHALFANT.B_B-BPL070.AT2.g3','RSN546_CHALFANT.B_B-SHE009.AT2.g3','RSN547_CHALFANT.B_B-ZAK270.AT2.g3','RSN589_WHITTIER.A_A-ALH180.AT2.g3','RSN590_WHITTIER.A_A-ALT000.AT2.g3','RSN592_WHITTIER.A_A-CAM009.AT2.g3','RSN593_WHITTIER.A_A-ARL180.AT2.g3','RSN594_WHITTIER.A_A-NHO180.AT2.g3','RSN595_WHITTIER.A_A-JAB207.AT2.g3','RSN598_WHITTIER.A_A-TUJ262.AT2.g3','RSN599_WHITTIER.A_A-FLO020.AT2.g3','RSN602_WHITTIER.A_A-BUE250.AT2.g3','RSN608_WHITTIER.A_A-WAT180.AT2.g3','RSN611_WHITTIER.A_A-CAS000.AT2.g3','RSN613_WHITTIER.A_A-BAD000.AT2.g3','RSN614_WHITTIER.A_A-BIR090.AT2.g3','RSN615_WHITTIER.A_A-DWN180.AT2.g3','RSN616_WHITTIER.A_A-FAI095.AT2.g3','RSN620_WHITTIER.A_A-GLP177.AT2.g3','RSN622_WHITTIER.A_A-COM140.AT2.g3','RSN625_WHITTIER.A_A-ING000.AT2.g3','RSN626_WHITTIER.A_A-116270.AT2.g3','RSN627_WHITTIER.A_A-BLD000.AT2.g3','RSN629_WHITTIER.A_A-CCN000.AT2.g3','RSN632_WHITTIER.A_A-CYP053.AT2.g3','RSN633_WHITTIER.A_A-VER083.AT2.g3','RSN634_WHITTIER.A_A-FLE144.AT2.g3','RSN635_WHITTIER.A_A-PEL000.AT2.g3','RSN637_WHITTIER.A_A-FIG058.AT2.g3','RSN638_WHITTIER.A_A-WST000.AT2.g3','RSN639_WHITTIER.A_A-OBR270.AT2.g3','RSN640_WHITTIER.A_A-GR2090.AT2.g3','RSN641_WHITTIER.A_A-STN020.AT2.g3','RSN642_WHITTIER.A_A-W70000.AT2.g3','RSN645_WHITTIER.A_A-OR2010.AT2.g3','RSN649_WHITTIER.A_A-BRC000.AT2.g3','RSN650_WHITTIER.A_A-RIM015.AT2.g3','RSN652_WHITTIER.A_A-DEL000.AT2.g3','RSN661_WHITTIER.A_A-ANG000.AT2.g3','RSN663_WHITTIER.A_A-MTW000.AT2.g3','RSN665_WHITTIER.A_A-NWH180.AT2.g3','RSN667_WHITTIER.A_A-STC090.AT2.g3','RSN671_WHITTIER.A_A-PKC000.AT2.g3','RSN672_WHITTIER.A_A-KAG045.AT2.g3','RSN683_WHITTIER.A_A-OLD000.AT2.g3','RSN690_WHITTIER.A_A-GRN180.AT2.g3','RSN691_WHITTIER.A_A-SMA270.AT2.g3','RSN692_WHITTIER.A_A-EJS048.AT2.g3','RSN694_WHITTIER.A_A-CO2092.AT2.g3','RSN695_WHITTIER.A_A-RO3000.AT2.g3','RSN697_WHITTIER.A_A-GLE180.AT2.g3','RSN698_WHITTIER.A_A-SYL000.AT2.g3','RSN699_WHITTIER.A_A-SAY045.AT2.g3','RSN700_WHITTIER.A_A-TAR000.AT2.g3','RSN705_WHITTIER.A_A-SOR225.AT2.g3','RSN718_SUPER.A_A-IVW090.AT2.g3','RSN942_NORTHR_ALH090.AT2.g3','RSN948_NORTHR_CAM009.AT2.g3','RSN949_NORTHR_ARL090.AT2.g3','RSN950_NORTHR_NHO180.AT2.g3','RSN951_NORTHR_JAB220.AT2.g3','RSN954_NORTHR_TUJ262.AT2.g3','RSN956_NORTHR_BPK090.AT2.g3','RSN962_NORTHR_WAT180.AT2.g3','RSN963_NORTHR_ORR090.AT2.g3','RSN964_NORTHR_CAS000.AT2.g3','RSN966_NORTHR_BAD000.AT2.g3','RSN967_NORTHR_BIR090.AT2.g3','RSN968_NORTHR_DWN090.AT2.g3','RSN970_NORTHR_FAI095.AT2.g3','RSN974_NORTHR_GLP177.AT2.g3','RSN976_NORTHR_COM140.AT2.g3','RSN981_NORTHR_ING000.AT2.g3','RSN982_NORTHR_JEN022.AT2.g3','RSN983_NORTHR_JGB022.AT2.g3','RSN984_NORTHR_116090.AT2.g3','RSN985_NORTHR_BLD090.AT2.g3','RSN988_NORTHR_CCN090.AT2.g3','RSN990_NORTHR_LAC090.AT2.g3','RSN991_NORTHR_CYP053.AT2.g3','RSN992_NORTHR_VER090.AT2.g3','RSN993_NORTHR_FLE144.AT2.g3','RSN995_NORTHR_PEL090.AT2.g3','RSN997_NORTHR_FIG058.AT2.g3','RSN998_NORTHR_WST000.AT2.g3','RSN999_NORTHR_OBR090.AT2.g3','RSN1001_NORTHR_GR2090.AT2.g3','RSN1003_NORTHR_STN020.AT2.g3','RSN1004_NORTHR_SPV270.AT2.g3','RSN1006_NORTHR_UCL090.AT2.g3','RSN1017_NORTHR_BRC000.AT2.g3','RSN1018_NORTHR_RIM015.AT2.g3','RSN1024_NORTHR_DEL000.AT2.g3','RSN1039_NORTHR_MRP090.AT2.g3','RSN1041_NORTHR_MTW000.AT2.g3','RSN1044_NORTHR_NWH090.AT2.g3','RSN1048_NORTHR_STC090.AT2.g3','RSN1052_NORTHR_PKC090.AT2.g3','RSN1070_NORTHR_GRN180.AT2.g3','RSN1072_NORTHR_SMA090.AT2.g3','RSN1076_NORTHR_EJS030.AT2.g3','RSN1077_NORTHR_STM090.AT2.g3','RSN1082_NORTHR_RO3000.AT2.g3','RSN1083_NORTHR_GLE170.AT2.g3','RSN1086_NORTHR_SYL090.AT2.g3','RSN1087_NORTHR_TAR090.AT2.g3','RSN1094_NORTHR_SOR225.AT2.g3','RSN1186_CHICHI_CHY014-N.AT2.g3','RSN1187_CHICHI_CHY015-N.AT2.g3','RSN1193_CHICHI_CHY024-E.AT2.g3','RSN1194_CHICHI_CHY025-E.AT2.g3','RSN1195_CHICHI_CHY026-E.AT2.g3','RSN1197_CHICHI_CHY028-E.AT2.g3','RSN1198_CHICHI_CHY029-E.AT2.g3','RSN1201_CHICHI_CHY034-N.AT2.g3','RSN1202_CHICHI_CHY035-E.AT2.g3','RSN1203_CHICHI_CHY036-E.AT2.g3','RSN1204_CHICHI_CHY039-E.AT2.g3','RSN1205_CHICHI_CHY041-E.AT2.g3','RSN1206_CHICHI_CHY042-E.AT2.g3','RSN1208_CHICHI_CHY046-E.AT2.g3','RSN1209_CHICHI_CHY047-N.AT2.g3','RSN1211_CHICHI_CHY052-N.AT2.g3','RSN1227_CHICHI_CHY074-E.AT2.g3','RSN1231_CHICHI_CHY080-E.AT2.g3','RSN1234_CHICHI_CHY086-E.AT2.g3','RSN1235_CHICHI_CHY087-E.AT2.g3','RSN1236_CHICHI_CHY088-E.AT2.g3','RSN1238_CHICHI_CHY092-N.AT2.g3','RSN1240_CHICHI_CHY094-N.AT2.g3','RSN1244_CHICHI_CHY101-E.AT2.g3','RSN1245_CHICHI_CHY102-E.AT2.g3','RSN1246_CHICHI_CHY104-N.AT2.g3','RSN1258_CHICHI_HWA005-N.AT2.g3','RSN1261_CHICHI_HWA009-E.AT2.g3','RSN1262_CHICHI_HWA011-E.AT2.g3','RSN1264_CHICHI_HWA013-E.AT2.g3','RSN1265_CHICHI_HWA014-E.AT2.g3','RSN1268_CHICHI_HWA017-E.AT2.g3','RSN1269_CHICHI_HWA019-E.AT2.g3','RSN1270_CHICHI_HWA020-E.AT2.g3','RSN1277_CHICHI_HWA028-E.AT2.g3','RSN1278_CHICHI_HWA029-E.AT2.g3','RSN1279_CHICHI_HWA030-E.AT2.g3','RSN1280_CHICHI_HWA031-E.AT2.g3','RSN1281_CHICHI_HWA032-E.AT2.g3','RSN1282_CHICHI_HWA033-E.AT2.g3','RSN1283_CHICHI_HWA034-E.AT2.g3','RSN1284_CHICHI_HWA035-E.AT2.g3','RSN1285_CHICHI_HWA036-E.AT2.g3','RSN1286_CHICHI_HWA037-E.AT2.g3','RSN1287_CHICHI_HWA038-E.AT2.g3','RSN1288_CHICHI_HWA039-E.AT2.g3','RSN1289_CHICHI_HWA041-E.AT2.g3','RSN1290_CHICHI_HWA043-E.AT2.g3','RSN1294_CHICHI_HWA048-N.AT2.g3','RSN1295_CHICHI_HWA049-N.AT2.g3','RSN1296_CHICHI_HWA050-N.AT2.g3','RSN1297_CHICHI_HWA051-N.AT2.g3','RSN1300_CHICHI_HWA055-N.AT2.g3','RSN1301_CHICHI_HWA056-E.AT2.g3','RSN1302_CHICHI_HWA057-E.AT2.g3','RSN1303_CHICHI_HWA058-E.AT2.g3','RSN1304_CHICHI_HWA059-E.AT2.g3','RSN1350_CHICHI_ILA067-E.AT2.g3','RSN1380_CHICHI_KAU054-E.AT2.g3','RSN1482_CHICHI_TCU039-E.AT2.g3','RSN1483_CHICHI_TCU040-E.AT2.g3','RSN1488_CHICHI_TCU048-E.AT2.g3','RSN1489_CHICHI_TCU049-E.AT2.g3','RSN1490_CHICHI_TCU050-E.AT2.g3','RSN1491_CHICHI_TCU051-E.AT2.g3','RSN1492_CHICHI_TCU052-E.AT2.g3','RSN1493_CHICHI_TCU053-E.AT2.g3','RSN1494_CHICHI_TCU054-E.AT2.g3','RSN1495_CHICHI_TCU055-E.AT2.g3','RSN1496_CHICHI_TCU056-E.AT2.g3','RSN1497_CHICHI_TCU057-E.AT2.g3','RSN1498_CHICHI_TCU059-E.AT2.g3','RSN1499_CHICHI_TCU060-E.AT2.g3','RSN1500_CHICHI_TCU061-E.AT2.g3','RSN1501_CHICHI_TCU063-E.AT2.g3','RSN1503_CHICHI_TCU065-E.AT2.g3','RSN1504_CHICHI_TCU067-E.AT2.g3','RSN1505_CHICHI_TCU068-E.AT2.g3','RSN1506_CHICHI_TCU070-E.AT2.g3','RSN1507_CHICHI_TCU071-E.AT2.g3','RSN1508_CHICHI_TCU072-E.AT2.g3','RSN1509_CHICHI_TCU074-E.AT2.g3','RSN1510_CHICHI_TCU075-E.AT2.g3','RSN1511_CHICHI_TCU076-E.AT2.g3','RSN1512_CHICHI_TCU078-E.AT2.g3','RSN1513_CHICHI_TCU079-E.AT2.g3','RSN1515_CHICHI_TCU082-E.AT2.g3','RSN1517_CHICHI_TCU084-E.AT2.g3','RSN1520_CHICHI_TCU088-N.AT2.g3','RSN1521_CHICHI_TCU089-E.AT2.g3','RSN1528_CHICHI_TCU101-E.AT2.g3','RSN1529_CHICHI_TCU102-E.AT2.g3','RSN1530_CHICHI_TCU103-E.AT2.g3','RSN1531_CHICHI_TCU104-E.AT2.g3','RSN1532_CHICHI_TCU105-E.AT2.g3','RSN1533_CHICHI_TCU106-E.AT2.g3','RSN1534_CHICHI_TCU107-E.AT2.g3','RSN1535_CHICHI_TCU109-E.AT2.g3','RSN1536_CHICHI_TCU110-E.AT2.g3','RSN1538_CHICHI_TCU112-E.AT2.g3','RSN1539_CHICHI_TCU113-E.AT2.g3','RSN1540_CHICHI_TCU115-E.AT2.g3','RSN1541_CHICHI_TCU116-E.AT2.g3','RSN1542_CHICHI_TCU117-E.AT2.g3','RSN1543_CHICHI_TCU118-E.AT2.g3','RSN1545_CHICHI_TCU120-E.AT2.g3','RSN1546_CHICHI_TCU122-E.AT2.g3','RSN1547_CHICHI_TCU123-E.AT2.g3','RSN1549_CHICHI_TCU129-E.AT2.g3','RSN1550_CHICHI_TCU136-N.AT2.g3','RSN1551_CHICHI_TCU138-N.AT2.g3','RSN1552_CHICHI_TCU140-N.AT2.g3','RSN1553_CHICHI_TCU141-N.AT2.g3','RSN1554_CHICHI_TCU145-N.AT2.g3','RSN1581_CHICHI_TTN031-E.AT2.g3','RSN1586_CHICHI_TTN041-N.AT2.g3','RSN1605_DUZCE_DZC180.AT2.g3'};
dt_list=[0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.0029,0.0029,0.0024,0.0029,0.0029,0.0024,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.02,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.02,0.005,0.005,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.02,0.01,0.01,0.005,0.02,0.01,0.01,0.01,0.02,0.02,0.02,0.01,0.02,0.01,0.02,0.01,0.02,0.01,0.01,0.02,0.02,0.01,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.004,0.005,0.005,0.005,0.005,0.005,0.005,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.004,0.004,0.005,0.005,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.004,0.004,0.004,0.004,0.004,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.004,0.004,0.004,0.004,0.004,0.005,0.004,0.005];


% ======================== Program ESPECTRO.m =========================== %

%    Program for the analysis of strong motion data, it calculates:

%    - velocity and displacement time histories by trapezoidal integration

%    - Fourier spectrum via FFT

%    - Displacement, pseudo-velocity and pseudo-acceleration response
%      spectrum for any values of damping ratio. Equation of motion is
%      solved using the duhamel integral.

%    - Husid plot, significant duration and bracketed duration

%               LUIS A. MONTEJO (lumontv@yahoo.com.ar)

%         uptades available at www.geocities.com/lumontv/eng

%      DEPARTMENT OF CIVIL, CONSTRUCTION AND ENVIROMENTAL ENGINEERING

%                  NORTH CAROLINA STATE UNIVERSITY

%                       last updated: 02/20/2007

% ======================================================================= %

neq=length(EQ_List);

for i=1:neq
    % Input data:
    clearvars -except EQ_List neq dt_list i
    name =char(EQ_List(i)) ;                                % identifies actual work, the output file will be name.xls

    dt   = dt_list(i);                                  % time step of accelerogram [sec]
    zi   = [0.05];                    % vector with the damping ratios
    nom  = char(EQ_List(i));          % name of earthquake file
    go   = 1;                         % acceleration of gravity in the same untis that
                                      % your eq. record is, ig your record is in g's
                                      % then input 1

    Tmax = 10;			              % final period for spectra [sec]
    Ti   = 0.001;                     % initial period for spectra [sec]
    dT   = 0.1;                       % interval for natural periods


    % ==============================================================
    % ============================================ END OF INPUT DATA
































    % ========================================================================
     % ========================================================================



    nzi  = length(zi);                      % number of damping ratios
    g    = 981;                             % acceleration of gravity [cm/s^2]
    addpath('C:\ConditionDependentPBEE\GroundmotionSelection\Response Spectrum Analysis\Mainshocks_RS');            % directory with the accelerograms
    terr = load ([nom]);             % read earthquake data file load ([nom,'.g3']
    [nr,nc]  = size(terr); 	                % columns and rows of data file
    nt       = nr*nc;		     		    % original number of data points
    xg(1:nt) = (g/go)*terr';                % copy accelerogram in a vector
    np = length(xg);
    Tu   = (np-1) * dt;	    	        % final time of accelerogram
    t    = 0: dt: Tu;				        % vector with sampled times

    %==========================================================================
    %========== velocities and displacements via trapezoidal integration:

    vel   = dt*cumtrapz(xg);
    despl = dt*cumtrapz(vel); 

%     figure;subplot(3,1,1); plot( t,xg/g, 'b','LineWidth',1); 
%     grid on; axis tight; 
%     title(['acceleration, velocity and displacement time histories for ',nom],'FontSize',16);
%     ylabel('accel. [%g]','FontSize',16);
% 
%     subplot(3,1,2); plot( t,vel,'b', 'LineWidth',1); grid on;
%     axis tight; ylabel('vel. [cm/s]','FontSize',16);
% 
%     subplot(3,1,3); plot( t,despl,'b' ,'LineWidth',1); grid on;
%     axis tight; xlabel('time [s]','FontSize',16); ylabel('displ. [cm]','FontSize',16);

    %==========================================================================
    %========== response spectrum via Duhamel integral:

    T    = Ti: dT: Tmax;     	                  % vector with natural periods

    nper   = length(T);						          % number of natural periods
    SD     = zeros(nzi,nper);				          % rel. displac. spectrum
    PSV    = zeros(nzi,nper);				          % pseudo-vel. spectrum
    PSA    = zeros(nzi,nper);				          % pseudo-acc. spectrum

    for k = 1 : nzi
        for j = 1 : nper
            wn       = 2*pi/T(j);
            ub       = duhamel(wn,zi(k),1,dt,nt,0,0,-xg);
            SD(k,j)    = max( abs(ub) );
        end
        PSV(k,:) = (2*pi./T) .* SD(k,:);                   % pseudo-vel. spectrum
        PSA(k,:) = (2*pi./T).^2 .* SD(k,:);  	          % pseudo-accel. spectrum
    end

%     figure; plot(T,SD); grid on; axis tight;title(['Displacement Spectrum of ',nom],'FontSize',16);
%     ylabel('displacement [cm]','FontSize',16); xlabel('period [s]','FontSize',16);
% 
%     figure; plot(T,PSV); grid on; axis tight;title(['Pseudo-Velocity Spectrum of ',nom],'FontSize',16);
%     ylabel('PSV [cm/s]','FontSize',16); xlabel('period [s]','FontSize',16);
% 
%     figure; plot(T,PSA./g); grid on; axis tight;title(['Pseudo-Acceleration Spectrum of ',nom],'FontSize',16);
%     ylabel('PSA [g]','FontSize',16); xlabel('period [s]','FontSize',16);


    %==========================================================================
    %============= Husid plot, significant duration:

    Ia = zeros(1,nt);
    for n = 1 : nt
        Ia(n) = pi/(2*g)*dt*trapz( xg(1:n).^2 );
    end
    AI  = Ia(np);
    Ia  = Ia/Ia(np);
    x1  = 0.05*ones(1,np);
    x2  = 0.95*ones(1,np);
    [xmin,imin] = min( abs(Ia-0.05) );
    [xmax,imax] = min( abs(Ia-0.95) );
    dur = t(imax) - t(imin);

%     figure;  plot( t,Ia, t,x1, t,x2, t(imin),Ia(imin),'o', t(imax),Ia(imax),'o','LineWidth',2,...
%                     'MarkerEdgeColor','m',...
%                     'MarkerSize',8);
%     grid on; axis tight; title(['Husid plot of ',nom],'FontSize',16);
%     xlabel('Time [s]','FontSize',16); ylabel('Normalized intensity','FontSize',16);
%     text(2,0.07,'5%','FontSize',16); text(2,0.97,'95%','FontSize',16); 
%     text(Tu/2.7,0.55,['significant duration: ',num2str(dur),' s'],'FontSize',16)

%     pc=find(abs(xg/g)>0.05);
%     last=length(pc);
%     tc1=t(pc(1)); acc1 = xg(pc(1)); 
%     tc2=t(pc(last)); acc2 = xg(pc(last));
%     duration=tc2-tc1;
% 
%     cs=ones(1,np);
%     cs=0.05.*cs;
% 
%     figure;plot( t,abs(xg/g), 'b',t,cs,'r',tc1,abs(acc1)/g,'o',tc2,abs(acc2)/g,'o','LineWidth',1,...
%                     'MarkerEdgeColor','r',...
%                     'MarkerSize',8); 
%     grid on; axis tight; 
%     title(['bracketed duration of ',nom],'FontSize',16);
%     ylabel('abs(accel.) [%g]','FontSize',16);xlabel('Time [s]','FontSize',16);
%     text(Tu/2.5,0.95*max(abs(xg))/g,['bracketed duration: ',num2str(duration),' s'],'FontSize',16)

    %==========================================================================
    %============= Fourier spectrum via FFT:


    wny = pi/dt;					          % Nyquist frequency: rad/sec
    dw  = 2*pi / Tu;                          % frequency interval: rad/sec
    w   = 0.00001: dw: nt/2*dw;               % vector with frequencies in rad/sec
    f   = w/(2*pi);                           % vector with frequencies in cycles/sec                  
    Amp = dt * abs( fft(xg) );                % calculate the FT of the earthquake
    Amp = Amp(1:length(f));

%     figure; plot( f, Amp ); grid on; 
%     axis tight; title(['Fourier spectrum of ',nom],'FontSize',16)
%     xlabel('Frequency f [Hz]','FontSize',16); ylabel('Amplitude','FontSize',16);

    % =========================================================================
    % =========================================================================


    THS = [t' xg'./g vel' despl' Ia'];
    PGA = max(abs(xg./g));
    PGV = max(abs(vel));
    PGD = max(abs(despl));


    fid = fopen([name,'.csv'],'w');

    fprintf(fid, ' \n');
    fprintf(fid, 'PGA:,  %4.2f g\n',PGA);
    fprintf(fid, 'PGV:,  %4.2f cm/s\n',PGV);
    fprintf(fid, 'PGD:,  %4.2f cm\n',PGD);
    fprintf(fid, 'Arias Intensity:,  %5.2f cm/s\n',AI);
    fprintf(fid, 'Significant duration:,  %5.2f s\n',dur);
%     fprintf(fid, 'Bracketed  duration:,  %5.2f s\n',duration);
    fprintf(fid, ' \n');
    for k =1:nzi
        res = [T'  SD(k,:)'  PSV(k,:)'  PSA(k,:)'/g];
        fprintf(fid, ' \n');
        fprintf(fid, 'Damping ratio for spectra:  %3.2f\n\n',zi(k));
        fprintf(fid, ' \n');
        fprintf(fid, 'period [s]\tSD. [cm]\tPSV [cm/s]\tPSA [g]\n'); 
        fprintf(fid, '%3.4f\t, %4.3f\t, %4.3f\t, %4.3f\n',res');
        fprintf(fid, ' \n');

    end
    fprintf(fid, 'acceleration, velocity, displacement and normalized Arias intensity time histories:\n');
    fprintf(fid, ' \n');
    fprintf(fid, 'time [s]\taccel. [g]\tvel. [cm/s]\tdispl. [cm]\tAI\n'); 
    fprintf(fid, '%3.3f\t%1.7f\t%4.7f\t%4.7f\t%1.7f\n',THS');
    fprintf(fid, ' \n');
    fclose(fid);
end