%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clearvars; clc;close all;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload(fullfile('Nasa','NasaThermalDatabase.mat'),{'Diesel','O2','N2','CO2','H2O'});
%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 0;
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355;
CaIVC = -135;
CaEVO = 149;
CaEVC = -344;
CaSOI = -3.2;
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams
%% Fuel and data file
Fuel     = FuelMix();                                   % edit the table in Functions/FuelMix.m
FullName = fullfile('Data','ExampleDataSet.txt');       % measurement file for this fuel
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
%% Plotting 
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
pp = plot(Ca,p/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')
%% pV-diagram
V = CylinderVolume(Ca(:,iselect),Cyl);
f2 = figure(2);
set(f2,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice
subplot(2,1,1)
plot(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                      % Matter of taste
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram, cycle 10'})
subplot(2,1,2)
loglog(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0.5 50]);
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('pV-diagram (log-log), cycle 10')





%% ===================== Engine and fuel results =====================
% --- Inputs (change per experiment)
RPM        = 1500;      % engine speed [rpm]
mdot_fuel  = NaN;       % measured fuel mass flow [kg/s]; NaN = estimate from heat release
gammaC     = 1.35;      % constant ratio of specific heats (first estimate)

% --- Volume and mean cycle
Vall  = CylinderVolume(Ca, Cyl);             % volume for all cycles [m^3]
Vd    = pi/4*Cyl.Bore^2*Cyl.Stroke;          % displacement volume [m^3]
Cam   = Ca(:,1);                             % crank angle of one cycle [deg]
pm    = mean(p, 2);                          % mean pressure over all cycles [Pa]
Vm    = CylinderVolume(Cam, Cyl);            % volume of one cycle [m^3]

% --- Work, IMEP and power
Wcyc = zeros(1, Ncycles);
for i = 1:Ncycles
    Wcyc(i) = trapz(Vall(:,i), p(:,i));      % closed integral of p dV [J]
end
W     = mean(Wcyc);   dW    = std(Wcyc);     % work per cycle and its scatter [J]
IMEP  = W/Vd;         dIMEP = dW/Vd;         % [Pa]
P     = W*RPM/120;                           % indicated power, 4-stroke [W]

% --- Apparent heat release (constant gamma), filtered mean cycle
pf    = smoothdata(pm, 'sgolay', 21);        % Savitzky-Golay, 21 points = 4.2 deg
dCa   = Cam(2) - Cam(1);
dpdCa = gradient(pf, dCa);                   % [Pa/deg]
dVdCa = gradient(Vm, dCa);                   % [m^3/deg]
aROHR = gammaC/(gammaC-1)*pf.*dVdCa + 1/(gammaC-1)*Vm.*dpdCa;   % [J/deg]

iw    = Cam >= CaIVC & Cam <= CaEVO;         % closed part of the cycle (valves shut)
Caw   = Cam(iw);
aHR   = cumtrapz(Caw, aROHR(iw));            % cumulative heat release [J]
[~, imax]      = max(aHR);
[aHRmin, imin] = min(aHR(1:imax));
Qhr   = aHR(imax) - aHRmin;                  % heat released by combustion [J]
frac  = (aHR - aHRmin)/Qhr;                  % burned fraction [-]
CA10  = Caw(find(frac(imin:imax) >= 0.10, 1) + imin - 1);
CA50  = Caw(find(frac(imin:imax) >= 0.50, 1) + imin - 1);
CA90  = Caw(find(frac(imin:imax) >= 0.90, 1) + imin - 1);

% --- Fuel
if isnan(mdot_fuel)
    mf = Qhr/Fuel.LHV;                       % estimate: all released heat comes from the fuel
    mdot_fuel = mf*RPM/120;
    fprintf('No fuel flow given: fuel mass estimated from heat release (lower bound).\n');
else
    mf = mdot_fuel*120/RPM;                  % fuel mass per cycle [kg]
end
Qfuel     = mf*Fuel.LHV;                     % chemical energy of the fuel per cycle [J]
eta       = W/Qfuel;                         % indicated efficiency [-]
Qloss     = Qfuel - W;                       % energy not converted into work [J]
mair_st   = mf*Fuel.AFRst;                   % stoichiometric air per cycle [kg]
mCO2      = mf*Fuel.wC*44.01/12.011;         % CO2 from complete combustion per cycle [kg]
mdot_CO2  = mCO2*RPM/120;                    % [kg/s]
isfc      = (mdot_fuel/P)/(g/kWhr);          % indicated specific fuel consumption [g/kWh]
isCO2     = (mdot_CO2/P)/(g/kWhr);           % indicated specific CO2 [g/kWh]

% --- Print
fprintf('\n=========== Results: %s ===========\n', Fuel.Name);
fprintf('Work per cycle        W      = %8.1f +/- %.1f J\n', W, dW);
fprintf('IMEP                         = %8.2f +/- %.2f bar\n', IMEP/bara, dIMEP/bara);
fprintf('Indicated power       P      = %8.2f kW\n', P/1000);
fprintf('Heat released         Qhr    = %8.1f J\n', Qhr);
fprintf('Fuel energy           Qfuel  = %8.1f J\n', Qfuel);
fprintf('Energy not to work    Qloss  = %8.1f J\n', Qloss);
fprintf('Fuel mass per cycle   mf     = %8.2f mg\n', mf*1e6);
fprintf('Stoich. air per cycle        = %8.1f mg\n', mair_st*1e6);
fprintf('Indicated efficiency  eta    = %8.1f %%\n', eta*100);
fprintf('isfc                         = %8.0f g/kWh\n', isfc);
fprintf('isCO2                        = %8.0f g/kWh\n', isCO2);
fprintf('CA10 / CA50 / CA90           = %.1f / %.1f / %.1f deg\n', CA10, CA50, CA90);


degtime = 1/(6*RPM)*1000;                     % duration of 1 deg crank angle [ms]
fprintf('Fuel use                     = %8.2f kg/h = %.2f L/h\n', mdot_fuel*3600, mdot_fuel*3600/Fuel.rho*1000);
fprintf('CO2 emitted                  = %8.2f kg/h\n', mdot_CO2*3600);
fprintf('Fuel per cycle               = %8.1f mg\n', mf*1e6);
fprintf('Of 100 J of fuel energy      = %8.0f J becomes work, %.0f J is lost\n', 100*eta, 100*(1-eta));
fprintf('Ignition delay (SOI to CA10) = %8.1f deg = %.2f ms\n', CA10 - CaSOI, (CA10 - CaSOI)*degtime);
fprintf('Burn duration (CA10 to CA90) = %8.1f deg = %.2f ms\n', CA90 - CA10, (CA90 - CA10)*degtime);
fprintf('CA50 after TDC               = %8.1f deg (handbook: should be closer to TDC)\n', CA50);
