% Ignorieren Sie die folgende Zeile. Überall wo ZuVeraendernderAusdruck im folgenden auftaucht 
% sind Änderungen durch Sie erforderlich.
ZuVeraendernderAusdruck=nan;

% Konstanten
G_K_max = ZuVeraendernderAusdruck;  % Maximale Kaliumleitfähigkeit(mS)
G_Na_max = ZuVeraendernderAusdruck; % Maximale Natriumleitfähigkeit (mS)
G_L = ZuVeraendernderAusdruck;      % Leckleitfähigkeit (mS)
C_m = ZuVeraendernderAusdruck;      % Membrankapazität (µF)

% Spannungen
E_K = ZuVeraendernderAusdruck;  % Nernstsspannung für Kalium (mV)
E_Na = ZuVeraendernderAusdruck; % Nernstsspannung für Natrium (mV)
E_L = ZuVeraendernderAusdruck;  % Nernstsspannung für Leckstrom (mV)

% initial values
n_0 = 0.32;
m_0 = ZuVeraendernderAusdruck;
h_0 = ZuVeraendernderAusdruck;
Vm_0 = ZuVeraendernderAusdruck;

% time parameters
deltat = 0.0005; %time step (ms)
tend = 10; %end of calculation (ms)
timesteps = 0:deltat:tend;

% Stimulusstromparameter
istim = ZuVeraendernderAusdruck; % Stimulusamplitude (µA)
tstim = ZuVeraendernderAusdruck; % Stimulusdauer (ms)

% Preallokation der Matrix für Raten, Gates, Ströme und Spannungen
nTimesteps = length(timesteps);
alpha_n = zeros(nTimesteps,1); 
beta_n = zeros(nTimesteps,1);
n = zeros(nTimesteps,1);
ndot = zeros(nTimesteps,1);
I_K = zeros(nTimesteps,1);
alpha_m = zeros(nTimesteps,1); 
beta_m = zeros(nTimesteps,1);
alpha_h = zeros(nTimesteps,1); 
beta_h = zeros(nTimesteps,1);
m = zeros(nTimesteps,1);
mdot = zeros(nTimesteps,1);
h = zeros(nTimesteps,1);
hdot = zeros(nTimesteps,1);
I_Na = zeros(nTimesteps,1);
I_L = zeros(nTimesteps,1);
Vm = zeros(nTimesteps,1);
Vmdot = zeros(nTimesteps,1);

%% Simulation
% Berechnung der Gatingvariablen, Ströme und Spannung

for 
    
    % Stimulusstrom
    
    % Kaliumstrom    
    
    % Natriumstrom    
    
    % Leckstrom    
    
    % Spannung Vm
    
            
end    
    


%% Plotten
subplot(3,1,1)
plot(timesteps,Vm)
xlabel('Time (ms)'), ylabel('Transmembrane voltage (mV)')

subplot(3,1,2)
plot(timesteps,n), hold on, plot(timesteps, m, 'g'), hold on, plot(timesteps,h,'r'),legend('n','m','h')
xlabel('Time (ms)')
ylabel('Gating variables')

subplot(3,1,3)
plot(timesteps,I_K), hold on, plot(timesteps, I_Na,'r'), hold on, plot(timesteps,I_L,'g'), legend('I_K','I_{Na}','I_L')
xlabel('Time (ms)')
ylabel('Current density (\muA/cm^2)') 
