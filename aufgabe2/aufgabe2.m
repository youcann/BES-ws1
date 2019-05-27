
% Ignorieren Sie die folgende Zeile. Überall wo ZuVeraendernderAusdruck im folgenden auftaucht 
% sind Änderungen durch Sie erforderlich.
ZuVeraendernderAusdruck=nan;

% Konstanten
G_K_max = 36; % Maximale Kalium-Leitfähigkeit (mS)
G_Na_max = 120; % Maximale Natrium-Leitfähigkeit (mS)

% Spannungen (Nernstspannungen können als konstant betrachtet werden)
E_K = -88;         % Nernstsspannung für Kalium (mV) 
E_Na = 50;        % Nernstsspannung für Natrium (mV)
VpreStep = -65;    % Spannung vor Sprung (mV)
VpostStep = -110;   % Spannung nach Sprung (mV)

clampVoltages = -100:10:40; % Protokoll für Voltageclamp (mV)

% Initialisiserung
n_0 = 0.32;
m_0 = 0.053;
h_0 = 0.6;

% Zeit
deltat = 0.0005; %Zeitschritt (ms)
tend = 10; %Ende der Berechnung (ms)
timesteps = 0:deltat:tend;

% Preallokation der Matrix für die Raten, Gates, Ströme, und Spannung
nVoltages = length(clampVoltages); % Anzahl an clamp voltages
nTimesteps = length(timesteps);
alpha_n = zeros(nVoltages, nTimesteps); 
beta_n = zeros(nVoltages, nTimesteps); 
n = zeros(nVoltages, nTimesteps); 
ndot = zeros(nVoltages, nTimesteps);
I_K = zeros(nVoltages, nTimesteps); 
alpha_m = zeros(nVoltages, nTimesteps); 
beta_m = zeros(nVoltages, nTimesteps); 
m = zeros(nVoltages, nTimesteps); 
mdot = zeros(nVoltages, nTimesteps); 
alpha_h = zeros(nVoltages, nTimesteps); 
beta_h = zeros(nVoltages, nTimesteps); 
h = zeros(nVoltages, nTimesteps); 
hdot = zeros(nVoltages, nTimesteps); 
I_Na = zeros(nVoltages, nTimesteps);
%Startwerte
n(:,1)=n_0*ones(1,nVoltages);
m(:,1)=m_0*ones(1,nVoltages);
h(:,1)=h_0*ones(1,nVoltages);

Vm=zeros(nVoltages,nTimesteps);

%% Simulation
% Berechne Raten, Gates (Euler-1-Schritt), I_K und I_Na
for iVoltage=1:nVoltages
    for iTimestep=1:nTimesteps
        %Vm
        if (timesteps(iTimestep)<=0.3)
            Vm_current=VpreStep;
        elseif (timesteps(iTimestep)>0.3 && timesteps(iTimestep)<=7.3)
            Vm_current=clampVoltages(iVoltage);    
        elseif (timesteps(iTimestep)>7.3 && timesteps(iTimestep)<=10)
             Vm_current=VpostStep;   
        else
             Vm_current=NaN;   
        end
        Vm(iVoltage,iTimestep)=Vm_current;
        
        % Kalium
        I_K(iVoltage,iTimestep) = G_K_max * n(iVoltage,iTimestep).^4 * (Vm_current-E_K);

        % Natrium
        I_Na(iVoltage,iTimestep) = G_Na_max * m(iVoltage,iTimestep).^3 * h(iVoltage,iTimestep) * (Vm_current-E_Na);

        %Uebergangsraten
        alpha_n(iVoltage,iTimestep) = 0.01*(-(Vm_current+55)/(exp(-(Vm_current+55)/(10))-1));
        if(isnan(alpha_n(iVoltage,iTimestep))) alpha_n(iVoltage,iTimestep) = 1; end
        beta_n(iVoltage,iTimestep) = 0.125*exp(-(Vm_current+65)/(80));
        if(isnan(beta_n(iVoltage,iTimestep))) beta_n(iVoltage,iTimestep) = 1; end
     
        alpha_m(iVoltage,iTimestep) = 0.1*(-(Vm_current+40)/(exp(-(Vm_current+40)/(10))-1));
        if(isnan(alpha_m(iVoltage,iTimestep))) alpha_m(iVoltage,iTimestep) = 1; end
        beta_m(iVoltage,iTimestep) = 4*exp(-(Vm_current+65)/(18));
        if(isnan(beta_m(iVoltage,iTimestep))) beta_m(iVoltage,iTimestep) = 1; end
     
        alpha_h(iVoltage,iTimestep) = 0.07*exp(-(Vm_current+65)/(20));
        if(isnan(alpha_h(iVoltage,iTimestep))) alpha_h(iVoltage,iTimestep) = 1; end
        beta_h(iVoltage,iTimestep) = 1/(exp(-(Vm_current+35)/(10))+1);
        if(isnan(beta_h(iVoltage,iTimestep))) beta_h(iVoltage,iTimestep) = 1; end

        %nmh dots
        ndot(iVoltage,iTimestep) = alpha_n(iVoltage,iTimestep) .* (1-n(iVoltage,iTimestep)) - beta_n(iVoltage,iTimestep) .* n(iVoltage,iTimestep);
        mdot(iVoltage,iTimestep) = alpha_m(iVoltage,iTimestep) .* (1-m(iVoltage,iTimestep)) - beta_m(iVoltage,iTimestep) .* m(iVoltage,iTimestep);
        hdot(iVoltage,iTimestep) = alpha_h(iVoltage,iTimestep) .* (1-h(iVoltage,iTimestep)) - beta_h(iVoltage,iTimestep) .* h(iVoltage,iTimestep);
        
        if(iTimestep < nTimesteps)
            n(iVoltage,iTimestep+1) = deltat * ndot(iVoltage,iTimestep) + n(iVoltage,iTimestep);
            m(iVoltage,iTimestep+1) = deltat * mdot(iVoltage,iTimestep) + m(iVoltage,iTimestep);
            h(iVoltage,iTimestep+1) = deltat * hdot(iVoltage,iTimestep) + h(iVoltage,iTimestep);
        end    
        
    end
end


%% Plotten
% Plotten der Raten und Gates für eine Spannung (Index 14) und Plotten der Ströme für alle Spannungen

% Kalium
subplot(3,1,1)
title('14. Voltagestep')
plot(timesteps,alpha_n(14,:))
hold on
plot(timesteps,beta_n(14,:),'r')
legend('alpha_n','beta_n')
xlabel('Time (ms)')
ylabel('Rate constant (1/ms)')

subplot(3,1,2)
title('14. Voltagestep')
plot(timesteps, n(14,:))
xlabel('Time (ms)')
ylabel('Gating variable n')

subplot(3,1,3)
plot(timesteps,I_K)
xlabel('Time (ms)')
ylabel('I_{K} current density (\muA)')

% Natrium
figure()
subplot(3,1,1)
plot(timesteps,alpha_m(14,:))
hold on, plot(timesteps,beta_m(14,:),'r')
plot(timesteps,alpha_h(14,:))
plot(timesteps,beta_h(14,:))
legend('alpha_m','beta_m', 'alpha_h', 'beta_h')
xlabel('Time (ms)'),
ylabel('Rate constant (1/ms)')

subplot(3,1,2)
plot(timesteps, m(14,:))
hold on
plot(timesteps, h(14,:))
legend('m','h')
xlabel('Time (ms)')
ylabel('Gating variables m, h')

subplot(3,1,3)
plot(timesteps,I_Na)
xlabel('Time (ms)')
ylabel('I_{Na} current density (\muA)')
    
    
