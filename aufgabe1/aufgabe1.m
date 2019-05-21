% Ignorieren Sie die folgende Zeile. Überall wo ZuVeraendernderAusdruck im folgenden auftaucht 
% sind Änderungen durch Sie erforderlich.
ZuVeraendernderAusdruck=nan;

% Parameter
clampVoltages = linspace(-100,40,15); % Berechne die Gates für diese Spannungen(mV) 

% Preallokierung des Vektors für Ratenkonstanten
nVoltages = length(clampVoltages); % Anzahl an Clampspannungen
alpha_n = zeros(nVoltages, 1); 
beta_n = zeros(nVoltages, 1);
alpha_m = zeros(nVoltages, 1);
beta_m = zeros(nVoltages, 1);
alpha_h = zeros(nVoltages, 1);
beta_h = zeros(nVoltages, 1);


%% Simulation 

% Berechnung der Ratenkonstanten für jede Clampspannung

for iVoltage = 1:nVoltages
     alpha_n(iVoltage) = 0.01*(-(clampVoltages(iVoltage)+55)/(exp(-(clampVoltages(iVoltage)+55)/(10))-1));
     beta_n(iVoltage) = 0.125*exp(-(clampVoltages(iVoltage)+65)/(80));
     
     alpha_m(iVoltage) = 0.1*(-(clampVoltages(iVoltage)+40)/(exp(-(clampVoltages(iVoltage)+40)/(10))-1));        
     beta_m(iVoltage) = 4*exp(-(clampVoltages(iVoltage)+65)/(18));
     
     alpha_h(iVoltage) = 0.07*exp(-(clampVoltages(iVoltage)+65)/(20));
     beta_h(iVoltage) = 1/(exp(-(clampVoltages(iVoltage)+35)/(10))+1);
end
alpha_n(isnan(alpha_n))=1;
beta_n(isnan(beta_n))=1;
alpha_m(isnan(alpha_m))=1;
beta_m(isnan(beta_m))=1;
alpha_h(isnan(alpha_h))=1;
beta_h(isnan(beta_h))=1;
 
%% Plotting

figure(); 
subplot(3,1,1), plot(clampVoltages, alpha_n, 'r-x'); % plot data points (crosses) and interpolating solid line in red
hold on, plot (clampVoltages, beta_n, '-x');
legend('alpha_n', 'beta_n');
xlabel('Transmembrane voltage (mV)')
ylabel('Rate constant (1/ms)')

subplot(3,1,2), plot(clampVoltages, alpha_m, 'r-x');
hold on, plot (clampVoltages, beta_m, '-x');
legend('alpha_m', 'beta_m');
xlabel('Transmembrane voltage (mV)')
ylabel('Rate constant (1/ms)')

subplot(3,1,3),plot(clampVoltages, alpha_h, 'r-x');
hold on, plot (clampVoltages, beta_h, '-x');
legend('alpha_h', 'beta_h');
xlabel('Transmembrane voltage (mV)')
ylabel('Rate constant (1/ms)')
