 

% Daten einlesen
load('RatenMischung.mat')

% Ignorieren Sie die folgende Zeile. Überall wo ZuVeraendernderAusdruck im folgenden auftaucht 
% sind Änderungen durch Sie erforderlich.
ZuVeraendernderAusdruck= nan;  

% Clampspannungen
clampVoltages = linspace(-100,40,15); % Loesung einfuegen .... Festlegen der Voltageclampspannungen in mV

% Ausgabe des Formates der eingelesenen Raten
groesse=sprintf('Groesse der eingelesenen Datei:\n Zeilen: %i Spalten: %i', size(Raten,1), size(Raten,2));
disp(groesse)


% Daten plotten
plot(clampVoltages, Raten)
% Legende hinzufügen
legend

% Klare Verläufe zuordnen
alpha_m=Raten(:,ZuVeraendernderAusdruck);
beta_m=Raten(:,ZuVeraendernderAusdruck);

% Position der Lücke suchen. Hierbei können die Funktionen isnan() und find() hilfreich sein.
loc_nans=ZuVeraendernderAusdruck;  %Index der Lücke im Wertebereich

% Mittelwert aus Werten vor und hinter Lücke berechnen
Mittel= ZuVeraendernderAusdruck;  %Lösung einfügen.... Mitte berechnen
% Lücke schließen
alpha_m(loc_nans)= Mittel;


% Restliche Verläufe zuordnen
alpha_h=Raten(:,ZuVeraendernderAusdruck);
beta_h=Raten(:,ZuVeraendernderAusdruck);
alpha_n=Raten(:,ZuVeraendernderAusdruck);
beta_n=Raten(:,ZuVeraendernderAusdruck);


% Verläufe getrennt plotten
figure;
subplot()
%Lösung einfügen

legend('alpha_n','beta_n')
xlabel('Spannung in mV')
ylabel('Ratenkonstante in 1/ms')
% Lösung einfügen

legend('alpha_m','beta_m')
xlabel('Spannung in mV')
ylabel('Ratenkonstante in 1/ms')
% Lösung einfügen

legend('alpha_h','beta_h')
xlabel('Spannung in mV')
ylabel('Ratenkonstante in 1/ms')

% Speichern der Daten im gewünschten Format
n.alpha=ZuVeraendernderAusdruck;
n.beta=ZuVeraendernderAusdruck;
% Lösung ergänzen