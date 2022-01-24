clear all; close all; clc;

%------------------------------Cindy Iskandar------------------------------
%----------------------------------161862----------------------------------
%--------------------------------TC MCC 2019-------------------------------
%---------------------------Responsable : R.Ghosn--------------------------

p=tf('p');

%% Intro 

disp('Ce programme sert a determiner les caracteristiques de la MCC et de la commander grace a deux boucles (vitesse & courant) en cascades.');
disp('        ');
disp('La machine etudiee est la MCC 11-1.');
disp('_______________________________________');

%% Donnees constructeur

Pn=336e3;                           %en W 
Un=600;                             %en V
In=560;                             %en A
Ufn=360;                            %en A
r_n=0.92;                           %rendement nominal
Nn=840;                             %en tr/min
Omega_n=(2*pi/60)*Nn;               %en rad/s
Cn=3513;                            %en N.m
R_115=0.081;                        %en Ohms
L=0.62e-3;                          %en H
J=6.9;                              %en Kg.m2
fn=50;                              %en Hz
wn=2*pi*fn;                         %en rad/s

R=R_115;
%% Inducteur

%D'apres le catalogue de Leroysomer on trouve Pfn=5.5kW

Pfn=5.5e3;                          %en W

%On a que Pfn=Ufn^2/Rf

Rf=(Ufn^2)/Pfn;                     %en Ohms
Ifn=sqrt(Pfn/Rf);                   %en A
Lf =170;                            %en H (estimation)
phi_0=(Un-R_115*In)/Omega_n;        %en Wb
Mfd=phi_0/Ifn;                      %en H

%% Bilan de puissance 

%Calcul des pertes Joules nominales:
%On a Pjn=R(115)*In^2 on prend R(115) car Worst Case.

Pjn=R_115*(In)^2;                   %en W

%Calcul des pertes Fer nominales et des pertes Mecaniques nominales
%D'apres le bilan de puissances on a : Un*In = PJn+Pfen+Pmn+Cn*Omega_n
%D'ou :

Pmn=(3/5)*(Un*In-r_n*Un*In-Pjn);    %en W
Pfen=(2/3)*Pmn;                     %en W


%% Commande de la MCC

ksi_n=5;
to_elec=L/R_115;
kconv=Un/ksi_n;


%% Coef de frottements visqueux 

f=Pmn/(Omega_n)^2;

%% Elements de la machine 

disp('         ');
disp('Les parametres de la machine sont les suivants:  ');
disp('________________________________________________________');
disp('         ');
disp('Parametres de l''inducteur : ');
disp('         ');
fprintf('Pfn = %0.5f W\n', Pfn);
fprintf('Ifn = %0.5f A\n', Ifn);
fprintf('Ufn = %0.5f V\n', Ufn);
fprintf('phi_0 = %0.5f Wb\n', phi_0);
fprintf('Rf = %0.5f Ohms\n', Rf);
fprintf('Lf = %0.5f H\n', Lf);
fprintf('Mfd = %0.5f H\n', Mfd);
disp('         ');
disp('         ');
disp('Parametres de l''induit : ');
disp('         ');
fprintf('R = %0.5f Ohms\n', R_115);
fprintf('L = %0.5f H\n', L);
disp('         ');
disp('         ');
disp('Bilan de puissance : ');
disp('         ');
fprintf('Puissance electrique nominale : Pn = %0.5f W\n', Pn);
fprintf('Pertes joules nominales : Pjn = %0.5f W\n', Pjn);
fprintf('Puissance electromagnetique nominale : Pemn = %0.5f V\n', (Un*In-Pjn));
fprintf('Pertes fer nominales : Pfen = %0.5f W\n', Pfen);
fprintf('Pertes mecanique nominales : Pmn = %0.5f W\n', Pmn);
fprintf('Puissance utile (receuillie sur l''arbre) : Pu = %0.5f W\n', Cn*Omega_n);
disp('         ');
disp('___________________________________________________________');
disp('         ');


%% Couple resistant

a1=phi_0*In-f*Omega_n;
a2=(phi_0*In-f*Omega_n)/Omega_n;
a3=(phi_0*In-f*Omega_n)/Omega_n^2;

run=1;

while(run==1)
    while(isempty(run) || isstring(run) || (run~=1 && run~=2))
        run=input('Reponse : ');
        disp('              ');
    end
    disp('_______________________________________________________________');
    Tsim=20;
    disp('          ');
    disp('Voulez vous effectuer la simulation avec : ');
    disp('          ');
    disp('1. Le modele de la MCC ? ');
    disp('          ');
    disp('2. Les fonctions de transfert ? ');
    disp('          ');
    reponse=input('Reponse : ');
    while (isempty(reponse) || isstring(reponse) || (reponse~=1 && reponse~=2))
        disp('          ');
        reponse=input('Reponse : ');
    end 
    disp('          ');
    disp('_______________________________________________________________');
    if (reponse==1)
        disp('          ');
        disp('Choisissez la forme du couple resistant : ');
        disp('          ');
        disp('1. Couple resistant constant : Cr = Kc1 = cte ');
        disp('          ');
        disp('2. Couple resistant proportionnel a la vitesse : Cr = Kc2 * Omega ');
        disp('          ');
        disp('3. Couple resistant proportionnel au carre de la vitesse : Cr = Kc3 * Omega^2 ');
        disp('          ');
        x=input('Reponse : ');
        disp('          ');

        while(isempty(x) || isstring(x) || (x~=1 && x~=2 && x~=3))
            x=input('Reponse : ');
            disp('          ');
        end

        if (x==1)
            fprintf('1. Sachant que Kc1 < %0.5f N.m :', a1);
            kc1=input('Kc1 = ');
            kc2=0; kc3=0;
            while(isempty(kc1) || isstring(kc1) || kc1<0 || kc1 > a1)
                fprintf('1. Sachant que Kc1 < %0.5f N.m :', a1);
                kc1=input('Kc1 = ');
            end
        end
     
        if (x==2)
            disp('          ');
            fprintf('1. Sachant que Kc2 < %0.5f N.m.s/rad :', a2);
            kc2=input('Kc2 = ');
            kc1=0; kc3=0;
            while(isempty(kc2) || isstring(kc2) || kc2<0 || kc2 > a2)
                fprintf('1. Sachant que Kc2 < %0.5f N.m.s/rad :', a2);
                kc2=input('Kc2 = ');
            end
        end 

        if (x==3)
            disp('          ');
            fprintf('1. Sachant que Kc3 < %0.5f N.m :', a3);
            kc3=input('Kc3 = ');
            kc1=0; kc2=0;
            while (isempty(kc3) || isstring(kc3) || kc3<0 || kc3 > a3)
                fprintf('1. Sachant que Kc3 < %0.5f N.m :', a3);
                kc3=input('Kc3 = ');
            end
       end 

        disp ('         ');
        disp('_______________________________________________________________');
        
        k=kc2+f;

        %% Resistance de l'induit variable

        to=1;
        alfa=3.98e-3;
        delta_T=75;
        R0=R_115/(1+alfa*delta_T);
        disp('              ');
        disp('Voulez vous : ');
        disp('              ');
        disp('1. Une resistance de l''induit constante ?');
        disp('              ');
        disp('2. Une resistance de l''induit variable en fonction de la temperature ?');
        disp('              ');
        y=input('Reponse : ');
        disp('              ');

        while(isempty(y) || isstring(y) || (y~=1 && y~=2))
            y=input('Reponse : ');
            disp('              ');
        end 

        if (y==1)
            variable = 1;
        end 

        if (y==2)
            variable = 0;
        end 
    end
    
%     if (reponse==2)
%         disp('          ');
%         disp('Sachant que le couple resistant s''ecrit de la forme : Cr = kc * Omega ');
%         disp('          ');
%         fprintf('Choisir kc sachant que kc < %0.5f N.m.s/rad ',a2);
%         kc2=input('Kc = ');
%         while(isempty(kc2) || isstring(kc2) || kc2<0 || kc2 > a2)
%                 fprintf('Choisir kc sachant que Kc < %0.5f N.m.s/rad :', a2);
%                 kc2=input('Kc = ');
%         end
%         k=f+kc2; 
%     end
    %% Boucle de courant

    %Hi=(J*p+k)/(L*J*p^2+(R*J+k*L)*p+phi_0^2+k*R);
    Hi=(J*p+f)/(L*J*p^2+(R*J+f*L)*p+phi_0^2+f*R);

    %Parametres du filtre de Butterworth

    Rfi=150e3;                  %en Ohms
    Cfi=22e-9;                  %en F
    Hfi=1/(1+2*Rfi*Cfi*p+2*(Rfi*Cfi*p)^2);

    ks=1/240;                   %gain du shunt
    ka=ksi_n/(ks*In);           %mise a echelle

    %Calculons la fct de transfert du courant en BO:

    Hboi=ks*ka*kconv*Hi;
    Tibfc=to_elec/2.8;
    %Tibfc=to_elec/4;
    Kci=1/Tibfc;

    %% Boucle de Vitesse

    %to_mec=J/k;
    to_mec=J/f;
    
    %puisque la puissance de la machine est de 336 kW on veut que la
    %constante de temps mecanique soit de 3 s donc:
    to_mec_mod=to_mec/20;
    
    %Hvi=phi_0/(J*p+k);
    Hvi=phi_0/(J*p+f);
    ktachy=90/Omega_n;

    %ke*ktachy*wn=ksi_n
    ke=ksi_n/(ktachy*Omega_n);
    kwref=ksi_n/Nn;

    kcourant=In/ksi_n;
    Hbov=kcourant*Hvi*ktachy*ke;

    Tvbfc=to_mec_mod/4;
    Kcv=1/Tvbfc;
    
%% Simulation
    if (reponse==1)
        if (x==1)
            %Kci=Kci*0.25;
            Kci=Kci*0.15/3.75;
            %Kcv=Kcv*6;
            Kcv=Kcv*2.5;
        else
            Kci=Kci/3.5;
            Kcv=Kcv*7;
        end 
        sim('trial28732792392.slx');
        figure()
        plot(t,R_variable,'linewidth',2);
        hold on
        plot(t,R_ref,'linewidth',2);
        legend('Resistance variable','Resistance a 115 degrees');
        xlabel('Temps (s)');
        ylabel('Resistance (Ohms)'),grid;
        title('Variation de la resistance de l''induit en fonction du temps et de la Temperature');
    end
    
    if (reponse==2)
       sim('Boucle_Fct_Transfert.slx');
       figure()
       plot(temps,current,'linewidth',2);
       hold on
       plot(temps,current_ref,'linewidth',2);
       axis([0.9 1.5 0 800]);
       legend('Courant','Courant ref');
       xlabel('Temps (s)');
       ylabel('Courant de l''induit (A)');
       title('Variation du Courant de l''induit en fonction du temps');
       figure()
       plot(temps,speed,'linewidth',2);
       hold on
       plot(temps,speed_ref,'linewidth',2);
       legend('Omega','Omega ref');
       xlabel('Temps (s)');
       ylabel('Omega (rad/s)');
       title('Variation de la Vitesse en fonction du temps');
       figure()
       subplot(2,1,1)
       plot(temps,speed,'linewidth',2),grid;
       hold on
       plot(temps,speed_ref,'linewidth',2);
       legend('Omega','Omega ref');
       xlabel('Temps (s)');
       ylabel('Omega (rad/s)');
       title('Variation de la Vitesse en fonction du temps');
       subplot(2,1,2)
       plot(temps,current,'linewidth',2),grid;
       hold on
       plot(temps,current_ref,'linewidth',2);
       legend('Courant','Courant ref');
       xlabel('Temps (s)');
       ylabel('Courant de l''induit (A)');
       title('Variation du Courant de l''induit en fonction du temps');
    end 
    
    if (reponse==1)
        if (x==1)
            figure()
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            figure()
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
            figure()
            subplot(2,1,1)
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            axis([0.9 5 0 100]);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            subplot(2,1,2)
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
            figure()
            subplot(2,1,1)
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            subplot(2,1,2)
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
        end 
        if (x==2)
            if (kc2==0)
               figure()
                subplot(2,1,1)
                plot(t,vitesse,'linewidth',2),grid;
                hold on
                plot(t,vitesse_ref,'linewidth',2);
                axis([0.9 5 0 100]);
                legend('Omega','Omega ref');
                xlabel('Temps (s)');
                ylabel('Omega (rad/s)');
                title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
                subplot(2,1,2)
                plot(t,courant,'linewidth',2),grid;
                hold on
                plot(t,courant_ref,'linewidth',2);
                axis([0 5 0 800]);
                legend('Courant','Courant ref');
                xlabel('Temps (s)');
                ylabel('Courant (A)');
                title('Allure du courant en fct du temps avec le modele reel de la MCC '); 
            end
            figure()
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            figure()
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            axis([0 6.5 0 800]);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
            
            figure()
            subplot(2,1,1)
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            axis([0.9 5 0 100]);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            subplot(2,1,2)
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            axis([0.9 1.5 0 800]);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
            
            figure()
            subplot(2,1,1)
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            subplot(2,1,2)
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
        end
        if (x==3)
            figure()
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            figure()
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
            
            figure()
            subplot(2,1,1)
            plot(t,vitesse,'linewidth',2),grid;
            hold on
            plot(t,vitesse_ref,'linewidth',2);
            axis([0.9 5 0 100]);
            legend('Omega','Omega ref');
            xlabel('Temps (s)');
            ylabel('Omega (rad/s)');
            title('Allure de la vitesse en fct du temps avec le modele reel de la MCC ');
            subplot(2,1,2)
            plot(t,courant,'linewidth',2),grid;
            hold on
            plot(t,courant_ref,'linewidth',2);
            axis([0.9 1.5 0 800]);
            legend('Courant','Courant ref');
            xlabel('Temps (s)');
            ylabel('Courant (A)');
            title('Allure du courant en fct du temps avec le modele reel de la MCC ');
        end
    end
    
    disp('_______________________________________________________________');
    disp('               ');
    disp('Voulez vous : ');
    disp('               ');
    disp('1. Lancer le programme ?');
    disp('               ');
    disp('2. Quitter le programme ?');
    disp('               ');
    run=input('Reponse : ');
    disp('               ');
    while(isempty(run) || isstring(run) || (run~=1 && run~=2))
        run=input('Reponse : ');
        disp('               ');
    end 
end
