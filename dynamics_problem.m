% This code is used to solve under carriage problem
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% 22 - 5 - 2016
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
format shortG
%% DATA
        %% Beam-12
         theta(1)=90; % degree
         Length(1)=1; % m
         Area(1)=0.01; % m^2
         Inertia(1)=0.00001; % m^4
         mu(1)=50; % kg/unit length
         E(1)=2.06e11; % N/m^2
        %% Beam-23
         theta(2)=90; % degree
         Length(2)=1; % m
         Area(2)=0.01; % m^2
         Inertia(2)=0.00001; % m^4
         mu(2)=50; % kg/unit length
         E(2)=2.06e11; % N/m^2
        %% Beam-24
         theta(3)=90+45; % degree
         Length(3)=1*sqrt(2); % m
         Area(3)=0.01; % m^2
         Inertia(3)=0.00001; % m^4
         mu(3)=50; % kg/unit length
         E(3)=2.06e11; % N/m^2
        %% Damping parameter
        % C = alpha * M + beta * K
        alpha=1;
        beta=0.00005;
        %% Forces
        % Force=[f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12]
        Force=[1000, 5000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        %% Response Time
        Time=linspace(0,0.025,10000);
        %% Plotting Control
        % Plot = 0 --> No Plotting data    Plot = 1 --> Plotting data
        Plot=1;  
        %% Save figures control
        % Save = 0 --> No Saving figures    Save = 1 --> Saving figures
        % The figures will saved in running folder directory
        Save=1;
%% Mass Matrix
        %% global mass matrix of each beam
        for i=1:length(Area)
            Cx=cosd(theta(i));
            Cy=sind(theta(i));
            M=mu(i)*Length(i);
            ML=M*Length(i);
            ML2=ML*Length(i);
            m{i}=[M*Cx^2/3+13/35*M*Cy^2,        -4/105*M*Cx*Cy,                         -11/210*ML*Cy,       M*Cx^2/6+9/70*M*Cy^2,           4/105*M*Cx*Cy,                           13/420*ML*Cy ; ...
                      -4/105*M*Cx*Cy,                         M*Cy^2/3+13/35*M*Cx^2,         11/210*ML*Cx,        4/105*M*Cx*Cy,                            M*Cy^2/6+9/70*M*Cx^2,          -13/420*ML*Cx ; ...
                      -11/210*ML*Cy,                           11/210*ML*Cx,                             ML2/105,                  -13/420*ML*Cy,                            13/420*ML*Cx,                             -1/140*ML2 ; ...
                      M*Cx^2/6+9/70*M*Cy^2,          4/105*M*Cx*Cy,                           -13/420*ML*Cy,        M*Cx^2/3+13/35*M*Cy^2,        -4/105*M*Cx*Cy,                          11/210*ML*Cy ; ...
                      4/105*M*Cx*Cy,                           M*Cy^2/6+9/70*M*Cx^2,          13/420*ML*Cx,          -4/105*M*Cx*Cy,                         M*Cy^2/3+13/35*M*Cx^2,         -11/210*ML*Cx ; ...
                      13/420*ML*Cy,                             -13/420*ML*Cx,                           -1/140*ML2,              11/210*ML*Cy,                             -11/210*ML*Cx,                           ML2/105 ];
        end
        %% global stiffness matrix of under carriage
        clear M;
        M=zeros(12,12);
        M(1:6,1:3)=m{1}(1:6,1:3);
        M(1:3,4:6)=m{1}(1:3,4:6);
        M(7:9,4:9)=m{2}(4:6,1:6);
        M(4:6,7:9)=m{2}(1:3,4:6);
        M(10:12,4:6)=m{3}(4:6,1:3);
        M(10:12,10:12)=m{3}(4:6,4:6);
        M(4:6,10:12)=m{3}(1:3,4:6);
        M(4:6,4:6)=m{1}(4:6,4:6)+m{2}(1:3,1:3)+m{3}(1:3,1:3);
        %% mass matrix of free nodes (1,2,3,4,5,6)
        M_free=M(1:6,1:6);
        %% reduced mass matrix nodes (1,2,4,5)
        AM=[M_free(3,1:2),M_free(3,4:5);M_free(6,1:2),M_free(6,4:5)];
        BM=[M_free(3,3),M_free(3,6);M_free(6,3),M_free(6,6)];
        CM=[M_free(1:2,1:2),M_free(1:2,4:5);M_free(4:5,1:2),M_free(4:5,4:5)];
        DM=[M_free(1:2,3),M_free(1:2,6);M_free(4:5,3),M_free(4:5,6)];
        M_reduced=CM-DM*inv(BM)*AM;
%% Stiffness Matrix
        %% global stiffness matrix of each beam
        for i=1:length(Area)
            cx=cosd(theta(i));
            cy=sind(theta(i));
            A=E(i)*Area(i)/Length(i);
            B=E(i)*Inertia(i)/Length(i);
            C=B/Length(i);
            D=C/Length(i);
            k{i}=[A*cx^2+12*D*cy^2,     (A-12*D)*cx*cy,                 -6*C*cy,          -A*cx^2-12*D*cy^2,           (-A+12*D)*cx*cy,                -6*C*cy ; ...
                     (A-12*D)*cx*cy,             A*cy^2+12*D*cx^2,         6*C*cx,            (-A+12*D)*cx*cy,                -A*cy^2-12*D*cx^2            6*C*cx ; ...
                     -6*C*cy,                         6*C*cx,                               4*B,                  6*C*cy,                                -6*C*cx,                               2*B ; ...
                     -A*cx^2-12*D*cy^2,    (-A+12*D)*cx*cy,               6*C*cy,            A*cx^2+12*C*cy^2,           (A-12*D)*cx*cy,                   6*C*cy ; ...
                     (-A+12*D)*cx*cy,          -A*cy^2-12*D*cx^2          -6*C*cx,          (A-12*D)*cx*cy,                   A*cy^2+12*D*cx^2,           -6*C*cx ; ...
                     -6*C*cy,                         6*C*cx,                               2*B,                 6*C*cy,                                 -6*C*cx,                               4*B];
        end
        %% global stiffness matrix of under carriage
        K=zeros(12,12);
        K(1:6,1:3)=k{1}(1:6,1:3);
        K(1:3,4:6)=k{1}(1:3,4:6);
        K(7:9,4:9)=k{2}(4:6,1:6);
        K(4:6,7:9)=k{2}(1:3,4:6);
        K(10:12,4:6)=k{3}(4:6,1:3);
        K(10:12,10:12)=k{3}(4:6,4:6);
        K(4:6,10:12)=k{3}(1:3,4:6);
        K(4:6,4:6)=k{1}(4:6,4:6)+k{2}(1:3,1:3)+k{3}(1:3,1:3);
        %% siffness matrix of free nodes (1,2,3,4,5,6)
        K_free=K(1:6,1:6);
        %% reduced siffness matrix nodes (1,2,4,5)
        AK=[K_free(3,1:2),K_free(3,4:5);K_free(6,1:2),K_free(6,4:5)];
        BK=[K_free(3,3),K_free(3,6);K_free(6,3),K_free(6,6)];
        CK=[K_free(1:2,1:2),K_free(1:2,4:5);K_free(4:5,1:2),K_free(4:5,4:5)];
        DK=[K_free(1:2,3),K_free(1:2,6);K_free(4:5,3),K_free(4:5,6)];
        K_reduced=CK-DK*inv(BK)*AK;
%% Damping Matrix
        %% global damping matrix of under carriage
        C=alpha*M+beta*K;
        %% damping matrix of free nodes (1,2,3,4,5,6)
        C_free=alpha*M_free+beta*K_free;
        %% reduced damping matrix
        C_reduced=alpha*M_reduced+beta*K_reduced;
%% Forces
        %% forces vector of free nodes (1,2,3,4,5,6)
        Force_free=Force(1:6);
        %% reduced forces vector
        Force_reduced=[Force(1:2),Force(4:5)]';
%% EigenValues (lamda) EigenVectors
[Eigen_Vectors, Lamda]=eig(inv(K_reduced)*M_reduced);
%% Natural Frequancy
Wn=diag(sqrt(inv(Lamda)));
%% Orthonormal Eigen Vectors
for i=1:length(Lamda(1,:))
    Alpha(i)=sqrt(1./(Eigen_Vectors(:,i)'*M_reduced*Eigen_Vectors(:,i)));
    Orthonormal_Eigen_Vectors(:,i)=Alpha(i)*Eigen_Vectors(:,i);
end
%% Transformation Axes (q)
Q=Orthonormal_Eigen_Vectors'*Force_reduced;
M_bar=Orthonormal_Eigen_Vectors'*M_reduced*Orthonormal_Eigen_Vectors;
K_bar=Orthonormal_Eigen_Vectors'*K_reduced*Orthonormal_Eigen_Vectors;
C_bar=Orthonormal_Eigen_Vectors'*C_reduced*Orthonormal_Eigen_Vectors;
%% Damping Coefficient (Zeta)
Zeta=diag(C_bar)./Wn/2;
%% Damping Frequancy (Wd)
Wd=Wn.*sqrt(1-Zeta.^2);
%% Displacement in Transformation Axes (q)
for i=1:length(Zeta)
    q(i,:)=Q(i)./(M_bar(i,i)*Wn(i)^2)*(1-exp(-Zeta(i)*Wn(i)*Time).*(cos(Wd(i)*Time)+Zeta(i)./sqrt(1-Zeta(i)^2).*sin(Wd(i)*Time)));
end
%% Velocity in Transformation Axes (q_dot)
for i=1:length(Zeta)
    q_dot(i,:)=Q(i)./(M_bar(i,i)*Wn(i)^2)*(Zeta(i)*Wn(i)*exp(-Zeta(i)*Wn(i)*Time).*(cos(Wd(i)*Time)+Zeta(i)./sqrt(1-Zeta(i)^2).*sin(Wd(i)*Time)) ...
                                                                   +exp(-Zeta(i)*Wn(i)*Time).*(Wd(i)*sin(Wd(i)*Time)-Wd(i)*Zeta(i)./sqrt(1-Zeta(i)^2).*cos(Wd(i)*Time)));
end
%% Acceleration in Transformation Axes (q_dot2)
for i=1:length(Zeta)
    q_dot2(i,:)=Q(i)./(M_bar(i,i)*Wn(i)^2)*(-exp(-Zeta(i)*Wn(i)*Time)*Zeta(i)^2*Wn(i)^2.*(cos(Wd(i)*Time)+Zeta(i)/sqrt(1-Zeta(i)^2)*sin(Wd(i)*Time)) ...
                                                                 +2*exp(-Zeta(i)*Wn(i)*Time)*Zeta(i)*Wn(i).*(Zeta(i)*Wd(i)/sqrt(1-Zeta(i)^2)*cos(Wd(i)*Time)-Wd(i)*sin(Wd(i)*Time)) ...
                                                                 +exp(-Zeta(i)*Wn(i)*Time).*(Wd(i)^2*cos(Wd(i)*Time)+Zeta(i)*Wd(i)^2/sqrt(1-Zeta(i)^2)*sin(Wd(i)*Time)));
end
%% Displacement in Initial Axes (Z)
Z=Orthonormal_Eigen_Vectors*q;
%% Velocity in Initial Axes (Z_dot)
Z_dot=Orthonormal_Eigen_Vectors*q_dot;
%% Acceleration in Initial Axes (Z_dot2)
Z_dot2=Orthonormal_Eigen_Vectors*q_dot2;
%% Plotting
if Plot ==1
        %% EigenVectors
        set(0,'defaultfigureposition',[0 50 1700 630]);
        figure('Name','Eigen Vectors','NumberTitle','off');
        set(gcf,'color','w')
        hold all;
        for i=1:length(diag(Lamda))
            plot(1:length(diag(Lamda)),Eigen_Vectors(:,i),'linewidth',2)
            LegenD{i}=['@ \Lambda = ' num2str(Lamda(i,i))];
        end
        grid on;
        title('Eigen Vectors','fontsize',18)
        ylabel('Eigen Vector Values','fontsize',18)
        xlabel('Length','fontsize',18)
        legend(LegenD)
        %% Displacement in Transformation Axes
        figure('Name','Displacement in Transformation Axes','NumberTitle','off');
        set(gcf,'color','w')
        subplot(4,1,1)
        plot(Time,q(1,:),'linewidth',2)
        title('Displacement in Transformation Axes','fontsize',18)
        grid on;
        ylabel('q_1(t) (m)','fontsize',18)
        subplot(4,1,2)
        plot(Time,q(2,:),'linewidth',2)
        grid on;
        ylabel('q_2(t) (m)','fontsize',18)
        subplot(4,1,3)
        plot(Time,q(3,:),'linewidth',2)
        grid on;
        ylabel('q_3(t) (m)','fontsize',18)
        subplot(4,1,4)
        plot(Time,q(4,:),'linewidth',2)
        grid on;
        ylabel('q_4(t) (m)','fontsize',18)
        xlabel('Time (sec)','fontsize',18)
        %% Velocity in Transformation Axes
        set(0,'defaultfigureposition',[0 50 1700 630]);
        figure('Name','Velocity in Transformation Axes','NumberTitle','off');
        set(gcf,'color','w')
        subplot(4,1,1)
        plot(Time,q_dot(1,:),'linewidth',2)
        title('Velocity in Transformation Axes','fontsize',18)
        grid on;
        ylabel('q_1^o(t) (m)','fontsize',18)
        subplot(4,1,2)
        plot(Time,q_dot(2,:),'linewidth',2)
        grid on;
        ylabel('q_2^o(t) (m)','fontsize',18)
        subplot(4,1,3)
        plot(Time,q_dot(3,:),'linewidth',2)
        grid on;
        ylabel('q_3^o(t) (m)','fontsize',18)
        subplot(4,1,4)
        plot(Time,q_dot(4,:),'linewidth',2)
        grid on;
        ylabel('q_4^o(t) (m)','fontsize',18)
        xlabel('Time (sec)','fontsize',18)
        %% Acceletration in Transformation Axes
        set(0,'defaultfigureposition',[0 50 1700 630]);
        figure('Name','Acceletration in Transformation Axes','NumberTitle','off');
        set(gcf,'color','w')
        subplot(4,1,1)
        plot(Time,q_dot(1,:),'linewidth',2)
        title('Acceletration in Transformation Axes','fontsize',18)
        grid on;
        ylabel('q_1^o^o(t) (m)','fontsize',18)
        subplot(4,1,2)
        plot(Time,q_dot2(2,:),'linewidth',2)
        grid on;
        ylabel('q_2^o^o(t) (m)','fontsize',18)
        subplot(4,1,3)
        plot(Time,q_dot2(3,:),'linewidth',2)
        grid on;
        ylabel('q_3^o^o(t) (m)','fontsize',18)
        subplot(4,1,4)
        plot(Time,q_dot2(4,:),'linewidth',2)
        grid on;
        ylabel('q_4^o^o(t) (m)','fontsize',18)
        xlabel('Time (sec)','fontsize',18)
        %% Displacement in Real Axes
        figure('Name','Displacement in Real Axes','NumberTitle','off');
        set(gcf,'color','w')
        subplot(4,1,1)
        plot(Time,Z(1,:),'linewidth',2)
        title('Displacement in Real Axes','fontsize',18)
        grid on;
        ylabel('u_1(t) (m)','fontsize',18)
        subplot(4,1,2)
        plot(Time,Z(2,:),'linewidth',2)
        grid on;
        ylabel('v_1(t) (m)','fontsize',18)
        subplot(4,1,3)
        plot(Time,Z(3,:),'linewidth',2)
        grid on;
        ylabel('u_2(t) (m)','fontsize',18)
        subplot(4,1,4)
        plot(Time,Z(4,:),'linewidth',2)
        grid on;
        ylabel('v_2(t) (m)','fontsize',18)
        xlabel('Time (sec)','fontsize',18)
        %% Velocity in Real Axes
        figure('Name','Velocity in Real Axes','NumberTitle','off');
        set(gcf,'color','w')
        subplot(4,1,1)
        plot(Time(1:length(Z_dot(1,:))),Z_dot(1,:),'linewidth',2)
        title('Velocity in Real Axes','fontsize',18)
        grid on;
        ylabel('u_1^o(t) (m)','fontsize',18)
        subplot(4,1,2)
        plot(Time(1:length(Z_dot(2,:))),Z_dot(2,:),'linewidth',2)
        grid on;
        ylabel('v_1^o(t) (m)','fontsize',18)
        subplot(4,1,3)
        plot(Time(1:length(Z_dot(3,:))),Z_dot(3,:),'linewidth',2)
        grid on;
        ylabel('u_2^o(t) (m)','fontsize',18)
        subplot(4,1,4)
        plot(Time(1:length(Z_dot(4,:))),Z_dot(4,:),'linewidth',2)
        grid on;
        ylabel('v_2^o(t) (m)','fontsize',18)
        xlabel('Time (sec)','fontsize',18)
        %% Acceletration in Real Axes
        figure('Name','Acceletration in Real Axes','NumberTitle','off');
        set(gcf,'color','w')
        subplot(4,1,1)
        plot(Time(1:length(Z_dot2(1,:))),Z_dot2(1,:),'linewidth',2)
        title('Acceletration in Real Axes','fontsize',18)
        grid on;
        ylabel('u_1^o^o(t) (m)','fontsize',18)
        subplot(4,1,2)
        plot(Time(1:length(Z_dot2(2,:))),Z_dot2(2,:),'linewidth',2)
        grid on;
        ylabel('v_1^o^o(t) (m)','fontsize',18)
        subplot(4,1,3)
        plot(Time(1:length(Z_dot2(3,:))),Z_dot2(3,:),'linewidth',2)
        grid on;
        ylabel('u_2^o^o(t) (m)','fontsize',18)
        subplot(4,1,4)
        plot(Time(1:length(Z_dot2(4,:))),Z_dot2(4,:),'linewidth',2)
        grid on;
        ylabel('v_2^o^o(t) (m)','fontsize',18)
        xlabel('Time (sec)','fontsize',18)
end
%% Save Figures
if Save ==1 && Plot ==1
    for S=1:7
        figure(S);
        saveas(gcf, [num2str(S) '.png']);
    end
end