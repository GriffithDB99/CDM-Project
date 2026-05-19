% =========================================================================
% FOUR-BAR LINKAGE WITH TRIANGULAR COUPLER
% ANÁLISIS CINEMÁTICO CON MATLAB ONLINE
%
% Mecanismo:
% O2A  = Manivela de entrada
% ABP  = Acoplador triangular
% BC   = Barra secundaria
% CO3  = Balancín/pivote de salida
%
% Resolución:
% - Posición mediante Newton-Raphson
% - Velocidad mediante matriz Jacobiana
% - Aceleración mediante matriz Jacobiana
% - Animación y gráficas
% =========================================================================

clc; clear; close all;

% -------------------------------------------------------------------------
% 1. PARÁMETROS GEOMÉTRICOS
% -------------------------------------------------------------------------

% Longitudes [cm]
L_O2O3 = 30;      % Ground link O2O3
L_O2A  = 10;      % Crank O2A
L_AB   = 25;      % Base del triángulo AB
L_AP   = 15;      % Lado AP del triángulo
L_BP   = 15;      % Lado BP del triángulo
L_BC   = 20;      % Barra secundaria BC
L_O3C  = 8;       % Barra pivotante O3C

% -------------------------------------------------------------------------
% NOTA AL PROFESOR:
% El enunciado da BC y O3C, pero no da el ángulo entre ellas.
% Para cerrar el mecanismo, se introduce una distancia efectiva O3B.
% Esto equivale a considerar que BC y O3C forman un balancín compuesto
% rígido. Sin este parametro no hemos conseguido que Matlab compute cxorrectamente.
% -------------------------------------------------------------------------

L_O3B = 25;       % [cm] Hipótesis geométrica para cerrar el mecanismo

% Pivotes fijos
O2 = [0; 0];
O3 = [L_O2O3; 0];

% -------------------------------------------------------------------------
% 2. PARÁMETROS TEMPORALES Y MOVIMIENTO DE ENTRADA
% -------------------------------------------------------------------------

dt = 0.01;          % [s] Paso de tiempo
t_end = 10;         % [s] Tiempo final
t_vec = 0:dt:t_end; % Para obtener solo una rotación en el tiempod definido
n = length(t_vec);

theta0 = 0;         % [rad] Ángulo inicial de la manivela
omega = 2*pi/t_end; % [rad/s] Una vuelta completa en t_end segundos
theta_ddot = 0;     % [rad/s^2] Velocidad angular constante

% Newton-Raphson
tol = 1e-10;
max_iter = 50;

% -------------------------------------------------------------------------
% 3. CÁLCULO DE LA POSICIÓN INICIAL
% -------------------------------------------------------------------------

theta = theta0;

% Punto A inicial
A0 = O2 + L_O2A*[cos(theta); sin(theta)];

% -------------------------------------------------------------------------
% Intersección de círculos para obtener B inicial
% Círculo 1: centro A0, radio AB
% Círculo 2: centro O3, radio O3B
% -------------------------------------------------------------------------

C0 = A0;
r0 = L_AB;
C1 = O3;
r1 = L_O3B;

d = norm(C1 - C0);
a = (r0^2 - r1^2 + d^2)/(2*d);
h = sqrt(max(r0^2 - a^2, 0));

ex = (C1 - C0)/d;
ey = [-ex(2); ex(1)];
Pmid = C0 + a*ex;

B_sol1 = Pmid + h*ey;
B_sol2 = Pmid - h*ey;

% Elegimos la solución superior
if B_sol1(2) >= B_sol2(2)
    B0 = B_sol1;
else
    B0 = B_sol2;
end

% -------------------------------------------------------------------------
% Intersección de círculos para obtener C inicial
% Círculo 1: centro O3, radio O3C
% Círculo 2: centro B0, radio BC
% -------------------------------------------------------------------------

C0 = O3;
r0 = L_O3C;
C1 = B0;
r1 = L_BC;

d = norm(C1 - C0);
a = (r0^2 - r1^2 + d^2)/(2*d);
h = sqrt(max(r0^2 - a^2, 0));

ex = (C1 - C0)/d;
ey = [-ex(2); ex(1)];
Pmid = C0 + a*ex;

C_sol1 = Pmid + h*ey;
C_sol2 = Pmid - h*ey;

% Elegimos la solución más cercana a la zona izquierda de O3
if C_sol1(1) <= C_sol2(1)
    C0_point = C_sol1;
else
    C0_point = C_sol2;
end

% -------------------------------------------------------------------------
% Intersección de círculos para obtener P inicial
% Círculo 1: centro A0, radio AP
% Círculo 2: centro B0, radio BP
% -------------------------------------------------------------------------

C0 = A0;
r0 = L_AP;
C1 = B0;
r1 = L_BP;

d = norm(C1 - C0);
a = (r0^2 - r1^2 + d^2)/(2*d);
h = sqrt(max(r0^2 - a^2, 0));

ex = (C1 - C0)/d;
ey = [-ex(2); ex(1)];
Pmid = C0 + a*ex;

P_sol1 = Pmid + h*ey;
P_sol2 = Pmid - h*ey;

% Elegimos la solución izquierda para que el triángulo sea visible
if P_sol1(1) <= P_sol2(1)
    P0 = P_sol1;
else
    P0 = P_sol2;
end

% Vector de coordenadas generalizadas
% q = [xA; yA; xB; yB; xC; yC; xP; yP]
q = [A0; B0; C0_point; P0];

% -------------------------------------------------------------------------
% 4. PREASIGNACIÓN DE HISTORIALES
% -------------------------------------------------------------------------

A_pos = zeros(n,2);    B_pos = zeros(n,2);
C_pos = zeros(n,2);    P_pos = zeros(n,2);
Gtri_pos = zeros(n,2); Gbc_pos = zeros(n,2);

A_vel = zeros(n,2);    B_vel = zeros(n,2);
C_vel = zeros(n,2);    P_vel = zeros(n,2);
Gtri_vel = zeros(n,2); Gbc_vel = zeros(n,2);

A_acc = zeros(n,2);    B_acc = zeros(n,2);
C_acc = zeros(n,2);    P_acc = zeros(n,2);
Gtri_acc = zeros(n,2); Gbc_acc = zeros(n,2);

err_pos = zeros(n,1);
err_vel = zeros(n,1);
err_acc = zeros(n,1);

% -------------------------------------------------------------------------
% 5. BUCLE PRINCIPAL DE SIMULACIÓN CINEMÁTICA
% -------------------------------------------------------------------------

for i = 1:n

    t = t_vec(i);

    % Movimiento impuesto de la manivela
    theta = theta0 + omega*t;
    theta_dot = omega;

    % ---------------------------------------------------------------------
    % 5.1 NEWTON-RAPHSON PARA POSICIÓN
    % ---------------------------------------------------------------------

    error_NR = 1;
    iter = 0;

    while error_NR > tol && iter < max_iter

        % Extraer coordenadas
        A = q(1:2);
        B = q(3:4);
        C = q(5:6);
        P = q(7:8);

        % Vector de restricciones Phi
        Phi = zeros(8,1);

        % Restricciones de la manivela O2A
        Phi(1) = A(1) - (O2(1) + L_O2A*cos(theta));
        Phi(2) = A(2) - (O2(2) + L_O2A*sin(theta));

        % Restricción de longitud AB
        Phi(3) = (B-A)'*(B-A) - L_AB^2;

        % Restricción de cierre del balancín compuesto O3B
        Phi(4) = (B-O3)'*(B-O3) - L_O3B^2;

        % Restricción de longitud BC
        Phi(5) = (C-B)'*(C-B) - L_BC^2;

        % Restricción de longitud O3C
        Phi(6) = (C-O3)'*(C-O3) - L_O3C^2;

        % Restricciones del triángulo AP y BP
        Phi(7) = (P-A)'*(P-A) - L_AP^2;
        Phi(8) = (P-B)'*(P-B) - L_BP^2;

        % Matriz Jacobiana Phi_q
        Phi_q = zeros(8,8);

        % Phi 1
        Phi_q(1,1) = 1;

        % Phi 2
        Phi_q(2,2) = 1;

        % Phi 3: AB
        Phi_q(3,1:2) = -2*(B-A)';
        Phi_q(3,3:4) =  2*(B-A)';

        % Phi 4: O3B
        Phi_q(4,3:4) = 2*(B-O3)';

        % Phi 5: BC
        Phi_q(5,3:4) = -2*(C-B)';
        Phi_q(5,5:6) =  2*(C-B)';

        % Phi 6: O3C
        Phi_q(6,5:6) = 2*(C-O3)';

        % Phi 7: AP
        Phi_q(7,1:2) = -2*(P-A)';
        Phi_q(7,7:8) =  2*(P-A)';

        % Phi 8: BP
        Phi_q(8,3:4) = -2*(P-B)';
        Phi_q(8,7:8) =  2*(P-B)';

        % Corrección de Newton-Raphson
        dq = -Phi_q \ Phi;
        q = q + dq;

        error_NR = norm(dq);
        iter = iter + 1;

    end

    % Recalcular restricciones con la posición ya convergida
    A = q(1:2);
    B = q(3:4);
    C = q(5:6);
    P = q(7:8);

    Phi = zeros(8,1);

    Phi(1) = A(1) - (O2(1) + L_O2A*cos(theta));
    Phi(2) = A(2) - (O2(2) + L_O2A*sin(theta));
    Phi(3) = (B-A)'*(B-A) - L_AB^2;
    Phi(4) = (B-O3)'*(B-O3) - L_O3B^2;
    Phi(5) = (C-B)'*(C-B) - L_BC^2;
    Phi(6) = (C-O3)'*(C-O3) - L_O3C^2;
    Phi(7) = (P-A)'*(P-A) - L_AP^2;
    Phi(8) = (P-B)'*(P-B) - L_BP^2;

    Phi_q = zeros(8,8);

    Phi_q(1,1) = 1;
    Phi_q(2,2) = 1;

    Phi_q(3,1:2) = -2*(B-A)';
    Phi_q(3,3:4) =  2*(B-A)';

    Phi_q(4,3:4) = 2*(B-O3)';

    Phi_q(5,3:4) = -2*(C-B)';
    Phi_q(5,5:6) =  2*(C-B)';

    Phi_q(6,5:6) = 2*(C-O3)';

    Phi_q(7,1:2) = -2*(P-A)';
    Phi_q(7,7:8) =  2*(P-A)';

    Phi_q(8,3:4) = -2*(P-B)';
    Phi_q(8,7:8) =  2*(P-B)';

    % ---------------------------------------------------------------------
    % 5.2 ANÁLISIS DE VELOCIDAD
    % Phi_q*q_dot = -Phi_t
    % ---------------------------------------------------------------------

    Phi_t = zeros(8,1);

    % Derivadas temporales de las restricciones de entrada
    Phi_t(1) = L_O2A * theta_dot * sin(theta);
    Phi_t(2) = -L_O2A * theta_dot * cos(theta);

    q_dot = Phi_q \ (-Phi_t);

    A_d = q_dot(1:2);
    B_d = q_dot(3:4);
    C_d = q_dot(5:6);
    P_d = q_dot(7:8);

    % ---------------------------------------------------------------------
    % 5.3 ANÁLISIS DE ACELERACIÓN
    % Phi_q*q_ddot = RHS_acc
    % ---------------------------------------------------------------------

    RHS_acc = zeros(8,1);

    % Restricciones de entrada
    RHS_acc(1) = -L_O2A*(theta_ddot*sin(theta) + theta_dot^2*cos(theta));
    RHS_acc(2) =  L_O2A*(theta_ddot*cos(theta) - theta_dot^2*sin(theta));

    % Términos centrípetos de las restricciones de distancia
    RHS_acc(3) = -2*((B_d - A_d)'*(B_d - A_d));
    RHS_acc(4) = -2*(B_d'*B_d);

    RHS_acc(5) = -2*((C_d - B_d)'*(C_d - B_d));
    RHS_acc(6) = -2*(C_d'*C_d);

    RHS_acc(7) = -2*((P_d - A_d)'*(P_d - A_d));
    RHS_acc(8) = -2*((P_d - B_d)'*(P_d - B_d));

    q_ddot = Phi_q \ RHS_acc;

    A_dd = q_ddot(1:2);
    B_dd = q_ddot(3:4);
    C_dd = q_ddot(5:6);
    P_dd = q_ddot(7:8);

    % ---------------------------------------------------------------------
    % 5.4 CENTROS DE GRAVEDAD
    % ---------------------------------------------------------------------

    Gtri = (A + B + P)/3;     % Centroide del triángulo ABP
    Gbc = (B + C)/2;          % Centro de gravedad de BC

    Gtri_d = (A_d + B_d + P_d)/3;
    Gbc_d = (B_d + C_d)/2;

    Gtri_dd = (A_dd + B_dd + P_dd)/3;
    Gbc_dd = (B_dd + C_dd)/2;

    % ---------------------------------------------------------------------
    % 5.5 GUARDAR RESULTADOS
    % ---------------------------------------------------------------------

    A_pos(i,:) = A';
    B_pos(i,:) = B';
    C_pos(i,:) = C';
    P_pos(i,:) = P';

    Gtri_pos(i,:) = Gtri';
    Gbc_pos(i,:) = Gbc';

    A_vel(i,:) = A_d';
    B_vel(i,:) = B_d';
    C_vel(i,:) = C_d';
    P_vel(i,:) = P_d';

    Gtri_vel(i,:) = Gtri_d';
    Gbc_vel(i,:) = Gbc_d';

    A_acc(i,:) = A_dd';
    B_acc(i,:) = B_dd';
    C_acc(i,:) = C_dd';
    P_acc(i,:) = P_dd';

    Gtri_acc(i,:) = Gtri_dd';
    Gbc_acc(i,:) = Gbc_dd';

    % Errores numéricos
    err_pos(i) = norm(Phi);
    err_vel(i) = norm(Phi_q*q_dot + Phi_t);
    err_acc(i) = norm(Phi_q*q_ddot - RHS_acc);

end

% -------------------------------------------------------------------------
% 6. MAGNITUDES PARA REPRESENTACIÓN GRÁFICA
% -------------------------------------------------------------------------

% Desplazamiento lineal respecto a la posición inicial
A_disp = sqrt(sum((A_pos - A_pos(1,:)).^2, 2));
B_disp = sqrt(sum((B_pos - B_pos(1,:)).^2, 2));
C_disp = sqrt(sum((C_pos - C_pos(1,:)).^2, 2));
P_disp = sqrt(sum((P_pos - P_pos(1,:)).^2, 2));
Gtri_disp = sqrt(sum((Gtri_pos - Gtri_pos(1,:)).^2, 2));
Gbc_disp = sqrt(sum((Gbc_pos - Gbc_pos(1,:)).^2, 2));

% Velocidades
A_vel_mag = sqrt(sum(A_vel.^2, 2));
B_vel_mag = sqrt(sum(B_vel.^2, 2));
C_vel_mag = sqrt(sum(C_vel.^2, 2));
P_vel_mag = sqrt(sum(P_vel.^2, 2));
Gtri_vel_mag = sqrt(sum(Gtri_vel.^2, 2));
Gbc_vel_mag = sqrt(sum(Gbc_vel.^2, 2));

% Aceleraciones
A_acc_mag = sqrt(sum(A_acc.^2, 2));
B_acc_mag = sqrt(sum(B_acc.^2, 2));
C_acc_mag = sqrt(sum(C_acc.^2, 2));
P_acc_mag = sqrt(sum(P_acc.^2, 2));
Gtri_acc_mag = sqrt(sum(Gtri_acc.^2, 2));
Gbc_acc_mag = sqrt(sum(Gbc_acc.^2, 2));

% -------------------------------------------------------------------------
% 7. GRÁFICAS DE DESPLAZAMIENTO, VELOCIDAD Y ACELERACIÓN
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.75]);

subplot(3,1,1);
plot(t_vec, A_disp, 'LineWidth', 1.3); hold on;
plot(t_vec, B_disp, 'LineWidth', 1.3);
plot(t_vec, C_disp, 'LineWidth', 1.3);
plot(t_vec, P_disp, 'LineWidth', 1.3);
plot(t_vec, Gtri_disp, '--', 'LineWidth', 1.3);
plot(t_vec, Gbc_disp, '--', 'LineWidth', 1.3);
grid on;
title('Linear Displacement / Desplazamiento lineal');
ylabel('Displacement [cm]');
legend('A','B','C','P','CG ABP','CG BC','Location','best');

subplot(3,1,2);
plot(t_vec, A_vel_mag, 'LineWidth', 1.3); hold on;
plot(t_vec, B_vel_mag, 'LineWidth', 1.3);
plot(t_vec, C_vel_mag, 'LineWidth', 1.3);
plot(t_vec, P_vel_mag, 'LineWidth', 1.3);
plot(t_vec, Gtri_vel_mag, '--', 'LineWidth', 1.3);
plot(t_vec, Gbc_vel_mag, '--', 'LineWidth', 1.3);
grid on;
title('Velocity Magnitude / Magnitud de velocidad');
ylabel('Velocity [cm/s]');
legend('A','B','C','P','CG ABP','CG BC','Location','best');

subplot(3,1,3);
plot(t_vec, A_acc_mag, 'LineWidth', 1.3); hold on;
plot(t_vec, B_acc_mag, 'LineWidth', 1.3);
plot(t_vec, C_acc_mag, 'LineWidth', 1.3);
plot(t_vec, P_acc_mag, 'LineWidth', 1.3);
plot(t_vec, Gtri_acc_mag, '--', 'LineWidth', 1.3);
plot(t_vec, Gbc_acc_mag, '--', 'LineWidth', 1.3);
grid on;
title('Acceleration Magnitude / Magnitud de aceleración');
ylabel('Acceleration [cm/s^2]');
xlabel('Time [s]');
legend('A','B','C','P','CG ABP','CG BC','Location','best');

% -------------------------------------------------------------------------
% 8. GRÁFICAS DE ERRORES NUMÉRICOS
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.15 0.15 0.7 0.6]);

subplot(3,1,1);
plot(t_vec, err_pos, 'LineWidth', 1.5);
grid on;
title('Position Constraint Error / Error de posición');
ylabel('norm(\Phi)');

subplot(3,1,2);
plot(t_vec, err_vel, 'LineWidth', 1.5);
grid on;
title('Velocity Constraint Error / Error de velocidad');
ylabel('norm(error)');

subplot(3,1,3);
plot(t_vec, err_acc, 'LineWidth', 1.5);
grid on;
title('Acceleration Constraint Error / Error de aceleración');
ylabel('norm(error)');
xlabel('Time [s]');

% -------------------------------------------------------------------------
% 9. TRAYECTORIA DEL PUNTO P
% -------------------------------------------------------------------------

figure('Color','w');
plot(P_pos(:,1), P_pos(:,2), 'r', 'LineWidth', 1.8); hold on;
plot(P_pos(1,1), P_pos(1,2), 'go', 'MarkerFaceColor','g');
plot(P_pos(end,1), P_pos(end,2), 'ko', 'MarkerFaceColor','k');
grid on; axis equal;
title('Trajectory of Coupler Point P / Trayectoria del punto P');
xlabel('x [cm]');
ylabel('y [cm]');
legend('Trajectory P','Initial P','Final P','Location','best');

% -------------------------------------------------------------------------
% 10. ANIMACIÓN DEL MECANISMO
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.15 0.1 0.7 0.75]);

skip = 8;   % Para que MATLAB Online no vaya demasiado lento

for i = 1:skip:n

    A = A_pos(i,:)';
    B = B_pos(i,:)';
    C = C_pos(i,:)';
    P = P_pos(i,:)';
    Gtri = Gtri_pos(i,:)';
    Gbc = Gbc_pos(i,:)';

    clf;

    % Base fija O2O3
    plot([O2(1) O3(1)], [O2(2) O3(2)], 'k--', 'LineWidth', 1.5); 
    hold on;

    % Manivela O2A
    plot([O2(1) A(1)], [O2(2) A(2)], '-or', ...
        'LineWidth', 3, 'MarkerFaceColor','r');

    % Acoplador triangular ABP
    patch([A(1) B(1) P(1)], [A(2) B(2) P(2)], [0.8 0.8 0.8], ...
        'FaceAlpha', 0.65, 'EdgeColor', 'k', 'LineWidth', 1.5);

    % Barra secundaria BC
    plot([B(1) C(1)], [B(2) C(2)], '-ob', ...
        'LineWidth', 3, 'MarkerFaceColor','b');

    % Barra pivotante CO3
    plot([C(1) O3(1)], [C(2) O3(2)], '-om', ...
        'LineWidth', 3, 'MarkerFaceColor','m');

    % Punto P
    plot(P(1), P(2), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 8);

    % Centros de gravedad
    plot(Gtri(1), Gtri(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(Gbc(1), Gbc(2), 'ks', 'MarkerSize', 8, 'LineWidth', 1.5);

    % Trayectoria de P
    plot(P_pos(1:i,1), P_pos(1:i,2), 'r:', 'LineWidth', 1.2);

    % Etiquetas
    text(O2(1)-2, O2(2)-2, 'O_2');
    text(O3(1)+1, O3(2)-2, 'O_3');
    text(A(1)+0.7, A(2)+0.7, 'A');
    text(B(1)+0.7, B(2)+0.7, 'B');
    text(C(1)+0.7, C(2)+0.7, 'C');
    text(P(1)+0.7, P(2)+0.7, 'P');

    axis equal;
    grid on;
    xlim([-20 55]);
    ylim([-35 45]);
    xlabel('x [cm]');
    ylabel('y [cm]');

    title(['Four-Bar Linkage with Triangular Coupler | t = ', ...
        num2str(t_vec(i),'%.2f'), ' s']);

    drawnow;

end
