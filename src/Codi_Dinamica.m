function Dynamic_Triangular_Coupler_Online
% =========================================================================
% FOUR-BAR LINKAGE WITH TRIANGULAR COUPLER
% ANÁLISIS DINÁMICO CON MATLAB ONLINE
%
% Método:
% - Multiplicadores de Lagrange
% - Integración con ode45
% - Energía cinética, potencial y total
% - Fuerzas y momentos equivalentes en A, B y C
% - Errores numéricos de posición, velocidad, aceleración y energía
%
% Simplificación dinámica:
% El triángulo ABP se considera como una barra recta AB,
% manteniendo la masa del acoplador triangular.
% =========================================================================

clc; close all;

% -------------------------------------------------------------------------
% 1. PARÁMETROS GEOMÉTRICOS [m]
% -------------------------------------------------------------------------

L_O2O3 = 0.30;     % Ground link O2O3
L_O2A  = 0.10;     % Crank O2A
L_AB   = 0.25;     % Coupler AB, simplificación del triángulo ABP
L_BC   = 0.20;     % Secondary link BC
L_O3C  = 0.08;     % Pivot link O3C

% Hipótesis geométrica para cerrar el balancín compuesto
L_O3B  = 0.25;     % Effective compound rocker length O3B

% Pivotes fijos
O2 = [0; 0];
O3 = [L_O2O3; 0];

% -------------------------------------------------------------------------
% 2. PARÁMETROS FÍSICOS
% -------------------------------------------------------------------------

m_crank   = 1.5;   % Masa manivela O2A [kg]
m_coupler = 6.0;   % Masa acoplador ABP simplificado como AB [kg]
m_BC      = 3.0;   % Masa barra secundaria BC [kg]
m_O3C     = 1.2;   % Masa barra O3C [kg]

g = 9.81;          % Gravedad [m/s^2]

% -------------------------------------------------------------------------
% 3. PARÁMETROS TEMPORALES
% -------------------------------------------------------------------------

dt = 0.01;
tfin = 10;
tspan = 0:dt:tfin;

% -------------------------------------------------------------------------
% 4. MATRIZ DE MASAS CONSISTENTE
% -------------------------------------------------------------------------

% Coordenadas:
% q = [xA; yA; xB; yB; xC; yC; alpha]

M = zeros(7,7);

M_crank   = consistentMass(m_crank);
M_coupler = consistentMass(m_coupler);
M_BC      = consistentMass(m_BC);
M_O3C     = consistentMass(m_O3C);

% Manivela O2A: O2 fijo, A móvil
M(1:2,1:2) = M(1:2,1:2) + M_crank(3:4,3:4);

% Barra AB
M(1:4,1:4) = M(1:4,1:4) + M_coupler;

% Barra BC
M(3:6,3:6) = M(3:6,3:6) + M_BC;

% Barra O3C: O3 fijo, C móvil
M(5:6,5:6) = M(5:6,5:6) + M_O3C(1:2,1:2);

% -------------------------------------------------------------------------
% 5. VECTOR DE FUERZAS EXTERNAS
% -------------------------------------------------------------------------

Q = zeros(7,1);

% Peso de la manivela aplicado en A
Q(2) = Q(2) - 0.5*m_crank*g;

% Peso del acoplador AB aplicado entre A y B
Q(2) = Q(2) - 0.5*m_coupler*g;
Q(4) = Q(4) - 0.5*m_coupler*g;

% Peso de la barra BC aplicado entre B y C
Q(4) = Q(4) - 0.5*m_BC*g;
Q(6) = Q(6) - 0.5*m_BC*g;

% Peso de la barra O3C aplicado en C
Q(6) = Q(6) - 0.5*m_O3C*g;

% -------------------------------------------------------------------------
% 6. CONDICIONES INICIALES
% -------------------------------------------------------------------------

alpha0 = deg2rad(45);

A0 = O2 + L_O2A*[cos(alpha0); sin(alpha0)];

% B inicial: intersección entre círculo de centro A y círculo de centro O3
[B1, B2] = circleIntersections(A0, L_AB, O3, L_O3B);

if B1(2) >= B2(2)
    B0 = B1;
else
    B0 = B2;
end

% C inicial: intersección entre círculo de centro O3 y círculo de centro B
[C1, C2] = circleIntersections(O3, L_O3C, B0, L_BC);

if C1(1) <= C2(1)
    C0 = C1;
else
    C0 = C2;
end

q0 = [A0; B0; C0; alpha0];

% El sistema parte desde reposo
qp0 = zeros(7,1);

% Vector de estado inicial
y0 = [q0; qp0];

% -------------------------------------------------------------------------
% 7. INTEGRACIÓN NUMÉRICA
% -------------------------------------------------------------------------

options = odeset('RelTol',1e-7,'AbsTol',1e-9);

[tout, yout] = ode45(@deriv, tspan, y0, options);

% -------------------------------------------------------------------------
% 8. POSTPROCESADO
% -------------------------------------------------------------------------

n = length(tout);

Ekin_total = zeros(n,1);
Epot_total = zeros(n,1);
Etot       = zeros(n,1);

T_links = zeros(n,4);
V_links = zeros(n,4);
E_links = zeros(n,4);

err_pos = zeros(n,1);
err_vel = zeros(n,1);
err_acc = zeros(n,1);

force_A = zeros(n,3);
force_B = zeros(n,3);
force_C = zeros(n,3);

moment_A = zeros(n,1);
moment_B = zeros(n,1);
moment_C = zeros(n,1);

A_hist = zeros(n,2);
B_hist = zeros(n,2);
C_hist = zeros(n,2);

for i = 1:n

    q  = yout(i,1:7)';
    qp = yout(i,8:14)';

    A = q(1:2);
    B = q(3:4);
    C = q(5:6);

    A_hist(i,:) = A';
    B_hist(i,:) = B';
    C_hist(i,:) = C';

    [qdd, lambda, Phi, Phiq, gamma] = solveDynamics(q, qp);

    % ---------------------------------------------------------------------
    % Energías
    % ---------------------------------------------------------------------

    [Tlink, Vlink] = energyLinks(q, qp);

    T_links(i,:) = Tlink;
    V_links(i,:) = Vlink;
    E_links(i,:) = Tlink + Vlink;

    Ekin_total(i) = sum(Tlink);
    Epot_total(i) = sum(Vlink);
    Etot(i)       = Ekin_total(i) + Epot_total(i);

    % ---------------------------------------------------------------------
    % Errores numéricos
    % ---------------------------------------------------------------------

    err_pos(i) = norm(Phi);
    err_vel(i) = norm(Phiq*qp);
    err_acc(i) = norm(Phiq*qdd + gamma);

    % ---------------------------------------------------------------------
    % Fuerzas de restricción
    % ---------------------------------------------------------------------
    % En el sistema:
    % M*qdd + Phi_q'*lambda = Q
    % Por tanto, la fuerza de restricción sobre las coordenadas es:
    % R = -Phi_q'*lambda

    Rgen = -Phiq'*lambda;

    FA = Rgen(1:2);
    FB = Rgen(3:4);
    FC = Rgen(5:6);

    force_A(i,:) = [FA(1), FA(2), norm(FA)];
    force_B(i,:) = [FB(1), FB(2), norm(FB)];
    force_C(i,:) = [FC(1), FC(2), norm(FC)];

    % Momentos equivalentes producidos por las fuerzas de reacción.
    % En una articulación ideal de revolución no existe momento resistente puro.
    % Aquí se representa el momento equivalente de la fuerza respecto a una referencia.

    moment_A(i) = cross2D(A - O2, FA);
    moment_B(i) = cross2D(B - A,  FB);
    moment_C(i) = cross2D(C - O3, FC);

end

energy_error = abs(Etot - Etot(1));

% -------------------------------------------------------------------------
% 9. GRÁFICAS DE ENERGÍA TOTAL
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.65]);

plot(tout, Ekin_total, 'LineWidth', 1.5); hold on;
plot(tout, Epot_total, 'LineWidth', 1.5);
plot(tout, Etot, 'k--', 'LineWidth', 1.8);

grid on;
title('Total System Energy / Energía total del mecanismo');
xlabel('Time [s]');
ylabel('Energy [J]');
legend('Kinetic Energy','Potential Energy','Total Energy','Location','best');

% -------------------------------------------------------------------------
% 10. ENERGÍA POR ESLABÓN
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.75]);

subplot(3,1,1);
plot(tout, T_links(:,1), 'LineWidth', 1.3); hold on;
plot(tout, T_links(:,2), 'LineWidth', 1.3);
plot(tout, T_links(:,3), 'LineWidth', 1.3);
plot(tout, T_links(:,4), 'LineWidth', 1.3);
grid on;
title('Kinetic Energy per Link / Energía cinética por eslabón');
ylabel('T [J]');
legend('O2A','AB','BC','O3C','Location','best');

subplot(3,1,2);
plot(tout, V_links(:,1), 'LineWidth', 1.3); hold on;
plot(tout, V_links(:,2), 'LineWidth', 1.3);
plot(tout, V_links(:,3), 'LineWidth', 1.3);
plot(tout, V_links(:,4), 'LineWidth', 1.3);
grid on;
title('Potential Energy per Link / Energía potencial por eslabón');
ylabel('V [J]');
legend('O2A','AB','BC','O3C','Location','best');

subplot(3,1,3);
plot(tout, E_links(:,1), 'LineWidth', 1.3); hold on;
plot(tout, E_links(:,2), 'LineWidth', 1.3);
plot(tout, E_links(:,3), 'LineWidth', 1.3);
plot(tout, E_links(:,4), 'LineWidth', 1.3);
grid on;
title('Total Energy per Link / Energía total por eslabón');
xlabel('Time [s]');
ylabel('E [J]');
legend('O2A','AB','BC','O3C','Location','best');

% -------------------------------------------------------------------------
% 11. FUERZAS EN LAS ARTICULACIONES A, B Y C
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.8 0.75]);

subplot(3,1,1);
plot(tout, force_A(:,1), 'LineWidth', 1.2); hold on;
plot(tout, force_A(:,2), 'LineWidth', 1.2);
plot(tout, force_A(:,3), 'k--', 'LineWidth', 1.5);
grid on;
title('Reaction Force at A / Fuerza de reacción en A');
ylabel('Force [N]');
legend('F_A_x','F_A_y','|F_A|','Location','best');

subplot(3,1,2);
plot(tout, force_B(:,1), 'LineWidth', 1.2); hold on;
plot(tout, force_B(:,2), 'LineWidth', 1.2);
plot(tout, force_B(:,3), 'k--', 'LineWidth', 1.5);
grid on;
title('Reaction Force at B / Fuerza de reacción en B');
ylabel('Force [N]');
legend('F_B_x','F_B_y','|F_B|','Location','best');

subplot(3,1,3);
plot(tout, force_C(:,1), 'LineWidth', 1.2); hold on;
plot(tout, force_C(:,2), 'LineWidth', 1.2);
plot(tout, force_C(:,3), 'k--', 'LineWidth', 1.5);
grid on;
title('Reaction Force at C / Fuerza de reacción en C');
xlabel('Time [s]');
ylabel('Force [N]');
legend('F_C_x','F_C_y','|F_C|','Location','best');

% -------------------------------------------------------------------------
% 12. MOMENTOS EQUIVALENTES EN A, B Y C
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.15 0.15 0.75 0.6]);

plot(tout, moment_A, 'LineWidth', 1.5); hold on;
plot(tout, moment_B, 'LineWidth', 1.5);
plot(tout, moment_C, 'LineWidth', 1.5);

grid on;
title('Equivalent Moment at Hinges / Momento equivalente en articulaciones');
xlabel('Time [s]');
ylabel('Moment [N·m]');
legend('M_A','M_B','M_C','Location','best');

% -------------------------------------------------------------------------
% 13. ERRORES NUMÉRICOS
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.75]);

subplot(4,1,1);
plot(tout, energy_error, 'LineWidth', 1.5);
grid on;
title('Total Energy Error / Error de energía total');
ylabel('|E - E_0| [J]');

subplot(4,1,2);
plot(tout, err_pos, 'LineWidth', 1.5);
grid on;
title('Position Constraint Error / Error de posición');
ylabel('norm(\Phi)');

subplot(4,1,3);
plot(tout, err_vel, 'LineWidth', 1.5);
grid on;
title('Velocity Constraint Error / Error de velocidad');
ylabel('norm(\Phi_q qdot)');

subplot(4,1,4);
plot(tout, err_acc, 'LineWidth', 1.5);
grid on;
title('Acceleration Constraint Error / Error de aceleración');
ylabel('norm(error)');
xlabel('Time [s]');

% -------------------------------------------------------------------------
% 14. ANIMACIÓN DINÁMICA DEL MECANISMO
% -------------------------------------------------------------------------

figure('Color','w','Units','normalized','Position',[0.15 0.1 0.7 0.75]);

skip = 20;

for i = 1:skip:n

    A = A_hist(i,:)';
    B = B_hist(i,:)';
    C = C_hist(i,:)';

    clf;

    % Ground
    plot([O2(1) O3(1)], [O2(2) O3(2)], 'k--', 'LineWidth', 1.5);
    hold on;

    % Crank O2A
    plot([O2(1) A(1)], [O2(2) A(2)], '-or', ...
        'LineWidth', 3, 'MarkerFaceColor','r');

    % Coupler AB
    plot([A(1) B(1)], [A(2) B(2)], '-ob', ...
        'LineWidth', 3, 'MarkerFaceColor','b');

    % Secondary link BC
    plot([B(1) C(1)], [B(2) C(2)], '-og', ...
        'LineWidth', 3, 'MarkerFaceColor','g');

    % Pivot link CO3
    plot([C(1) O3(1)], [C(2) O3(2)], '-om', ...
        'LineWidth', 3, 'MarkerFaceColor','m');

    % Effective closing link O3B
    plot([O3(1) B(1)], [O3(2) B(2)], ':k', 'LineWidth', 1.2);

    % Labels
    text(O2(1)-0.015, O2(2)-0.015, 'O_2');
    text(O3(1)+0.005, O3(2)-0.015, 'O_3');
    text(A(1)+0.005, A(2)+0.005, 'A');
    text(B(1)+0.005, B(2)+0.005, 'B');
    text(C(1)+0.005, C(2)+0.005, 'C');

    axis equal;
    grid on;
    xlim([-0.15 0.50]);
    ylim([-0.25 0.35]);
    xlabel('x [m]');
    ylabel('y [m]');
    title(['Dynamic Motion | t = ', num2str(tout(i),'%.2f'), ' s']);

    drawnow;

end

disp('Dynamic simulation finished successfully.');

% =========================================================================
% FUNCIONES INTERNAS
% =========================================================================

    function yp = deriv(~, y)

        q  = y(1:7);
        qp = y(8:14);

        [qdd, ~, ~, ~, ~] = solveDynamics(q, qp);

        yp = [qp; qdd];

    end

% -------------------------------------------------------------------------

    function [qdd, lambda, Phi, Phiq, gamma] = solveDynamics(q, qp)

        [Phi, Phiq, gamma] = constraintData(q, qp);

        A_mat = [M, Phiq'; 
                 Phiq, zeros(6,6)];

        b_vec = [Q; -gamma];

        sol = A_mat \ b_vec;

        qdd = sol(1:7);
        lambda = sol(8:end);

    end

% -------------------------------------------------------------------------

    function [Phi, Phiq, gamma] = constraintData(q, qp)

        A = q(1:2);
        B = q(3:4);
        C = q(5:6);
        alpha = q(7);

        Ap = qp(1:2);
        Bp = qp(3:4);
        Cp = qp(5:6);
        alphap = qp(7);

        Phi = zeros(6,1);

        % Restricciones geométricas
        Phi(1) = (A-O2)'*(A-O2) - L_O2A^2;
        Phi(2) = (B-A)'*(B-A) - L_AB^2;
        Phi(3) = (B-O3)'*(B-O3) - L_O3B^2;
        Phi(4) = (C-B)'*(C-B) - L_BC^2;
        Phi(5) = (C-O3)'*(C-O3) - L_O3C^2;

        % Restricción de conducción con alpha
        % Se cambia la ecuación para evitar problemas numéricos cerca de 90°
        if abs(tan(alpha)) < 1
            Phi(6) = A(2) - O2(2) - L_O2A*sin(alpha);
        else
            Phi(6) = A(1) - O2(1) - L_O2A*cos(alpha);
        end

        % Jacobiana
        Phiq = zeros(6,7);

        % Phi 1: O2A
        Phiq(1,1:2) = 2*(A-O2)';

        % Phi 2: AB
        Phiq(2,1:2) = -2*(B-A)';
        Phiq(2,3:4) =  2*(B-A)';

        % Phi 3: O3B
        Phiq(3,3:4) = 2*(B-O3)';

        % Phi 4: BC
        Phiq(4,3:4) = -2*(C-B)';
        Phiq(4,5:6) =  2*(C-B)';

        % Phi 5: O3C
        Phiq(5,5:6) = 2*(C-O3)';

        % Phi 6: driving constraint
        if abs(tan(alpha)) < 1
            Phiq(6,2) = 1;
            Phiq(6,7) = -L_O2A*cos(alpha);
        else
            Phiq(6,1) = 1;
            Phiq(6,7) = L_O2A*sin(alpha);
        end

        % Término gamma
        gamma = zeros(6,1);

        gamma(1) = 2*(Ap'*Ap);
        gamma(2) = 2*((Bp-Ap)'*(Bp-Ap));
        gamma(3) = 2*(Bp'*Bp);
        gamma(4) = 2*((Cp-Bp)'*(Cp-Bp));
        gamma(5) = 2*(Cp'*Cp);

        if abs(tan(alpha)) < 1
            gamma(6) = L_O2A*(alphap^2)*sin(alpha);
        else
            gamma(6) = L_O2A*(alphap^2)*cos(alpha);
        end

    end

% -------------------------------------------------------------------------

    function M_e = consistentMass(m)

        M_e = m * [1/3 0   1/6 0;
                   0   1/3 0   1/6;
                   1/6 0   1/3 0;
                   0   1/6 0   1/3];

    end

% -------------------------------------------------------------------------

    function [Tlink, Vlink] = energyLinks(q, qp)

        A = q(1:2);
        B = q(3:4);
        C = q(5:6);

        Ap = qp(1:2);
        Bp = qp(3:4);
        Cp = qp(5:6);

        v0 = [0; 0];

        % Energía cinética por eslabón
        T_crank   = 0.5 * [v0; Ap]' * M_crank   * [v0; Ap];
        T_coupler = 0.5 * [Ap; Bp]' * M_coupler * [Ap; Bp];
        T_BC      = 0.5 * [Bp; Cp]' * M_BC      * [Bp; Cp];
        T_O3C     = 0.5 * [Cp; v0]' * M_O3C     * [Cp; v0];

        % Centros de gravedad
        G_crank   = 0.5*(O2 + A);
        G_coupler = 0.5*(A + B);
        G_BC      = 0.5*(B + C);
        G_O3C     = 0.5*(O3 + C);

        % Energía potencial por eslabón
        V_crank   = m_crank   * g * G_crank(2);
        V_coupler = m_coupler * g * G_coupler(2);
        V_BC      = m_BC      * g * G_BC(2);
        V_O3C     = m_O3C     * g * G_O3C(2);

        Tlink = [T_crank, T_coupler, T_BC, T_O3C];
        Vlink = [V_crank, V_coupler, V_BC, V_O3C];

    end

% -------------------------------------------------------------------------

    function [P1, P2] = circleIntersections(C0, r0, C1, r1)

        d = norm(C1 - C0);

        if d > r0 + r1
            error('No intersection: circles are too far apart.');
        elseif d < abs(r0 - r1)
            error('No intersection: one circle is inside the other.');
        elseif d == 0
            error('Infinite intersections: same circle center.');
        end

        a = (r0^2 - r1^2 + d^2)/(2*d);
        h = sqrt(max(r0^2 - a^2, 0));

        ex = (C1 - C0)/d;
        ey = [-ex(2); ex(1)];

        Pmid = C0 + a*ex;

        P1 = Pmid + h*ey;
        P2 = Pmid - h*ey;

    end

% -------------------------------------------------------------------------

    function mz = cross2D(r, f)

        mz = r(1)*f(2) - r(2)*f(1);

    end

end