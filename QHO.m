%% Class definition for Quantum Harmonic Oscilator (QHO)
% Implements a class for the QHO
% Solves the equation and plots the solution

%% Declare class variables
%

classdef QHO

    properties
        
        % Fundamental constants
        hbar = 1.054571e-34;                        % Plank's constant, J*s
        m_e = 1.62661e-27;                          % mass of an electron, Kg
        
        % Physical dimensions
        L = 2e-10;                                  % length of domain, m
        omega = 5.63212e14;                         % frequency of the potential, s^-1
        x_scale;                                    % characteristic length of the system
        
        % Other choices
        n_elect = 1;                                % Number of electrons in the system
        max_freq = 6;                              % maximum frequency of a plane wave in the basis set, 2 * n_PW + 1 will be used
        
        % CPU information for the solution
        fft_CPU = 0;                                % CPU time required to compute the fast Fourier transform (FFT) of the potential, s
        eig_CPU = 0;                                % CPU time required to compute the eigenvalues of the Hamiltonian matrix, s
        total_CPU = 0;
        
        % Solution information
        eig_vecs;
        eig_vals;
        
    end
    
    methods
        
        %% Initialize the class
        % constructor does nothing
        
        function qho = QHO(nmf)
            qho.x_scale = sqrt(qho.hbar / qho.m_e / qho.omega);            % used to non-dimensionalize, m
            qho.max_freq = nmf;
        end
        
        %% Solve the equation with eigenvalues
        % 
        function qho = solve(qho)
            
            before_solve = cputime;
            
            % Build Potential
            n_bais_vecs = 2 * qho.max_freq + 1;                                     % N + 1 basis functions are used
            n_fourier = 4 * qho.max_freq + 1;
            n_eng_levels = ceil(qho.n_elect / 2);                                    % number of energy levels occupied by electrons in the ground state
            x = linspace(-qho.L/2, qho.L/2 - qho.L/n_fourier, n_fourier);           % spatially resolve the domain
            V = qho.harm_pot(x);
            
            % Take fft of the potential
            before_fft = cputime;
            pot_freq = fft(V);
            pot_freq = circshift(pot_freq,[1, 2 * qho.max_freq]);
            qho.fft_CPU = cputime - before_fft;
            
            %% Build Hamiltonian Matrix
            
            PW_freqs = -qho.max_freq : qho.max_freq;
            
            % Build kinetic energy matrix
            Ham_KE = zeros(n_bais_vecs, n_bais_vecs);
            for k = 1:n_bais_vecs
                Ham_KE(k,k) = qho.hbar ^ 2 / 2 / qho.m_e * qho.L ^ -2 * 4 * pi^2 * PW_freqs(k)^2;
            end
            
            % Build potential energy matrix
            Ham_PE = zeros(n_bais_vecs, n_bais_vecs);
            for i = 1:n_bais_vecs
                for j = 1:n_bais_vecs
                    freqdiff = PW_freqs(i) - PW_freqs(j);
                    Ham_PE(i,j) = pot_freq( freqdiff + n_bais_vecs) / n_fourier;
                end
            end
            
            Ham = Ham_KE + Ham_PE;
            %Ham = real(Ham);
            
            % Solve the eigenvalues
            before_eig = cputime;
            [qho.eig_vecs, qho.eig_vals] = eig(Ham);       % returns all eigenvalues and eigenvectors
            % [Vecs, Vals] = eigs(Ham,n_eng_levels,'SM');       % This does not work...
            qho.eig_CPU = cputime - before_eig;
            
            qho.total_CPU = cputime - before_solve;
            
        end
        
        %% Plot the potential
        function qho = plot_pot(qho)

            n_pp = 100;                                % number of points to use in the plot
            x_vec = linspace(-qho.L/2, qho.L/2, n_pp);
            V_vec = qho.harm_pot(x_vec);
            
            set(0,'defaultlinelinewidth',1.5)
            set(0,'defaultaxeslinewidth',2)

            %clf
            plot(x_vec, V_vec)
            xlabel('Position (m)')
            ylabel('Potential (J)')
            ax = gca;
            ax.FontSize = 20;
        end
        
        %% Plot electron density in the domain
        
        function qho = plot_density(qho, toplot)
            
            eng_lev = 0;                                 % energy level, only ground state wavefunction can be computed anyway
            n_pp = 100;                                % number of points to use in the plot
            x_vec = linspace(-qho.L/2, qho.L/2, n_pp);
            
            % Analytical Solution
            anal_wf = zeros(1,n_pp);
            for i = 1:n_pp
                anal_wf(i) = qho.solwf(eng_lev, x_vec(i));
            end
            anal_pd = abs(anal_wf) .^ 2;
            
            % Approximate Solution
            PW_freqs = -qho.max_freq : qho.max_freq;
            n_bais_vecs = 2 * qho.max_freq + 1;
            phimat = zeros(n_pp, n_bais_vecs);
            for j = 1:n_bais_vecs
                phimat(:,j) = qho.phi(PW_freqs(j), x_vec);
            end
                
            wf_sln = qho.eig_vecs(:, eng_lev+1);            % extract wavefunction vector from the matrix of eigenvectors
            numer_wf = phimat * wf_sln;
            if numer_wf(round(n_pp/2)) < 0                   % Flip the sign of the wavefunction if it is negative
                numer_wf = - numer_wf;
            end
            numer_pd = abs(numer_wf) .^ 2;
            
            set(0,'defaultlinelinewidth',1.5)
            set(0,'defaultaxeslinewidth',2)

            %clf
            if strcmp(toplot, 'density')
                plot(x_vec, anal_pd, '--')
                hold on
                plot(x_vec, numer_pd, '-')
                hold off
                box('on')
                ylabel('Electron density (1/m)')
            elseif strcmp(toplot, 'wavefunction')
                plot(x_vec, real(anal_wf), x_vec, real(numer_wf))
                ylabel('Wavefunction (1/m^{1/2})')
            else
                disp('nothing')
            end
            xlabel('Position (m)')
            
            legend('Analytical', 'Numerical')
            legend('boxoff')
            ax = gca;
            ax.FontSize = 20;
            
        end
        
        %% Plot energy levels
        function qho = plot_eng_lvls(qho)

            set(0,'defaultlinelinewidth',1.5)
            set(0,'defaultaxeslinewidth',2)

            n_bais_vecs = 2 * qho.max_freq + 1;
            lvl_vec = 0 : n_bais_vecs-1;
            eng_lvl_anal = qho.qho_eng(lvl_vec);
            eng_lvl_numer = diag(qho.eig_vals);
            
            %clf
            
            plot(lvl_vec, eng_lvl_anal, 'o')
            hold on
            plot(lvl_vec, eng_lvl_numer, 'x')
            hold off
            box('on')
            xlabel('Energy level')
            ylabel('Energy (J)')
            legend('Analytical', 'Numerical')
            legend('boxoff')
            legend('Location', 'northwest')
            ax = gca;
            ax.FontSize = 20;
            
        end
        
        %% Harmonic potential
        % Return value given position
        function V = harm_pot(qho, x)
            V = 0.5 * qho.m_e * qho.omega^2 * x .^2;
        end
        
        %% Plane wave
        % Return value given frequency and position
        function wf = phi(qho, n, x)
            xp = x/qho.L + 1/2;
            wf = 1/sqrt(qho.L) * exp( 2 * pi * 1i * n * xp);
        end
        
        %% Analytical Wavefunction Solution
        % x - position
        % n - energy level, starting at 0
        function wf = solwf(qho, n, x)
            alpha = qho.m_e * qho.omega / qho.hbar;
            y = sqrt(alpha) * x;
            wf = (alpha / pi)^0.25 / sqrt(2^n * factorial(n)) * hermite(n,y) * exp(-y^2/2);
        end
        
        %% Analitical Energy
        % n - energy level, starting at 0
        function E_n = qho_eng(qho, n)
            E_n = qho.hbar * qho.omega * (n + 0.5);
        end
        
    end
    
end