function output = get_learning_speed(t_end)
   %keyboard
   % Params
   var_sig   = 100;              % Trialwise variance in the signal
   mu        = 1.4*sqrt(var_sig);% Amount of signal in each trial
   n         = 100;              % Total number of neurons
   phi_list  = 0:0.01:1;         % Set of phi values to investigate
   alpha     = 0.0001;           % Learning rate   
   delta     = 0.5;              % Reward prediction error strength
   if nargin < 1
   t_end     = 100;              % Number of trials in simulation
   end

   % Which phi's to plot
   phi_sel   = [0, 0.02, 0.2];

   % Pools
   m1 = [ones( n./2,1) ; zeros(n./2,1)];
   m2 = [zeros(n./2,1);  ones(n./2,1)];
   m1h = m1./norm(m1);
   m2h = m2./norm(m2);
   
   % Signal dimension
   v_s = (m1 - m2)./norm(m1 - m2);
   
   % Conserved high variance dimension (v "plus")
   v_p = (m1 + m2)./norm(m1 + m2);

   % Get a noise representative
   xi_perp = get_noise_rep(n, m1h, m2h);

   % Compute everything, for each choice of phi
   k = 1;
   for phi = phi_list
      % Generate the covariance matrix
      [S, a] = get_cov_mat(var_sig, phi, n);

      % Empirical pedeltandicular variance (albeit with a lazy v_p...)
      tot_var_e(k) = xi_perp'*S*xi_perp*(n-2) + v_p'*S*v_p;

      % Algebraic pedeltandicular variance
      tot_var_a(k) = n*a - 1*var_sig;
      
      % Algebraic one dimensional signal detection parameters
      % Signal separation distance
      d(:,k)   = get_d_dynamics(t_end, alpha, delta, mu);
      
      % In-signal-dimension noise standard deviations
      [sig(:,k), sig_star(:,k)] = get_sig_dynamics(t_end, var_sig, mu, tot_var_a(k), n, alpha, delta);
      
      % Accuracy as a function of these
      acc(:,k) = get_accuracy(d(:,k), sig(:,k));
      
      % D-prime statistic, b.c. who doesn't want to know?
      dp(:,k)  = d(:,k)./sig(:,k);

      % Increment the counter
      k = k+1;
   end
   
   % Plot the results
   output=plot_all(phi_list, tot_var_a, tot_var_e, dp, d, sig, acc, phi_sel, alpha, delta, mu, sig_star)

end

function xi = get_noise_rep(n, m1h, m2h)
   %  
   xi = mvnrnd(zeros(n,1), eye(n))';
   
   % Orthogonalize it to the signal...
   %
   % And since the projection onto ones(200,1) really matters
   % orthogonalize it to this as well, since I don't want to 
   % bother averaging it's contribution back to unity...
   xi = xi - m1h*xi'*m1h - m2h*xi'*m2h;
   
   % Normalize
   xi = xi./norm(xi);
end

function [S, a] = get_cov_mat(rowsum, phi, n)


   % Construct covariance matrix
   a   = 2*rowsum/(2+(n-2)*phi);
   m1 = [ones( n./2,1) ; zeros(n./2,1)];
   m2 = [zeros(n./2,1);  ones(n./2,1)];

   S = phi*a*m1*m1' + phi*a*m2*m2' + (1-phi)*a*eye(n);

end

function d = get_d_dynamics(t_end, alpha, delta, mu)

for t = 1:t_end
   ws = (t-1)*alpha*delta*mu;
   d(t) = 2*ws*mu;
end

end

function [sig, sig_star]= get_sig_dynamics(t_end, var_sig, mu, tot_var_perp, n, alpha, delta)

for t = 1:t_end
   ws = alpha*delta*(t-1)*mu;
   wp = alpha*delta*sqrt((t-1)*tot_var_perp);
   
   wa = alpha*delta*sqrt((t-1)*var_sig);
      
   sig_star(t) = sqrt(ws^2*var_sig);
   sig(t)      = sqrt(ws^2*var_sig + wp^2*tot_var_perp/(n-1));
end

end

function acc = get_accuracy(d, sigma)

   acc = 1 - normcdf(-d/2, 0, sigma);
   
   % The d = 0, sigma = 0 case should yield
   % 50% accuracy, since sigma is not truly zero
   % whereas d is.
   if d(1) == 0
      acc(1) = 0.5;
   end

end

function [output] = plot_all(phi_list, tot_var_a, tot_var_e, dp, d, sig, acc, phi_sel, alpha, delta, mu, sig_star)
   
   % Two misc. plot requirements 
   phi_inds = find(sum(phi_list' == phi_sel,2))';
   namelist = cell(1,length(phi_inds));
   
   % Plot stuff
   %figure()
   
   subplot(1,3,1)
   plot(phi_list, sqrt(tot_var_a)); hold on
   plot(phi_list(1:2:101), sqrt(tot_var_e(1:2:101)), 'o')
   title('Total Perpendicular Std.')
   grid on
   xlabel('\phi')
   ylabel('\sigma(\xi_\perp)')
   legend({'Algebraic','Empirical'})
   
%    subplot(2,3,2)
%    hold on
    k = 1;
    for i = phi_inds
      %plot(dp(:,i), '-');
      namelist{k} = ['\phi = ' sprintf('%0.2f', phi_list(i))];
      k = k+1;
    end
%    grid on;
%    title('d-prime Dynamics' )
%    xlabel('Time [steps]')
%    ylabel('d-prime')
%    legend(namelist,'Location', 'SouthEast')
   
%    subplot(2,2,3)
%    hold on;
%    for i = phi_inds
%       plot(sig(:,i))
%    end
%    grid on
%    xlabel('Time')
%    ylabel('Standard Deviation')
%    legend(namelist)
   
   subplot(1,3,3)
   for i = phi_inds
      plot(acc(:,i));   hold on
   end
   grid on
   title('Analytic Accuracy')
   xlabel('Time')
   ylabel('Percent Correct')
   legend(namelist, 'Location', 'SouthEast')
%    
%    subplot(2,1,2)
%    corder = get(gca,'ColorOrder');
%    
%    xs = -60:60;
%    
%    t   = 100;
%    f   = t*alpha*delta*mu;
%    for i = 1:length(phi_inds)
%    
%       plot(xs, normpdf(xs, -d(t,i)/f/2, sig(t,i)/f), 'Color', corder(i,:)); hold on;
%       plot(xs, normpdf(xs,  d(t,i)/f/2, sig(t,i)/f), 'Color', corder(i,:))
%    end
%    title('Signal Distributions Along Signal Dimension')
%    xlabel('\Delta y/(t\alpha\delta\mu)')
%    ylabel('PDF')
%    grid on

output.acc=acc;
output.phi_inds=phi_inds;
output.phi_list=phi_list;
output.tot_var_e=tot_var_e;
output.tot_var_a=tot_var_a;
output.namelist=namelist;



end