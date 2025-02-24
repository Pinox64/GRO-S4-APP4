function [y, t] = impulse_rational(sys, t)
% impulse_rational computes the impulse response using partial fraction 
% decomposition and a Pade approximation to handle any input delay.
%
% Usage:
%   [y, t] = impulse_rational(sys, t)
%
% Inputs:
%   sys : a transfer function (tf object) that may include an InputDelay.
%   t   : a time vector for evaluation.
%
% Outputs:
%   y   : the computed impulse response evaluated at the time vector t.
%
% Note: This function assumes that sys is a SISO system.
%
% Steps:
% 1. If sys has an InputDelay, approximate it using a Pade expansion.
% 2. Remove the InputDelay (set it to 0) and multiply sys by the Pade
%    approximation (delay_tf) so that the delay is incorporated in the polynomial.
% 3. Extract the numerator and denominator coefficients.
% 4. If the system is improper (numerator degree >= denominator degree),
%    perform polynomial division.
% 5. Use residue() to perform partial fraction decomposition.
% 6. Reconstruct the impulse response, properly handling repeated poles and 
%    any polynomial (feedthrough) terms.
% 7. If there is a direct feedthrough (system degree equal), add the appropriate
%    delta contribution at t=0.

    % Step 1: Handle InputDelay using a Pade approximation if needed.
    delay = sys.InputDelay;
    if delay > 0
        pade_order = 5;  % You may adjust the order as needed.
        [Nd, Dd] = pade(delay, pade_order);
        delay_tf = tf(Nd, Dd);
        sys.InputDelay = 0;  % Remove the InputDelay from the original system.
        sys = sys * delay_tf;  % Multiply to incorporate the delay.
    end

    % Step 2: Extract numerator and denominator coefficients (as vectors).
    [num, den] = tfdata(sys, 'v');

    % Step 3: Check for direct feedthrough (if degrees are equal).
    direct_feedthrough = (length(num) == length(den));

    % Step 4: Ensure the system is proper.
    if length(num) >= length(den)
        [Q, R] = deconv(num, den);  % Perform polynomial division.
        disp('System was improper; performing polynomial division.');
        num = R;  % Use the remainder for the residue calculation.
    end

    % Step 5: Compute partial fraction decomposition.
    [Rres, Pres, K] = residue(num, den);

    % Step 6: Reconstruct the impulse response from the partial fractions.
    y = zeros(size(t));
    i = 1;
    while i <= length(Pres)
        mult = sum(Pres == Pres(i));  % Count repeated poles.
        if mult == 1
            y = y + Rres(i) * exp(Pres(i) * t);
        else
            % For repeated poles, include terms with t^k.
            for k = 0:(mult-1)
                y = y + (Rres(i) / factorial(k)) * (t.^k) .* exp(Pres(i) * t);
            end
        end
        i = i + mult;
    end

    % Step 7: Add any polynomial (direct feedthrough) terms from K.
    if ~isempty(K)
        y = y + polyval(K, t);
    end

    % Step 8: Handle direct feedthrough by adding the impulse contribution at t = 0.
    if direct_feedthrough
        % Approximate the delta impulse contribution at t = 0.
        y(t == 0) = y(t == 0) + num(1) / den(1);
    end
end