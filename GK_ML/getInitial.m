function [ initial ] = getInitial()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
count = 1;
initial = zeros(26,1);
initial(count) = 0.537489130000000; count = count+1;    % C
initial(count) = 1; count = count+1;                    % D
initial(count) = 2.819586840000000; count = count+1;    % F
initial(count) = 0.169757100000000; count = count+1;    % G
initial(count) = 0.141539270000000; count = count+1;    % I
initial(count) = 0; count = count+1;                    % In
initial(count) = 5.661570770000001; count = count+1;    % K
initial(count) = 0.333333330000000; count = count+1;    % L
initial(count) = 1.415392690000000; count = count+1;    % N
initial(count) = 1.402779960000000; count = count+1;    % Ne
initial(count) = 0.012612740000000; count = count+1;    % Nn
initial(count) = 0.760019200000000; count = count+1;    % Pm
initial(count) = 1; count = count+1;                    % Q
initial(count) = 1.012601010000000; count = count+1;    % Rk
initial(count) = 1; count = count+1;                    % U
initial(count) = 0.848785500000000; count = count+1;    % Y
initial(count) = 0.848785500000000; count = count+1;    % Ym
initial(count) = 3.709888960000000; count = count+1;    % Z
initial(count) = 0.025000000000000; count = count+1;    % delta
initial(count) = 1.511020840000000; count = count+1;    % eta
initial(count) = 1.010101010000000; count = count+1;    % i
initial(count) = 1; count = count+1;                    % infl
initial(count) = 1; count = count+1;                    % inflstar
initial(count) = 0.003739780000000; count = count+1;    % nu
initial(count) = 4; count = count+1;                    % phi
initial(count) = 1.020101010000000;                     % x
end

