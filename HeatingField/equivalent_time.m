%% Authored by Joshua Soneson 2007
function[T]=equivalent_time(T,N,T0)
% Computes the integrand of the CEM43 thermal dose integral for 
% an N-vector of temperatures T.

diff = 43-T0;
for n=1:N
  if(T(n)<=diff)
    T(n) = 0.25^(diff-T(n));
  else
    T(n) = 0.5^(diff-T(n));
end; end
