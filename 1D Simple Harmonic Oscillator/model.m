function [frequency] = model(stiffness,mass)

%MODEL: 1-Degree of freedom simple harmonic oscillator system

frequency = sqrt(stiffness./mass);

end