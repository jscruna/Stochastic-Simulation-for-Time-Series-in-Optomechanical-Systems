using Distributed
@time begin
addprocs(5) #the number of cores used for the parallel computing in Julia programming language
fs=1e6; #initial frequency sample
reduc=100; # by using this factor we can manipulate the final resolution of the data
fs=fs/reduc; # frequency sample
dt=1/fs; #time interval
T=10.1; #total amount of time for the simulation
steps=round(Integer,T/dt); #numer of steps defined by taking an interval
#We also estimate the time of computation from here to evaluate the waiting time
@everywhere begin
    using Random
    using Statistics
    using LinearAlgebra
    function simulate_trajectory(steps::Int, dt::Float64, seed::Int)
        #We generate three seeds for our random variables with a separation for independence
        rng = MersenneTwister(1233+seed);
        rand1 = randn(rng, steps);
        rng = MersenneTwister(2243+seed);
        rand2 = randn(rng, steps);
        rng = MersenneTwister(3253+seed);
        rand3 = randn(rng, steps);
        #Optomechanical parameters
        N_th=19; #Phonon number occupation laser
        hbar=6.626e-34/2/pi; #Planck constant
        k_B=1.380648e-23; #Boltzman constant
        c=299.792458e6; #Light speed
        lambda=1064e-9; #Wavelength of the laser
        omega_laser=c/lambda*2*pi;
        Temp=10.8e-3; #Temperature of the mechanical mode
        L=0.098; #Length of the optical cavity
        M=7.71e-6; #mass of the mechanical system
        FSR=c/L; finesse=1849;
        kappa=FSR/finesse*pi; #Optical loss
        kappa_in_over_all=0.2247;
        detune=0.0299*2;
        omega_m=2*pi*280; #mechanical resonance frequency
        gamma_m=2*pi*1.1; #mechanical dissipation
        n_th=8e5 #mean phonon occupation
        eta=0.92;#detector efficiency
        g=-2*pi*3.2*1e4 #optomechanical coupling
        Detune=0.02092; #normalized detunning
        kappa_p=2*pi*1.64e6; #optical loss
        k=kappa_p;
        Detune2=Detune*kappa_p;
        delta=Detune2;
        
        #Variables to consider the change in the system due to measurement
        #Brownian force
        n_x=2*gamma_m*(2*n_th+1)+16*g^2*(2*N_th+1)/((1+4*Detune^2)*kappa_p);
        #Measurement rate
        lambda_x= 64*g^2*eta*Detune^2/((2*eta*N_th+1)*kappa_p*(1+4*Detune^2)^2);
        #This value is related with the degree of squeezing level
        sigma_x= (-32*g^2*eta*Detune)*(2*N_th+1)/((1+4*Detune^2)^2*kappa_p*(2*eta*N_th+1));
        #Modified resonance mechanical frequency due to measurement
        omega_x=sqrt(sqrt(omega_m^4+2*sigma_x*omega_m^3+n_x*lambda_x*omega_m^2))
        #Modified mechanical dissipation due to measurment
        gamma_x=sqrt(gamma_m^2-2*omega_m*(omega_m+sigma_x)+2*omega_x^2)
        
        #We finished to set the parameters;
        T=10.1; #Measurement time
        fs=1/dt; #frequency sample
        size=T*fs; #Size for our vectors in the space state representation
        dt=1/fs;
        time=0:dt:T; #vector time
        steps=round(Integer,T/dt);
        delta1=delta*ones(steps,1);
        #Matrix for the sensing noise    
        Raux=(sqrt(eta)*(k^2-4*delta^2)/(k^2+4*delta^2))^2*(2*N_th+1) + (4*delta*k*sqrt(eta)/(k^2+4*delta^2))^2*(2*N_th+1) +(1-eta);
        #Matrix for the cross-correlation noise
        Saux=[0; 4*g*k*sqrt(k)*sqrt(eta)*(k^2-4*delta^2)/(k^2+4*delta^2)^2*(2*N_th+1) + 8*g*delta*sqrt(k)*4*delta*k*sqrt(eta)/(k^2+4*delta^2)^2*(2*N_th+1)];
        #Dynamic matrix for the mechanical system
        A=[0 omega_m;-omega_m -gamma_m];
        #Vector for position and momentum
        xorig2=zeros(2,steps+1);
        #Vector for photodetection
        obs=zeros(1,steps);
        #Initial vector value for position and momentum
        xorig2[:,1]=[0;0];
        #Variance in the mechanical system
        Qaux=[0 0;0 n_x];
        #Effective matrix according to theory
        Qf=Qaux-Saux*(Raux)^(-1)*Saux';
        noise_r=zeros(1,steps);# the definitions of this vector is important for size compatibility
        noise_r[1,1]=sqrt(Raux/dt)*rand3[1];
        C_matrix_aux=[-A' [1 0;0 1] zeros(2,2) zeros(2,1);zeros(2,2) -A' Qf  zeros(2,1); zeros(2,2) zeros(2,2) A Saux*Raux^(-1)*noise_r[1,1]; 0 0 0 0 0 0 0];
        phik=exp(A*dt);#exp(C_matrix_aux*dt)[1:2,1:2];#dynamics for a time t+dt
        qk=(exp(A*dt)'*exp(C_matrix_aux*dt)[3:4,5:6]);#variance for a time t+dt
        ck=cholesky((qk+qk')/2).L;#Cholesky applied according to the theory
        for rr in 1:steps
        noise_r[1,rr]=sqrt(Raux/dt)*rand3[rr];
        C_matrix_aux=[-A' [1 0;0 1] zeros(2,2) zeros(2,1);zeros(2,2) -A' Qf  zeros(2,1); zeros(2,2) zeros(2,2) A Saux*Raux^(-1)*noise_r[1,rr]; 0 0 0 0 0 0 0];
        uk=exp(C_matrix_aux*dt)[5:6,7];#input
        #state-space representation for our optomechanical system
        xorig2[:,rr+1]=phik*xorig2[:,rr]+ck*[rand1[rr];rand2[rr]]+uk; #mechanical system
        obs[1,rr]=[-8*delta*g*sqrt(k*eta)*(k^2+4*delta^2)^(-1), 0]'*xorig2[:,rr+1]+ noise_r[1,rr]; #photodetector equation
        end
        return xorig2[1,:]
    end 
end
#the number of trajectories
num_trajectories = 1000;
seeds = [1 + i for i in 1:num_trajectories];
#By using parallel computing, we compute trajectories for position (In this case the number of cores used in 5)
trajectories = pmap((seed) -> simulate_trajectory(steps, dt, seed), seeds);
#We obtain the variances at time t by considering the information of each trajectory
variances = mapreduce(i -> var(map(traj -> traj[i], trajectories)), vcat, 1:steps)
#We save the information to manipulate the data using MATLAB
using　MAT
matwrite("sim.mat",Dict("variances"=>variances))
end
#We save the information to manipulate the data using Python
using CSV
using DataFrames
data = DataFrame(
    Variances = variances
)
CSV.write("variances.csv", data)
data2 = DataFrame(
    Trajectories1 = trajectories[1],
    Trajectories2 = trajectories[2],
    Trajectories3 = trajectories[3]
)
CSV.write("time_series_sample.csv",data2)
nothing