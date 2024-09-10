# Filename: Github_Memory_Num_MoveBirth.jl
# Author: Yifei Li
# Harbin institute of Technology, Harbin, China, September 2024 
# The script contains the numerical method for solving memory-based reaction-diffusion equation

using DifferentialEquations
using Plots
using DelimitedFiles
function bc_model(du, u, h, para, t)
    N, D, p, tau, r = para
    hist = h(para, t - tau)
    for i=2:N-1
        du[i] = D*(u[i-1]-2*u[i]+u[i+1]-
                hf(u[i],u[i+1],p)*(hist[i+1]-hist[i])+
                hf(u[i],u[i-1],p)*(hist[i]-hist[i-1]))+r*u[i]*(1-u[i])
    end

    i=1
    du[i] = D*(u[N]-2*u[i]+u[i+1]-
            hf(u[i],u[i+1],p)*(hist[i+1]-hist[i])+
            hf(u[i],u[N],p)*(hist[i]-hist[N]))+r*u[i]*(1-u[i])
    i=N
    du[i] = D*(u[i-1]-2*u[i]+u[1]-
            hf(u[i],u[1],p)*(hist[1]-hist[i])+
            hf(u[i],u[i-1],p)*(hist[i]-hist[i-1]))+r*u[i]*(1-u[i])
end
function hf(a,b,p)
    return (2*p*a*(1-a)+2*p*b*(1-b))/2
end

h(para,t)=[zeros(1,80) ones(1,40) zeros(1,80)]

tau =20
lags = [tau]
p = 0.5
N=200
D=1/4
r=0.01
ts=200

para = (N, D, p, tau, r)
tspan = (0.0, ts)
u0=[zeros(1,80) ones(1,40) zeros(1,80)]

prob = DDEProblem(bc_model, u0, h, tspan, para; constant_lags = lags)
alg = MethodOfSteps(Tsit5())
sol = solve(prob, alg)
a=range(0, 200, length=200)
plot(a,vec(sol.u[end]))