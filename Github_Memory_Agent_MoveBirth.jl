# Filename: Github_Memory_Agent_MoveBirth.m
# Author: Yifei Li
# Harbin institute of Technology, Harbin, China, May 2024
# The script contains the simulation of the agent-based model generating data
# of the evolution process, which are consistent to the solutions of memory-based PDE model.

using DelimitedFiles
using Plots, SpecialFunctions, Random, LaTeXStrings, DifferentialEquations

function simulation(para)
    p, r, N, t, tau = para
    #initial distribution
    #1/0 in the i-st row means that the lattice is/isn't occupied at i-th step
    global Seed=zeros(Int64,t+2,N)
    Seed[2,81:120].=1
    global total=sum(Seed[2,:])
    global Seed_coor=Array{Int64}(undef, total)
    
    global count=0
    for i=1:N
        global count
        if Seed[2,i]>0
            count=count+1
            Seed_coor[count]=i;
        end
    end
    for i=2:t+1
        Seed[i+1,:]=Seed[i,:]
        local Temp=Seed[i,:]
        total=length(Seed_coor)
        for j=1:total
            index=rand(1:total)
            selected=Seed_coor[index]
            if selected==1
                left=N
            else
                left=selected-1
            end
            if selected==N
                right=1
            else
                right=selected+1
            end
            if i<tau+2
                before=1
            else
                before=i-tau
            end
            possibility_left=(1-(Seed[before,selected]-Seed[before,left])*p)/4
            possibility_right=(1+(Seed[before,right]-Seed[before,selected])*p)/4
            tt=rand(1)[1]
            if tt<possibility_left
                if Temp[left]<1
                    Seed[i+1,selected]=0
                    Temp[selected]=0
                    Seed[i+1,left]=1
                    Temp[left]=1
                    Seed_coor[index]=left
                end
            elseif tt<(possibility_left+possibility_right)
                if Temp[right]<1
                    Seed[i+1,selected]=0
                    Temp[selected]=0
                    Seed[i+1,right]=1
                    Temp[right]=1
                    Seed_coor[index]=right
                end
            end
        end
        for j=1:total
            index=rand(1:total)
            selected=Seed_coor[index]
            if selected==1
                left=N
            else
                left=selected-1
            end
            if selected==N
                right=1
            else
                right=selected+1
            end
            if i<tau+2
                before=1
            else
                before=i-tau
            end
            possibility_left=r/2
            possibility_right=r/2
            tt=rand(1)[1]
            if tt<possibility_left
                if Temp[left]<1
                    Seed[i+1,left]=1
                    Temp[left]=1
                    append!(Seed_coor,left)
                end
            elseif tt<(possibility_left+possibility_right)
                if Temp[right]<1
                    Seed[i+1,right]=1
                    Temp[right]=1
                    append!(Seed_coor,right)
                end
            end
        end
    end
    return Seed
end
#plot!(tspan,recordAd[:],label="Ad")
#parameters
#memory effect coefficient
p=0.5
#birth rate
r=0.01
#total amount of lattice
N=200
#total time steps
t=200
#time delay
tau=40;
para = (p, r, N, t, tau)
global Result_Seed=simulation(para)

for c=1:999
    global Result_Seed
    New=simulation(para)
    Result_Seed=Result_Seed+New
end
a=range(0, 200, length=200)
Result_Seed=Result_Seed/1000
plot(a,Result_Seed[end,:],label="Agents")