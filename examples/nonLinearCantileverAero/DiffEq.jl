# Include libreries (DifferentialEquations, BoundaryValueDiffEq, Plots, DataFrames, must be included into manifest)
using BoundaryValueDiffEq, Plots, FileIO, DataFrames, CSV
# Define problem parameters
L = 10;
d = L/100; 
J = pi * (d^4) / 32;
Iyy = J/2;
Izz = Iyy;
E = 1e9;

# Define wind paramets
c_d = 1.2;
vw = 5;
rhoAire = 1.225
dimCaracteristic = d
q0 = 1/2 * rhoAire * vw^2 * dimCaracteristic 
qz = q0 * c_d; 

# Define numerical method params
tspan = (0.0,L)
deltaX = 0.1
aboslutoTolerance = 1e-7
relativeTolerance = 1e-9

# Define differential equations of the problem
function nonLinearStaticCantilever!(du,u,p,t)
    Vz, My, θy, uz = u
    du[1] = -qz * cos(θy)^3
    du[2] = Vz
    du[3] = My / (E*Iyy)
    du[4] = θy
end

function bc2!(residual, u, p, t) # u[1] is the beginning of the time span, and u[end] is the ending
    residual[1] = u[end][1] + 0 # the sheer force solution at the end of the x span should be 0
    residual[2] = u[end][2] + 0 # the benging moment solution at the end of the x span should be 0
    residual[3] = u[1][3] + 0 # the benging moment solution at the end of the x span should be 0
    residual[4] = u[1][4] + 0 # the displacement at x = 0 span should be 0
end

# Register elapsed time
tiempo = @elapsed begin
# The MIRK4 solver is necessary for TwoPointBVProblem
bvp2 = TwoPointBVProblem(nonLinearStaticCantilever!, bc2!, [0,0,0,0], tspan)
sol = solve(bvp2, MIRK4(), dt = deltaX, abstol = aboslutoTolerance, retol = relativeTolerance, save_everystep = true, alg_hints=[:stiff]) 
end 
# Plot solutions
cd("output/")
lw = 5; plotDensity = 1000;
plotVz = plot(sol, plotdensity=plotDensity, vars=(1), fmt = :png, title ="Sheer force in axis z", xaxis="x (m)", yaxis="Vz (N)", label ="Vz",lc=[:red], linewidth = lw)
png(plotVz,"plotVz")
plotMy = plot(sol, plotdensity=plotDensity, vars=(2), fmt = :png, title ="Bending moment in axis y", xaxis="x (m)", yaxis="My (N.m)", label ="My", lc=[:black], linewidth=lw)
png(plotMy,"plotMy")
plotThetaY = plot(sol, plotdensity=plotDensity, vars=(3), fmt = :png, title ="Angle in axis y", xaxis="x (m)", yaxis="θy (rad)", label ="θy", lc=[:blue], linewidth=lw)
png(plotThetaY,"plotThetaY")
plotUz = plot(sol, plotdensity=plotDensity, vars=(4), fmt = :png, title ="Displacement in axis z", xaxis="x (m)", yaxis="uz (m)", label ="uz", lc=[:green], linewidth=lw)
png(plotUz,"plotUz")

# Export results functions
function createTxt(x; name)
    open(name, "w") do io
        [print(io, xi, " ,") for xi in x]
    end
end

function sliceMatrix(A)
    m, n = size(A)
    B = Array{Array{eltype(A), 1}, 1}(undef, m)
    
    for i = 1:m
        B[i] = A[i, :]
    end
    return B
end
# Export solution in CSV format
df_solution = DataFrame(sol)
CSV.write("solJDiffEq.csv", df_solution)

# Export solutions into .txt file
solArray = sliceMatrix(sol)
# createTxt(solArray[3]; name="solJDiffEq_thetaY.txt")
# createTxt(solArray[4]; name="solJDiffEq_uz.txt")
# createTxt(range(0, L, length = trunc(Int, round(L / deltaX +1))); name="solJDiffEq_xcords.txt")

defEnergy = sum( solArray[2].^2 / (2*E*Iyy)*deltaX )

cd("..")

return solArray[3][end], tiempo


