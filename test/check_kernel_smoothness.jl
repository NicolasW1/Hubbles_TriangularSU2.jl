using Hubbles_TriangularSU2
using Plots

# very preliminary

thr = Hubbles_TriangularSU2.float_threshold

x1 = range(-2*thr, 2*thr, length=500)
y1_1 = @. Hubbles_TriangularSU2.kernel_Kplus(0.05, 0.1, 0.1+x)
y1_2 = @. Hubbles_TriangularSU2.kernel_Kplus(0.05, 0.001, 0.001+x)
y1_3 = @. Hubbles_TriangularSU2.kernel_Kplus(0.01, 0.5, 0.5+x)
y1_4 = @. Hubbles_TriangularSU2.kernel_Kplus(0.01, -0.5, -0.5-x)

plot([(x1, y1_1), (x1, y1_2), (x1, y1_3), (x1, y1_4)], layout = 4, legend=false, ticks=false)

x2 = range(-1.0, 1.0, length=1000)
y2 = @. Hubbles_TriangularSU2.kernel_Kminus(10^-2, x2, 0.1)

scatter(x2, y2, legend=false)

plot(x1, y1_1, xscale=:log10, legend=false)
plot(x1, y1_2, xscale=:log10, legend=false)

