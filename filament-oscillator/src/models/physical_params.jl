# All lengths are given divided by a
# All times are given divided by T


# Single cilium geomerty
a = 1.0  # Sphere radius, reference
N = 5  # Number of discretisation spheres per cilium
sphere_space = 0.2  # Fraction of radii between spheres
L = N*(2 + sphere_space)*a  # Cilium length w/r to sphere radius

# Cilia system geomerty
M = 1  # Number of cilia
d = 50  # Amount of radii between cilia
φ = 0 # 12.5*π/180.0  # Angle of cilia beat plane

# Single cilium reference beat dynamics
T = 1.0  # Beat period, reference
θ_0 = π/2.1  # Maximum bending angle
f_eff = 0.3  # Effective stroke fraction
f_ψ = 0.85  # Fraction of recovery controlled by travelling wave
f_w = 0.4  # Travelling wavelength w/r to cilium length
orientation = π/2.0  # Normal to attachment plane

# Fluid dynamics
μ = 1.0  # Fluid viscosity, reference
κ_tilde = 1.2  # Ratio of bending stiffness to viscous drag
κ_b = κ_tilde*μ*T  # bending stiffness

