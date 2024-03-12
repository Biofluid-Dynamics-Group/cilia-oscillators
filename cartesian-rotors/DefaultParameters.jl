
# Dimensionless parameters
Λ = 0.1  # = λ*d/f_drive
num_beats = 120  # = final_time/T
R = 0.5  # = r_0/d
L = 2  # = l/d
A = 0.01  # = a/d

μ = 1e3  # kg/μm/s
d = 7  # μm
T = 1/33  # s

r_0 = R*d
a = A*d
gamma_0 = 6*π*μ*a
final_time = num_beats*T
f_drive = 2*π*r_0*gamma_0/T
λ = Λ*f_drive/d
η = λ
x_0 = [0, 0, d]
