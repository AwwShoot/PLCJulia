# Written by Alice and Zealand

using Printf

function sgn(N)
	if (N < 0)
		return -1
	elseif (N == 0)
		return 0
	else
		return 1
	end
end

input = readlines()

#println(input)

satt_i = split(input[1], " ")
satt_j = split(input[2], " ")
satt_k = split(input[3], " ")
satt_l = split(input[4], " ")

# println(satt_i)
# println(satt_j)
# println(satt_k)
# println(satt_l)

i_coords = Float64[]
j_coords = Float64[]
k_coords = Float64[]
l_coords = Float64[]

for i in range(1,3)
	println(i)
	curr_coord = parse(Float64, satt_i[i])
	push!(i_coords, curr_coord)

    curr_coord = parse(Float64, satt_j[i])
	push!(j_coords, curr_coord)

    curr_coord = parse(Float64, satt_k[i])
	push!(k_coords, curr_coord)

    curr_coord = parse(Float64, satt_l[i])
	push!(l_coords, curr_coord)
end

println(i_coords)
println(j_coords)
println(k_coords)
println(l_coords)

#x four distances
x_ji = j_coords[1] - i_coords[1]
x_ki = k_coords[1] - i_coords[1]
x_jk = j_coords[1] - k_coords[1]
x_lk = l_coords[1] - k_coords[1]

# y four distances 
y_ki = k_coords[2] - i_coords[2]
y_ji = j_coords[2] - i_coords[2]
y_lk = l_coords[2] - k_coords[2]
y_jk = j_coords[2] - k_coords[2]

# z four distances
z_ji = j_coords[3] - i_coords[3]
z_ki = k_coords[3] - i_coords[3]
z_jk = j_coords[3] - k_coords[3]
z_lk = l_coords[3] - k_coords[3]

# Speed of light (m/s) use meters per nanosecond since times are provided in nanoseconds
sLight = 0.299792458

# Iterate through the remaining lines of input (5 and onward) 
for lines in range(5, length(input))

    times = split(input[lines], " ")
# Timestamps    
    T_i = parse(Float64, times[1])
    T_j = parse(Float64, times[2])
    T_k = parse(Float64, times[3])
    T_l = parse(Float64, times[4])

    #=
	The distance between sattelites can be calculated by multiplying the difference between
	signal arrival times by the speed of light. Assuming all sattelites sent their signals at the same time,
	the time difference can only be explained by one being closer or further.
	
	This is the R value in the equation set.
	=#
    R_i = (sLight * T_i)
	R_j = (sLight * T_j)
	R_k = (sLight * T_k)
	R_l = (sLight * T_l)

    R_ij = R_i - R_j
    R_ik = R_i - R_k
    R_kj = R_k - R_j
    R_kl = R_k - R_l

    # set of quantities - Equations 1 -> 6 on the paper
    X_ijy = (R_ij * y_ki) - (R_ik * y_ji)
    X_ikx = (R_ik * x_ji) - (R_ij * x_ki)
    X_ikz = (R_ik * z_ji) - (R_ij * z_ki)
    X_kjy = (R_kj * y_lk) - (R_kl * y_jk)
    X_klx = (R_kl * x_jk) - (R_kj * x_lk)
    X_klz = (R_kl * z_jk) - (R_kj * z_lk)

    # the squared sums of the coordinates of the sattelites are used elsewhere. God knows why.
	S_i = (i_coords[1]^2) + (i_coords[2]^2) + (i_coords[3]^2)
	S_j = (j_coords[1]^2) + (j_coords[2]^2) + (j_coords[3]^2)
	S_k = (k_coords[1]^2) + (k_coords[2]^2) + (k_coords[3]^2)
	S_l = (l_coords[1]^2) + (l_coords[2]^2) + (l_coords[3]^2)
	
    # final definition of some incomprehensible mass of letters.
	R_ij2 = (R_ij^2) + S_i - S_j
	R_ik2 = (R_ik^2) + S_i - S_k
	R_kj2 = (R_kj^2) + S_k - S_j
	R_kl2 = (R_kl^2) + S_k - S_l

    # getting A from the given equations since they couldn't be bothered to move the vars themselves
	A = X_ikx / X_ijy
	B = X_ikz / X_ijy
	C = X_klx / X_kjy
	D = X_klz / X_kjy

    # I don't even know what's a variable anymore and what's merely part of another var
	E = ((R_ik * R_ij2) - (R_ij * R_ik2)) / (2 * X_ijy)
	F = ((R_kl * R_kj2) - (R_kj * R_kl2)) / (2 * X_kjy)



	# more following equations
	G = (D - B) / (A - C)
	H = (F - E) / (A - C)

    println(G)
	println(H)
	# I and J I guess
	I = (A * G) + B
	J = (A * H) + E
	
	
	println(I)
	println(J)

    # Once more I am truly unsure of what is and isn't a variabel. Why are these only using ki?
	K = R_ik2 + (2 * x_ki * H) + (2 * y_ki * J)
	L = 2 * ((x_ki * G) + (y_ki * I) + z_ki)
	

	R2_ik4 = 4 * (R_ik^2)

	G2 = G^2
	I2 = I^2
	L2 = L^2
	# Almost there...
	M = R2_ik4 * (G2 + I2 + 1) - L2
	N = ((8 * (R_ik^2)) * (G * (i_coords[1] - H) + (I * (i_coords[2] - J)) + i_coords[3])) + (2 * L * K)
	O = (R2_ik4 *(((i_coords[1] - H)^2) + ((i_coords[2] - J)^2) + (i_coords[3]^2))) - (K^2)
	
	println(M)
	println(O)
	
	
	# Finally, a formula I recognize
	Q = N + ((sgn(N))*(sqrt((N^2) - (4 * M * O))))
	


	# Finally, The coordinates.
	z_pos = Q / (2 * M)
	z_neg = O * 2 / Q
	x_pos = G * z_pos + H
	x_neg = G * z_neg + H
	y_pos = I * z_pos + J
	y_neg = I * z_neg + J
	
	# and the distances from the origin for some reason?
	pos_distance = (x_pos ^ 2 + y_pos ^ 2 + z_pos ^ 2) ^ 0.5
	neg_distance = (x_neg ^ 2 + y_neg ^ 2 + z_neg ^ 2) ^ 0.5
	
	@printf "g= %9.2e, h= %9.2e, j= %9.2e, m= %9.2e, o= %9.2e \n" G H J M O
	@printf "+) x= %10.2e, y= %10.2e, z= %10.2e; I= %10.2e\n" x_pos y_pos z_pos pos_distance
	@printf "+) x= %10.2e, y= %10.2e, z= %10.2e; I= %10.2e\n" x_neg y_neg z_neg neg_distance
    
end