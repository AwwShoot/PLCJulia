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

satt_i = readline()
satt_j = readline()
satt_k = readline()
satt_l = readline()

# we need to split the coordinates from the readlines
# so we can get each individual x, y, z value
i_unparsed = split(satt_i, " ")
#initialize empty array of parsed float numbers
i_coords = Float64[]

for i in i_unparsed
	print(i)
	curr_coord = parse(Float64, i)
	push!(i_coords, curr_coord)
end

# repeat the above steps for the rest of the satilites
j_unparsed = split(satt_j, " ")
j_coords = Float64[]

for j in j_unparsed
	curr_coord = parse(Float64, j)
	push!(j_coords, curr_coord)
end

k_unparsed = split(satt_k, " ")
k_coords = Float64[]

for k in k_unparsed
	curr_coord = parse(Float64, k)
	push!(k_coords, curr_coord)
end

l_unparsed = split(satt_l, " ")
l_coords = Float64[]

for l in l_unparsed
	curr_coord = parse(Float64, l)
	push!(l_coords, curr_coord)
end
	
# x four distances
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

# Speed of light (m/s) # use meters per nanosecond since times are provided in them
SPEED_OF_LIGHT = 299792458

# for now, im initiating an array of zeros where time would go in 
# the dimentions should be determined by how the input is given
# get following lines of times and assign them to an time array/list

for line in eachline(stdin) 
	# do equation stuff here
    #read line, split, assign each value
    times = line
    t_array = split(times, " ")
    T_i = parse(Float64, t_array[1])
    T_j = parse(Float64, t_array[2])
    T_k = parse(Float64, t_array[3])
    T_l = parse(Float64, t_array[4])

    #=
	The distance between sattelites can be calculated by multiplying the difference between
	signal arrival times by the speed of light. Assuming all sattelites sent their signals at the same time,
	the time difference can only be explained by one being closer or further.
	
	This is the R value in the equation set.
	=#
	R_i = (T_i * SPEED_OF_LIGHT)
	R_j = (T_j * SPEED_OF_LIGHT)
	R_k = (T_k * SPEED_OF_LIGHT)
	R_l = (T_l * SPEED_OF_LIGHT)

    R_ij = R_i - R_j
    R_ik = R_i - R_k
    R_kj = R_k - R_j
    R_kl = R_k - R_l

    # set of quantities - Equations 1 -> 6 on the paper
    x_ijy = (R_ij * y_ki) - (R_ik * y_ji)
    x_ikx = (R_ik * x_ji) - (R_ij * x_ki)
    x_ikz = (R_ik * z_ji) - (R_ij * z_ki)
    x_kjy = (R_kj * y_lk) - (R_kl * y_jk)
    x_klx = (R_kl * x_jk) - (R_kj * x_lk)
    x_klz = (R_kl * z_jk) - (R_kj * z_lk)

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
	A = x_ikx / x_ijy
	B = x_ikz / x_ijy
	C = x_klx / x_kjy
	D = x_klz / x_kjy
	
	print(A)
	print(B)
	print(C)
	print(D)
	# I don't even know what's a variable anymore and what's merely part of another var
	E = ((R_ik * R_ij2) - (R_ij * R_ik2)) / (2 * x_ijy)
	F = ((R_kl * R_kj2) - (R_kj * R_kl2)) / (2 * x_kjy)
	
	print(E)
	print(F)
	# more following equations
	G = (D - B) / (A - C)
	H = (F - E) / (A - C)
	
	print(G)
	print(H)
	# I and J I guess
	I = (A * G) + B
	J = (A * H) + E
	
	print(I)
	print(J)

	# Once more I am truly unsure of what is and isn't a variabele. Why are these only using ki?
	K = R_ik2 + (2 * x_ki * H) + (2 * y_ki * J)
	L = 2 * ((x_ki * G) + (y_ki * I) + z_ki)
	
	print(K)
	print(L)
	# Almost there...
	M = ((4 * (R_ik ^ 2)) * ((G ^ 2) + (I ^ 2) + 1)) - (L ^ 2)
	N = ((8 * (R_ik ^ 2)) * (G * (i_coords[1] - H) + (I * (i_coords[2] - J)) + i_coords[3])) + (2 * L * K)
	O = ((4 * (R_ik ^ 2)) * (((i_coords[1] - H) ^ 2) + ((i_coords[2] - J) ^ 2) + (i_coords[3] ^ 2))) - (K ^ 2)
	
	print(M)
	print(N)
	print(O)
	
	
	# Finally, a formula I recognize
	Q = N + ((sgn(N))*(sqrt((N^2) - (4 * M * O))))

	print(Q)
	
	# Finally, The coordinates.
	z_pos = (2 * O) / Q_pos - 2
	z_neg = (2 * O) / Q_neg - 2
	x_pos = G * z_pos + H
	x_neg = G * z_neg + H
	x_pos = I * z_pos + J
	x_neg = I * z_neg + J
	
	# and the distances from the origin for some reason?
	pos_distance = (x_pos ^ 2 + y_pos ^ 2 + z_pos ^ 2) ^ 0.5
	neg_distance = (x_neg ^ 2 + y_neg ^ 2 + z_neg ^ 2) ^ 0.5
	
	@printf "g= %9.2e, h= %9.2e, j= %9.2e, m= %9.2e, o= %9.2e \n" G H J M O
	@printf "+) x= %10i, y= %10i, z= %10i; I= %10i%n" x_pos y_pos z_pos pos_distance
	@printf "+) x= %10i, y= %10i, z= %10i; I= %10i%n" x_neg y_neg z_neg neg_distance
end