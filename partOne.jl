using Printf
i = readline()
j = readline()
k = readline()
l = readline()

# we need to split the coordinates from the readlines
# so we can get each individual x, y, z value
i_array = split(i, " ")
j_array = split(j, " ")
#k_array = split(k, " ")
#l_array = split(l, " ")

# x four distances
x_ji = parse(Float64, j_array[1]) - parse(Float64, i_array[1])
x_ki = parse(Float64, k_array[1]) - parse(Float64, i_array[1])
x_jk = parse(Float64, j_array[1]) - parse(Float64, k_array[1])
x_lk = parse(Float64, l_array[1]) - parse(Float64, k_array[1])

println(x_ji)
# y four distances 
y_ki = parse(Float64, k_array[2]) - parse(Float64, i_array[2])
y_ji = parse(Float64, j_array[2]) - parse(Float64, i_array[2])
y_lk = parse(Float64, l_array[2]) - parse(Float64, k_array[2])
y_jk = parse(Float64, j_array[2]) - parse(Float64, k_array[2])


# z four distances
z_ji = parse(Float64, j_array[3]) - parse(Float64, i_array[3])
z_ki = parse(Float64, k_array[3]) - parse(Float64, i_array[3])
z_jk = parse(Float64, j_array[3]) - parse(Float64, k_array[3])
z_lk = parse(Float64, l_array[3]) - parse(Float64, k_array[3])

# Speed of light (m/s)
c = 299792458

# for now, im initiating an array of zeros where time would go in 
# the dimentions should be determined by how the input is given
# get following lines of times and assign them to an time array/list

for line in eachline(stdin) 
	# do equation stuff here
    #read line, split, assign each value
    times = readline()
    t_array = split(times, " ")
    T_i = parse(Floart64, t_array[1])
    T_j = parse(Float64, t_arrat[2])
    T_k = parse(Float64, t_array[3])
    T_l = parse(Float64, t_array[4])

    # R values
    R_ij = T_i - T_j
    R_ik = T_i - T_k
    R_kj = T_k - T_j
    R_kl = T_k - T_l

    # set of quantities - Equations 1 -> 6 on the paper
    X_ijy = (R_ij * y_ki) - (R_ik * y_ji)
    X_ikx = (R_ik * x_ji) - (R_ij * x_ki)
    X_ikz = (R_ik * z_ji) - (R_ij * z_ki)
    X_kjy = (R_kj * y_lk) - (R_kl * y_jk)
    X_klx = (R_kl * x_jk) - (R_kj * x_lk)
    X_klz = (R_kl * z_jk) - (R_kj * z_lk)

	g = -952
	h = 401000000
	j = -22900000
	m = 25500000000000000000
	o = 446000000000000000000000000000000
	@printf "g= %9.2e, h= %9.2e, j= %9.2e, m= %9.2e, o= %9.2e \n" g h j m o
	@printf "+) x= %10i, y= %10i, z= %10i; I= %10i%n" pos_objectX pos_objectY pos_objectZ pos_distance
	@printf "+) x= %10i, y= %10i, z= %10i; I= %10i%n" neg_objectX neg_objectY neg_objectZ neg_distance
end
#print("next")
#print(coords1)
#print(coords2)
#print(coords3)
#print(coords4)

