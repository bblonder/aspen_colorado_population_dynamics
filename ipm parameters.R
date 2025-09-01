# biggest possible size
DBH_min = 5
DBH_max = 60
MAX_ITERATIONS = 300
# IPM size resolution
mesh_points = DBH_max

# number of resamples for uncertainty propagation in data
NUM_RESAMPLES = 10

dbh_range = seq(0, DBH_max, length.out=mesh_points)

# cores for parallel
NUM_CORES = 4
