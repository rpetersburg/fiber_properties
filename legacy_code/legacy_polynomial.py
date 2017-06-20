mesh_grid = mesh_grid_from_array(image)
if x0 is None or y0 is None:
    x0 = image.shape[1] / 2.0
    y0 = image.shape[0] / 2.0
# Place (0, 0) at fiber center
mesh_grid[0] -= x0
mesh_grid[1] -= y0

initial_guess = tuple(np.ones(deg/2 + 1))

opt_coeffs, cov_matrix = opt.curve_fit(polynomial_array,
                                       mesh_grid,
                                       image.ravel(),
                                       p0=initial_guess)

return polynomial_array(mesh_grid, *opt_coeffs).reshape(image.shape)

x_array = mesh_grid[0]
y_array = mesh_grid[1]
r_array = x_array**2 + y_array**2

polynomial_array = np.zeros_like(r_array)
for i in xrange(len(coeff)):
    polynomial_array += coeff[i] * r_array**i