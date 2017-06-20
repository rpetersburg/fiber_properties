    def getFiberCenterEdgeDetection(self):
        kernel = np.array([[0,1,0],
                           [1,0,1],
                           [0,1,0]])

        #edgeImage = ndi.convolve(self.imageArray, kernel, mode='constant', cval=0.0)
        gausImage = ndi.filters.gaussian_filter(self.imageArray, 50, mode='nearest')

        plt.subplot(211)
        plt.imshow(self.imageArray)
        plt.subplot(212)
        plt.imshow(gausImage)
        plt.show()


    def getFiberCenterCircleIteration(self):

      centerX, centerY = self.getMaxCircleCenter(radius, radius, self.width - radius, radius, self.height - radius, 64)

      # Gradually decrease the range of x & y and the iteration delta to hone in on center
      for interval, delta in [(64,32),(32,16),(16,8),(8,4),(4,2),(2,1)]:
        centerY, centerX = self.getMaxCircleCenter(radius, centerX - interval, centerX + interval,
                                                           centerY - interval, centerY + interval, delta)
      self.centerY_circle = centerY
      self.centerX_circle = centerX


    def getMaxCircleCenter(self, radius, xStart, xEnd, yStart, yEnd, delta):
        maxCircleSum = -1
        centerX = -1
        centerY = -1

        # Prevent problems at the edges of the image array
        if xStart - radius < 0:
          xStart = radius
        if xEnd + radius >= self.width:
          xEnd = self.width - radius - 1
        if yStart - radius < 0:
          yStart = radius
        if yEnd + radius >= self.height:
          yEnd = self.height - radius - 1

        # Construct a circle with center (x,y) and radius and sum over contained pixels
        for x in xrange(xStart, xEnd + 1, delta):
          for y in xrange(yStart, yEnd + 1, delta):
            circleSum = 0
            for i in xrange(x - radius, x + radius + 1):
              for j in xrange(y - radius, y + radius + 1):
                if (x-i)**2 + (y-j)**2 <= radius**2:
                  circleSum += self.imageArray[j,i].astype(float)
            circleSum /= radius**2 # Normalization to the area
            if circleSum > maxCircleSum:
              maxCircleSum = circleSum
              centerX = x
              centerY = y

        return centerY, centerX


  def getFiberCenterCircleIteration(self, radius):
    """
      getFiberCenterCircleIteration(radius)

      Uses golden mean method to maximize pixel sum over various circle centers

      returns: (centerY, centerX)
    """
    x = np.zeros(4).astype(int)
    x[0] = radius
    x[1] = int(round(x0 + (1 - self.phi) * (x3 - x0)))
    x[2] = int(round(x0 + self.phi * (x3 - x0)))
    x[3] = self.width - radius

    y = np.zeros(4).astype(int)
    y[0] = radius
    y[1] = int(round(y0 + (1 - self.phi) * (y3 - y0)))
    y[2] = int(round(y0 + self.phi * (y3 - y0)))
    y[3] = self.height - radius

    circleSum = np.zeros((2,2)).astype(int)
    for i in xrange(1,3):
      for j in xrange(1,3):
        circleSum[j,i] = self.getCircleSum(radius, x[i], y[j])

    while x3 - x0 > 2 or y3 - y0 > 2:
      maxIndex = np.unravel_index(np.argmax(circleSum), (2,2))

      if maxIndex[0] == 0:
        y3 = y2
        y2 = y1
        y1 = int(round(y0 + (1 - self.phi) * (y3 - y0)))
      else:
        y0 = y1
        y1 = y2
        y2 = y2 = int(round(y0 + self.phi * (y3 - y0)))  
      if maxIndex[1] == 0:
        x3 = x2
        x2 = x1
        x1 = int(round(x0 + (1 - self.phi) * (x3 - x0)))
      else: 
        x0 = x1
        x1 = x2
        x2 = int(round(x0 + self.phi * (x3 - x0)))

      # Replace the opposite corner circle sum
      circleSum[int(not bool(maxIndex))] = circleSum[maxIndex]
      for i in xrange(1,3):
        for j in xrange(1,3):
          if (j,i) != int(not bool(maxIndex))
   

    maxCircleSum = -1
    for x in xrange(x0, x3 + 1):
      for y in xrange(y0, y3 + 1):
        circleSum = self.getCircleSum(radius, x, y)
        print x,y,circleSum
        if circleSum > maxCircleSum:
          centerX = x
          centerY = y
          maxCircleSum = circleSum

    return centerY, centerX



def getCircleSum(self, radius, x, y, imageArray = None):
    """
      getCircleSum(radius, x, y)

      Sums all pixels within circle with radius centered at pixel (x,y) in self.ImageArray

      returns: circleSum
    """
    if imageArray is None:
      imageArray = self.imageArray

    circleSum = 0.0
    for i in xrange(x - radius, x + radius + 1):
      for j in xrange(y - radius, y + radius + 1):
        if (x-i)**2 + (y-j)**2 <= radius**2:
          circleSum += imageArray[j,i].astype(float)
    return circleSum


def setFiberCenterTophatMethod(self):
        """Finds fiber center using a Tophat Fit
        
        Sets:
            fiber_diameter_tophat
            center_y_tophat
            center_x_tophat
        """
        approx_center_x = self.getFiberCenter()[1]
        approx_center_y = self.getFiberCenter()[0]
        approx_radius = self.getFiberRadius()

        initial_guess = (approx_center_x, approx_center_y, approx_radius)

        bounds = [(approx_center_x - 10, approx_center_x + 10),
                  (approx_center_y - 10, approx_center_y + 10),
                  (approx_radius - 10, approx_radius + 10)]

        opt_parameters = opt.differential_evolution(self.removedCircleArraySum, bounds=bounds, tol=1).x

        self.tophat_fit = self.circleArray((self.mesh_grid), *opt_parameters) * 256

        self.center_x_tophat = opt_parameters[0]
        self.center_y_tophat = opt_parameters[1]
        self.fiber_diameter_tophat = opt_parameters[2] * 2

        #self.showImageArray(self.tophat_fit)
        self.plotOverlaidCrossSections(self.image_array, self.tophat_fit,
                                       self.center_y_tophat, self.center_x_tophat)

def removedCircleArraySum(self, var):
        """
        Args:
            var[0]: x0
            var[1]: y0
            var[2]: radius
        """
        return self.getArraySum(self.removeCircle(var[2], var[0], var[1], self.image_array)) + np.pi * var[2]**2

def tophatArray(self, mesh_grid, x0, y0, radius, amp):
        """Creates a 2D tophat function as a 1D array

        Args:
            mesh_grid: independent variables x and y separated into two arrays
                each with the proper dimensions (np.meshgrid)
            x0: center position x of tophat
            y0: center position y of tophat
            radius: radius of tophat
            amp: amplitude of tophat

        Returns:
            tophat_array: Ravelled tophat numpy array (single dimension) usable in
                Scipy.optimize.curve_fit method and properly reshaped by
                tophat_array.reshape(height, width)
        """
        res = 1 # Resolution element for more precise circular edges

        x_array = mesh_grid[0]
        y_array = mesh_grid[1]
        x0 = float(x0)
        y0 = float(y0)
        radius = float(radius)
        amp = float(amp)

        if res == 1:
            tophat_array = amp * ((x_array-x0)**2 + (y_array-y0)**2 <= (radius * 1.001)**2).astype(float)
        else:
            tophat_array = amp * ((x_array-x0)**2 + (y_array-y0)**2 < (radius-1)**2).astype(float)

            height, width = tophat_array.shape
            for x in xrange(width):
                for y in xrange(height):
                    if (x-x0)**2 + (y-y0)**2 >= (radius-1)**2 and (x-x0)**2 + (y-y0)**2 <= (radius+1)**2:
                        for i in xrange(res):
                            for j in xrange(res):
                                x_temp = x - 0.5 + 1 / (2.0*res) + i / float(res)
                                y_temp = y - 0.5 + 1 / (2.0*res) + j / float(res)
                                if (x_temp-x0)**2 + (y_temp-y0)**2 < radius**2:
                                    tophat_array[y, x] += amp / res**2

        return tophat_array

def getFiberCenterTophatMethod(self):
    if self.center_x_tophat is None:
        self.setFiberCenterTophatMethod()
    return self.center_y_tophat, self.center_x_tophat

    
def getFiberDiameterTophatMethod(self):
    if self.fiber_diameter_tophat is None:
        self.setFiberCenterTophatMethod()            
    return self.fiber_diameter_tophat



def polynomialFit(self, image_array=None):
    """Creates a polynomial fit for the image

    Sets:
        polynomial_fit: array of np.poly1d objects with 
    """
    if image_array is None:
        image_array = self.image_array

    height, width = image_array.shape

    y_array = np.arange(height)
    coeffs = np.polyfit(y_array, image_array, 103)

    polynomial_func = []
    for i in xrange(width):
        polynomial_func.append(np.poly1d(coeffs[:,i]))

    polynomial_fit = np.zeros_like(image_array)

    for i in xrange(width):
        polynomial_fit[:,i] = polynomial_func[i](y_array)

    self.showImageArray(polynomial_fit)
    self.plotOverlaidCrossSections(image_array, polynomial_fit, *self.getFiberCenter())

def getPolynomialFit2(self, image_array, deg=4):
    height, width = image_array.shape
    mesh_grid = self.getMeshGrid(image_array)
    x_array = mesh_grid[0].astype('float64')
    y_array = mesh_grid[1].astype('float64')
    num_terms = (deg+1)**2

    poly_tensor = np.zeros((height, width, num_terms))
    ij = itertools.product(range(deg+1), range(deg+1))
    for k, (i,j) in enumerate(ij):
        poly_tensor[:,:,k] = x_array**i * y_array**j

    opt_coeffs, resids, rk, s = np.linalg.lstsq(poly_tensor.reshape(num_terms, height*width).swapaxes(0,1),
                                                image_array.reshape(height*width))

    return self.polynomialArray(mesh_grid, *opt_coeffs).reshape(image_array.shape)


def setBacklitImage(self, *images):
      """Sets the backlit fiber image

      Args:
          *images: either a single numpy array or multiple image file names

      Returns:
          None
      """
      if isinstance(images[0], basestring):
          self.backlit_image_array = np.zeros((self.image_height, self.image_width))

          for image_string in images:
              self.backlit_image_array += self.convertImageToArray(image_string) \
                                          / float(len(images))
      else:
          self.backlit_image_array = images[0]


# Depricated centroiding method

    x_array = self.getRowSum(x_array)
    y_array = self.getColumnSum(y_array)

    row_sum = self.getRowSum(image_array_iso)
    column_sum = self.getColumnSum(image_array_iso)

    row_weighted_sum = 0
    row_weight = 0
    for i in xrange(self._image_width):
        row_weighted_sum += i * row_sum[i]
        row_weight += row_sum[i]
    centroid_column = row_weighted_sum/row_weight

    column_weighted_sum = 0
    column_weight = 0
    for i in xrange(self._image_height):
        column_weighted_sum += i * column_sum[i]
        column_weight += column_sum[i]
    centroid_row = column_weighted_sum / column_weight

    self._centroid_y = centroid_row
    self._centroid_x = centroid_column

def applyWindow():
  fft_window_x, fft_window_y = np.meshgrid(self.hann_poisson_window(width),
                                           self.hann_poisson_window(height)) 
  self.showImageArray(fft_window_x*fft_window_y)     
  self.plotOverlaidCrossSections(image_array,
                                 fft_window_x * fft_window_y * image_array.max(), 
                                 height/2, width/2)

  return image_array * fft_window_x * fft_window_y

# Alternate circle array method (too slow)
height, width = mesh_grid[0].shape

rad_ceil = int(np.ceil(radius)) + 1
res_array = np.arange(-0.5, 2 * rad_ceil - 0.5, 1.0 / res) + 0.5 / res
res_mesh_x, res_mesh_y = np.meshgrid(res_array, res_array)
x0_new = rad_ceil + x0 - int(round(x0))
y0_new = rad_ceil + y0 - int(round(y0))
circle_array = ((res_mesh_x-x0_new)**2 + (res_mesh_y-y0_new)**2 <= radius**2).astype('float64') / res**2
circle_array = NumpyArrayHandler.sum_chunk(NumpyArrayHandler.sum_chunk(circle_array, res, axis=0), res, axis=1)
top_pad = int(round(y0)) - rad_ceil
bottom_pad = height - (top_pad + 2 * rad_ceil)
left_pad = int(round(x0)) - rad_ceil
right_pad = width - (left_pad + 2 * rad_ceil)
circle_array = np.pad(circle_array, ((top_pad, bottom_pad), (left_pad, right_pad)), mode='constant', constant_values=0)


@staticmethod
def sum_chunk(x, chunk_size, axis=-1):
    shape = x.shape
    if axis < 0: 
        axis += x.ndim
    shape = shape[:axis] + (-1, chunk_size) + shape[axis+1:]
    x = x.reshape(shape)
    return x.sum(axis=axis+1)

def getFiberCenterRadiusMethod(self, tol=1, test_range=None, show_image=False):
    """Getter for the fiber center using the radius method

    See setFiberCenterRadiusMethod() for method details

    Args:
        show_image (optional): boolean for whether or not to show image of
            completed method

    Returns:
        center y (pixels), center x (pixels)
    """
    if self._center['radius']['x'] is None:
        self.setFiberCenterRadiusMethod(tol, test_range)
        
    if show_image:
        self.showOverlaidTophat(self._center['radius']['x'],
                                self._center['radius']['y'],
                                self._diameter['radius'] / 2.0,
                                tol=1)

    return self._center['radius']['y'], self._center['radius']['x']

def getFiberCenterCircleMethod(self, radius=None, tol=1, test_range=None, show_image=False):
    """Getter for the fiber center using the circle method

    See setFiberCenterCircleMethod() for method details

    Args:
        show_image (optional): boolean for whether or not to show image of
            completed method

    Returns:
        center y (pixels), center x (pixels)
    """
    if radius is None:
        radius = self.getFiberRadius()
    if self._center['circle']['x'] is not None:
        show_image = False

    self.setFiberCenterCircleMethod(radius, tol, test_range)

    if show_image:
        self.showOverlaidTophat(self._center['circle']['x'],
                                self._center['circle']['y'],
                                self._diameter['circle'] / 2.0,
                                tol=1)

    return self._center['circle']['y'], self._center['circle']['x']

def getFiberCenterEdgeMethod(self, show_image=False):
    """Getter for the fiber center using the edge method

    See setFiberCenterEdgeMethod() for method details

    Args:
        show_image (optional): boolean for whether or not to show image of
            completed method

    Returns:
        center y (pixels), center x (pixels)
    """
    if self._center['edge']['y'] is None:
        self.setFiberCenterEdgeMethod()

    if show_image:
        self.showOverlaidTophat(self._center['edge']['x'],
                                self._center['edge']['y'],
                                self._diameter['edge'] / 2.0)

    return self._center['edge']['y'], self._center['edge']['x']

def getFiberCenterGaussianMethod(self, show_image=False):
    """Getter for the fiber center using the gaussian method

    See setFiberCenterGaussianMethod() for method details

    Args:
        show_image (optional): boolean for whether or not to show image of
            completed method

    Returns:
        center y (pixels), center x (pixels)
    """
    if self._center['gaussian']['x'] is None:
        self.setFiberCenterGaussianMethod()

    if show_image:
        showImageArray(self._fit['gaussian'])
        plotOverlaidCrossSections(self.image, self._fit['gaussian'],
                                  self._center['gaussian']['y'], self._center['gaussian']['x'])

    return self._center['gaussian']['y'], self._center['gaussian']['x']


def setImageArray(self, image_input, calibration):
  """Sets image to be analyzed

  Args
  ----
  image_input : {None, 1D iterable, 2D iterable, string}
      The input used to set the image array. See
      NumpyArrayHandler.convertImageToArray() for details

  Sets
  ----
  image
  """
  uncorrected_image, output_obj = convertImageToArray(image_input, True)
  if uncorrected_image is None:
      return

  self.setImageInfo(output_obj)

  if calibration is None:
      calibration = Calibration(None, None, None)

  self.image = calibration.executeErrorCorrections(uncorrected_image,
                                                   self._image_info.subframe_x,
                                                   self._image_info.subframe_y,
                                                   self._image_info.exp_time)


class FRD_Input(object):
    """Container for FRD test input information

    Attributes
    ----------
    name : string
        Name to be used for saving plots
    input_files : list(list(string)) or list(string)
        Far field image file names grouped by input f-number
        e.g. [['in_2.5/im_000.fit', 'in_2.5/im_001.fit', ...],
              ['in_3.0/im_000.fit', 'in_3.0/im_001.fit', ...], ...]
        or   ['in_2.5/im_data.pkl', in_3.0/im_data.pkl', ...]
    output_files : list(list(string))
        Far field image file names grouped by output f_number
        e.g. [['out_2.5/im_000.fit', 'out_2.5/im_001.fit', ...],
              ['out_3.0/im_000.fit', 'out_3.0/im_001.fit', ...], ...]
        or   ['out_2.5/im_data.pkl', out_3.0/im_data.pkl', ...]
    input_fnums : list(float), optional
        Input focal ratios which have associated images. Must have same length
        as input_files
    cal_fnums : list(float), optional
        Output focal ratios which were used as calibration images. Must have
        same length as output_files

    Raises
    ------
    RuntimeError
        If the list lengths do not match up as described above
    """
    def __init__(self, name, input_files, output_files,
                 dark=None, ambient=None, flat=None,
                 input_fnums=[2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
                 cal_fnums=[3.0, 4.0, 5.0]):
        self.name = name
        self.input_files = input_files
        self.output_files = output_files
        self.dark = dark
        self.ambient = ambient
        self.flat = flat
        self.input_fnums = input_fnums
        self.cal_fnums = cal_fnums

        if len(input_files) != len(input_fnums):
            raise RuntimeError('Number of input files and input f numbers must match')
        if len(output_files) != len(cal_fnums):
            raise RuntimeError('Number of output files and calibration f numbers must match')