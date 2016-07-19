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