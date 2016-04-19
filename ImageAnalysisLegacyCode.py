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