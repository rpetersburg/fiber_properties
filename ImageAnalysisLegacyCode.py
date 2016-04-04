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