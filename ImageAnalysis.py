import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from PIL import Image

class ImageAnalysis():

  def __init__(self, imageString):
    self.threshold = 1

    self.imageString = imageString
    self.imageArray = self.convertImageToArray(imageString)
    self.darkImageArray = None
    self.flatFieldImageArray = None
    self.backlitImageArray = None
    
    self.height, self.width = self.imageArray.shape

    self.left = None
    self.right = None
    self.top = None
    self.bottom = None
    self.fiberDiameter = None

    self.centroidY = None
    self.centroidX = None
    self.centerY_circle = None
    self.centerX_circle = None
    self.centerY_edge = None
    self.centerX_edge = None

    # Golden Ratio for maximization tests
    self.phi = (5 ** 0.5 - 1) / 2

#=============================================================================#
#==== Image Initialization ===================================================#
#=============================================================================#

  def convertImageToArray(self, imageString):
    return np.array(Image.open(imageString))[:,:,0].astype(float)

  def setDarkImage(self, *imageStrings):
    if self.darkImageArray is not None:
      print 'Dark Images have already been set'
      return

    self.darkImageArray = np.zeros((self.height, self.width))
    for imageString in imageStrings:
      self.darkImageArray += self.convertImageToArray(imageString) / float(len(imageStrings))

    # Remove any dark current from each pixel
    self.imageArray -= self.darkImageArray

    # 
    for x in xrange(self.width):
      for y in xrange(self.height):
        if self.imageArray[y,x] < 0:
          self.imageArray[y,x] = 0.0

  def setFlatFieldImage(self, *imageStrings):
    if self.flatFieldImageArray is not None:
      print 'Flat Field Images have alraedy been set'
      return

    self.flatFieldImageArray = np.zeros((self.height, self.width))
    for imageString in imageStrings:
      self.flatFieldImageArray += self.convertImageToArray(imageString) / float(len(imageStrings))

  def executeErrorCorrections(self, imageArray = None):
    """
      executeErrorCorrections(imageArray = None)

      Uses Flat Field Image and Dark Image to correct for errors in the detector.

      returns: correctedImageArray (argument given)
      returns: None (no argument given)
    """
    if self.flatFieldImageArray is None or self.darkImageArray is None:
      if self.flatFieldImageArray is None:
        print 'Flat Field Image has not been set'
      if self.darkImageArray is None:
        print 'Dark Image has not been set'
      return

    if imageArray is None:
      self.imageArray = self.removeDarkImage(self.imageArray)
      self.imageArray = self.imageArray * np.mean(self.flatFieldImageArray) / (self.flatFieldImageArray - self.darkImageArray)
      return None
    else:
      imageArray = self.removeDarkImage(imageArray)
      return imageArray * np.mean(self.flatFieldImageArray) / (self.flatFieldImageArray - self.darkImageArray)

  def removeDarkImage(self, imageArray):
    imageArray -= self.darkImageArray

    # Prevent any pixels from becoming negative values
    for x in xrange(self.width):
      for y in xrange(self.height):
        if imageArray[y,x] < 0:
          imageArray[y,x] = 0.0

    return imageArray

  def setBacklitImage(self, imageString):
    self.backlitImageArray = self.convertImageToArray(imageString)

#=============================================================================#
#==== Image Analysis =========================================================#
#=============================================================================#

  def getFiberCentroid(self):
    """
      getFiberCentroid()

      returns: (centroidY, centroidX)
    """
    if self.centroidY is None or self.centroidX is None:
      rowSum = self.getRowSum()
      columnSum = self.getColumnSum()

      rowWeightedSum = 0
      rowWeight = 0
      for i in xrange(self.width):
        rowWeightedSum += i * rowSum[i]
        rowWeight += rowSum[i]
      centroidColumn = rowWeightedSum/rowWeight

      columnWeightedSum = 0
      columnWeight = 0
      for i in xrange(self.height):
        columnWeightedSum += i * columnSum[i]
        columnWeight += columnSum[i]
      centroidRow = columnWeightedSum/columnWeight

      self.centroidY = centroidRow
      self.centroidX = centroidColumn

    return self.centroidY, self.centroidX

  def getRowSum(self):
    rowSum = np.sum(self.imageArray, axis=0)
    return ((rowSum - np.min(rowSum)) / self.height).astype(float)

  def getColumnSum(self):
    columnSum = np.sum(self.imageArray, axis=1)
    return ((columnSum - np.min(columnSum)) / self.width).astype(float)

  def getFiberCenterCircleIteration(self, radius, imageArray = None):
    """
      getFiberCenterCircleIteration(radius)

      Uses golden mean method to find the minimum total pixel sum with a removed
        circle with radius centered at (x,y)

      returns: (centerY, centerX)
    """
    if imageArray is None:
      imageArray = self.imageArray

    # Create four "corners" to test center of the removed circle
    x = np.zeros(4).astype(float)
    x[0] = radius
    x[3] = self.width - radius
    x[1] = x[0] + (1 - self.phi) * (x[3] - x[0])
    x[2] = x[0] + self.phi * (x[3] - x[0])

    y = np.zeros(4).astype(float)
    y[0] = radius
    y[3] = self.height - radius
    y[1] = y[0] + (1 - self.phi) * (y[3] - y[0])
    y[2] = y[0] + self.phi * (y[3] - y[0])

    # Initialize the array sum to each corner
    arraySum = np.zeros((2,2)).astype(float)
    for i in xrange(2):
      for j in xrange(2):
        # Center at (y[j+1], x[i+1])
        removedCircleArray = self.removeCircle(radius, x[i+1], y[j+1], imageArray)
        arraySum[j, i] = self.getArraySum(removedCircleArray)
        print (j+1, i+1), (y[j+1], x[i+1]), arraySum[j, i]

    # Find the index of the minimum array sum corner
    minIndex = np.unravel_index(np.argmin(arraySum), (2,2))
    print arraySum[minIndex], (y[minIndex[0]+1], x[minIndex[1]+1]), (minIndex[0]+1, minIndex[1]+1)

    # Move the other corners to smaller search area
    if minIndex[0] == 0:
      y[3] = y[2]
      y[2] = y[1]
      y[1] = y[0] + (1 - self.phi) * (y[3] - y[0])
    else:
      y[0] = y[1]
      y[1] = y[2]
      y[2] = y[0] + self.phi * (y[3] - y[0])
    if minIndex[1] == 0:
      x[3] = x[2]
      x[2] = x[1]
      x[1] = x[0] + (1 - self.phi) * (x[3] - x[0])
    else: 
      x[0] = x[1]
      x[1] = x[2]
      x[2] = x[0] + self.phi * (x[3] - x[0])

    while x[3] - x[0] > 2 or y[3] - y[0] > 2:
      # Replace the opposite corner array sum 
      arraySum[1 - minIndex[1], 1 - minIndex[0]] = arraySum[minIndex]
      print (1 - minIndex[1] + 1, 1-minIndex[0] + 1), (y[1-minIndex[0]+1], x[1-minIndex[1]+1]), arraySum[1 - minIndex[1], 1 - minIndex[0]]

      # Recalculate new sums for other three corners
      for i in xrange(2):
        for j in xrange(2):
          if i != 1 - minIndex[1] or j != 1 - minIndex[0]:
            removedCircleArray = self.removeCircle(radius, x[i+1], y[j+1], imageArray)
            arraySum[j, i] = self.getArraySum(removedCircleArray)
            print (j+1,i+1), (y[j+1], x[i+1]), arraySum[j, i]

      minIndex = np.unravel_index(np.argmin(arraySum), (2,2)) # Tuple
      print arraySum[minIndex], (y[minIndex[0]+1], x[minIndex[1]+1]), (minIndex[0]+1, minIndex[1]+1)

      if minIndex[0] == 0:
        y[3] = y[2]
        y[2] = y[1]
        y[1] = y[0] + (1 - self.phi) * (y[3] - y[0])
      else:
        y[0] = y[1]
        y[1] = y[2]
        y[2] = y[0] + self.phi * (y[3] - y[0])
      if minIndex[1] == 0:
        x[3] = x[2]
        x[2] = x[1]
        x[1] = x[0] + (1 - self.phi) * (x[3] - x[0])
      else: 
        x[0] = x[1]
        x[1] = x[2]
        x[2] = x[0] + self.phi * (x[3] - x[0])  

    minCircleSum = self.getArraySum(imageArray)
    for xc in xrange(x[0], x[3] + 1):
      for yc in xrange(y[0], y[3] + 1):
        removedCircleArray = self.removeCircle(radius, xc, yc)
        circleSum = self.getArraySum(removedCircleArray)
        print xc, yc, circleSum
        if circleSum < minCircleSum:
          centerX = xc
          centerY = yc
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

  def removeCircle(self, radius, x, y, imageArray = None):
    """
      removeCircle(imageArray, radius, x, y)

      returns: Copy of 2D numpy array with circle removed from image
    """
    if imageArray is None:
      imageArray = self.imageArray

    # Round center to nearest pixel... for now...
    x = int(round(x))
    y = int(round(y))

    outputArray = np.copy(imageArray)
    for i in xrange(x - radius, x + radius + 1):
      for j in xrange(y - radius, y + radius + 1):
        if (x-i)**2 + (y-j)**2 <= radius**2:
          outputArray[j,i] = 0.0
    return outputArray

  def getArraySum(self, imageArray):
    return np.sum(imageArray)

  def circleArray(self, radius, x, y):
    """
      circleArray(radius, x, y)

      returns: 2D numpy array of integer values where points
        inside the circle are 1 and outside the circle are 0
    """
    circleArray = np.zeros((self.height, self.width))
    for i in xrange(x - radius, x + radius + 1):
      for j in xrange(y - radius, y + radius, + 1):
        if (x-i)**2 + (y-j)**2 <= radius**2:
          circleArray[j,i] = 1

  def getFiberCenterEdgeMethod(self):
    """
      getFiberCenterEdgeMethod( )

      The averages of the fiber edges gives the fiber center

      returns: centerY [float], centerX [float]
    """
    if self.centerY_edge is None or self.centerX_edge is None:      
      if self.left is None or self.top is None:
          self.setFiberEdges()

      self.centerY_edge = (self.top + self.bottom)/2.0
      self.centerX_edge = (self.left + self.right)/2.0

    return self.centerY_edge, self.centerX_edge

  def getFiberDiameter(self):
    if self.fiberDiameter is None:
      self.setFiberEdges()
    return self.fiberDiameter

  def setFiberEdges(self, threshold = None):
    """
      setFiberEdges( threshold = None )

      Sets the left, right, top, and bottom edges of the fiber by finding where
        the row and column sums cross the given threshold. Also sets the width
        of the fiber by the maximum of the horizontal and vertical lengths

      returns: None
    """
    if threshold is None:
      threshold = self.threshold

    rowSum = self.getRowSum()
    columnSum = self.getColumnSum()

    left = -1
    right = -1
    for index, value in enumerate(rowSum):
      if left < 0:
        if value > threshold:
          left = index
      else:
        if value > threshold:
          right = index
    top = -1
    bottom = -1
    for index, value in enumerate(columnSum):
      if top < 0:
        if value > threshold:
          top = index
      else:
        if value > threshold:
          bottom = index

    self.left = left
    self.right = right
    self.top = top
    self.bottom = bottom
    self.fiberDiameter = max(right - left, bottom - top)

#=============================================================================#
#===== Image Plotting ========================================================#
#=============================================================================#

  def plotHorizontalCrossSection(self, row):
    plt.plot(self.imageArray[row,:])
    plt.title('Horizontal Cross Section (row = %s)'%row)
    plt.xlabel('Pixel')

  def plotVerticalCrossSection(self, column):
    plt.plot(self.imageArray[:,column])
    plt.title('Vertical Cross Section (column = %s)'%column)
    plt.xlabel('Pixel')

  def plotCrossSections(self, row, column):
    plt.figure(1)
    plt.subplot(211)
    self.plotHorizontalCrossSection(row)
    plt.subplot(212)
    self.plotVerticalCrossSection(column)
    plt.show()

  def plotCrossSectionSums(self):
    plt.figure(1)
    plt.subplot(211)
    plt.plot(self.getRowSum())
    plt.title('Average for each Column')
    plt.xlabel('Column')
    plt.subplot(212)
    plt.plot(self.getColumnSum())
    plt.title('Average for each Row')
    plt.xlabel('Row')
    plt.show()

  def showImageArray(self, imageArray=None):
    if imageArray is None:
        imageArray = self.imageArray
    plt.figure(1)
    plt.imshow(imageArray)
    plt.show()

  def getVariation(self):
    centerRow, centerColumn = self.getFiberCenter()

  def iterationCounter(self, iteration):
    iteration += 1
    if iteration % 10 == 0:
        print iteration
    return iteration