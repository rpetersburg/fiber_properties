import numpy as np
import ImageAnalysis as IA

class FiberProperties():

  def __init__(self, inObject = None, nfObject = None, ffObject = None):
    self.inObject = inObject
    self.nfObject = nfObject
    self.ffObject = ffObject
    self.nearFieldScramblingGain = None
    self.farFieldScramblingGain = None

  def setInputObject(self, inObject):
    self.inObject = inObject

  def setNearFieldObject(self, nfObject):
    self.nfObject = nfObject

  def setFarFieldObject(self, ffObject):
    self.ffObject = ffObject

  def getNearFieldScramblingGain(self):
    if self.nearFieldScramblingGain is None:
      self.nearFieldScramblingGain = getScramblingGain(self.inObject, self.nfObject)
    return self.nearFieldScramblingGain

  def getFarFieldScramblingGain(self):
    if self.farFieldScramblingGain is None:
      self.farFieldScramblingGain = getScramblingGain(self.inObject, self.ffObject)
    return self.farFieldScramblingGain

  def getScramblingGain(self, inObject, outObject):
    inCentroidY, inCentroidX = inObject.getFiberCentroid()
    inCenterY, inCenterX = inObject.getFiberCenterEdgeMethod()
    inDiameter = inObject.getFiberDiameter()

    outCentroidY, outCentroidX = outObject.getFiberCentroid()
    outCenterY, outCenterX = outObject.getFiberCenterEdgeMethod()
    outDiameter = outObject.getFiberDiameter()

    deltaD_in = np.sqrt((inCentroidX - inCenterX)**2 + (inCentroidY - inCenterY)**2)
    deltaD_out = np.sqrt((outCentroidX - outCenterX)**2 + (outCentroidY - outCenterY)**2)

    scramblingGain = (deltaD_in / inDiameter) / (deltaD_out / outDiameter)

    return scramblingGain

  def getModalNoise(self):
    
