#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import sys
import math

# uses https://github.com/SolidCode/SolidPython
# Assumes SolidPython is in site-packages or elsewhwere in sys.path
from solid import *
from solid.utils import *
import numpy as np

SEGMENTS = 12

class Zome:
  def __init__(self, radius, n, strutLength, layers):
    self.radius = radius
    self.n = n # count of sides
    self.strutLength = strutLength

    self.vertices = [[0.0, 0.0, 0.0]]
    self.faces = []

    self.createTopFaces()
    for _ in range(layers):
      self.createNextLayer()    

  def createTopFaces(self):
    dPhi = 2 * math.pi / self.n
    # height of the top layer (apex to hieght of first circle of vertices)
    h = math.sqrt(self.strutLength * self.strutLength - self.radius * self.radius)
    for i in range(self.n):
      phi = i * dPhi
      self.vertices.append([
        self.radius * math.cos(phi),
        self.radius * math.sin(phi),
        -h])
    # Create the associated faces
    # Apex has vertix index 0
    # Face vertices order: Left bottom, top, right bottom
    offset = 1
    for i in range(self.n):
      nextI = (i + 1) % self.n
      self.faces.append([offset + i, 0, offset + nextI])

  def createNextLayer(self):
    # Faces made up of triangles poointing downwards
    for faceToMirror in self.faces[-self.n:]:
      npPoint0 = np.array(self.vertices[faceToMirror[0]])
      npPoint1 = np.array(self.vertices[faceToMirror[1]])
      npPoint2 = np.array(self.vertices[faceToMirror[2]])

      npNewVertex = npPoint0 - npPoint1 + npPoint2
      newVertexIndex = len(self.vertices)
      self.vertices.append(npNewVertex.tolist())

      # Vertex order: top left, top right, bottom
      self.faces.append([faceToMirror[0], faceToMirror[2], newVertexIndex])

    # Faces that fill the gaps
    offset = len(self.vertices) - self.n
    for i in range(self.n):
      nextI = (i + 1) % self.n
      # Vertex order: Left bottom, top, right bottom
      self.faces.append([offset + i, offset + nextI - self.n, offset + nextI])

  def createSolid(self):
    return polyhedron(self.vertices, self.faces)


def createEdge(startIndex, endIndex):
  return (startIndex, endIndex) if startIndex < endIndex else (endIndex, startIndex)

def createCylinders(vertices, faces, radius):
  # There should be only one strut per edge
  edgeKeysSet = set()

  objects = []
  for face in faces:
    for i in range(3):
      startVertexIndex = face[i]
      endVertexIndex = face[(i + 1) % 3]

      # Find the associated edge
      edgeKey = createEdge(startVertexIndex, endVertexIndex)
      # Has this edge already been handled?
      if edgeKey in edgeKeysSet:
        continue

      edgeKeysSet.add(edgeKey)

      startVertex = vertices[startVertexIndex]
      endVertex = vertices[endVertexIndex]

      objects.append(createCylinder(startVertex, endVertex, radius))
      #objects.push(cylinder({start: startPoint, end: endPoint, fn:fn, r: 0.05}));

  return union()(objects)

def createCylinder(startPoint, endPoint, radius):
  # Create transform that results in z axis pointing to endPoint - startPoint
  # See http://forum.openscad.org/Rods-between-3D-points-td13104.html

  npStartPoint = np.array(startPoint)
  npEndPoint = np.array(endPoint)
  # New z axis
  wDirection = npEndPoint - npStartPoint
  # normalize
  distance = np.linalg.norm(wDirection)

  #print(startPoint, endPoint, distance)
  w = wDirection / distance

  # The new x axis needs to be perpendicular to w.

  # Take a standard vector k which is not parallel to w, and form the cross
  # product w x k, which is guaranteed to be orthogonal to w.
  # The crux of the method is to take the vector k which is "the less parallel"
  # to w, i.e. that minimizes the dot product: take the index k corresponding to
  # the smallest |vk|, and you are on the safe side.
  u = findLeastParallel_(w)

  # The new y direction needs to be perpendicular to the new w (x) and u (z)
  # vectors.
  v = np.cross(w, u)

  matrix = np.transpose(np.vstack((u, v, w, npStartPoint)))
  matrix = np.vstack((matrix, [0,0,0,1]))

  strut = cylinder(r = radius, h = distance)

  return multmatrix(m = matrix.tolist())(strut)

def findLeastParallel_(w):
  standardVectors = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
  minDot = 1
  minIndex = 0
  for index, standardVector in enumerate(standardVectors):
    dot = abs(np.inner(w, standardVector))
    if dot < minDot:
      minDot = dot
      minIndex = index

  return standardVectors[minIndex]


if __name__ == '__main__':
  out_dir = sys.argv[1] if len(sys.argv) > 1 else os.curdir
  file_out = os.path.join(out_dir, 'zone.scad')

  zome = Zome(2.0, 21, 2.5, 11)
  
  print("strut length: ", zome.strutLength)

  showStruts = True
  o = None
  if showStruts:
    cylinderRadius = zome.strutLength * 0.1
    objects = createCylinders(zome.vertices, zome.faces, cylinderRadius)
    o = union()(objects)
  else:
    o = zome.createSolid()

  print("%(__file__)s: SCAD file written to: \n%(file_out)s" % vars())
  scad_render_to_file(o, file_out, file_header='$fn = %s;' % SEGMENTS, include_orig_code=True)
