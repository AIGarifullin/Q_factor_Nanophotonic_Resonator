from __future__ import division

import argparse
import math as mt
import numpy as np

import meep as mp
from meep import mpb

from geom import check_nonnegative, GeometricObject, Vector3



class EllipticalCylinder(GeometricObject):
    """
    A cylinder, with elliptical cross-section and finite height.

    **Properties:**

    + **`major_radius` [`number`]** — Major radius of the ellipse's cross-section. No default value.
    
    + **`minor_radius` [`number`]** — Minor radius of the ellipse's cross-section. No default value.

    + **`height` [`number`]** — Length of the cylinder along its axis. No default value.

    + **`axis` [`Vector3`]** — Direction of the cylinder's axis; the length of this vector
      is ignored. Defaults to `Vector3(x=0, y=0, z=1)`.
    """

    def __init__(self, major_radius, minor_radius, axis=Vector3(0, 0, 1), height=1e20, **kwargs):
        """
        Constructs a `EllipticalCylinder`.
        """
        self.axis = Vector3(*axis)
        self.major_radius = float(major_radius)
        self.minor_radius = float(minor_radius)
        self.height = float(height)
        super().__init__(**kwargs)

    @property
    def major_radius(self):
        return self._major_radius
    
    @property
    def minor_radius(self):
        return self._minor_radius

    @property
    def height(self):
        return self._height

    @major_radius.setter
    def major_radius(self, val):
        self._major_radius = check_nonnegative("EllipticalCylinder.major_radius", val)

    @minor_radius.setter
    def minor_radius(self, val):
        self._minor_radius = check_nonnegative("EllipticalCylinder.minor_radius", val)
        
    @height.setter
    def height(self, val):
        self._height = check_nonnegative("EllipticalCylinder.height", val)