# point3d.py

import math

class Point3D:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def distance_to(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return math.sqrt(dx**2 + dy**2 + dz**2)

    def __repr__(self):
        return f"Point3D({self.x}, {self.y}, {self.z})"
