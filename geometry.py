class Point:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __str__(self):
        return '({0}, {1}, {2})'.format(self.x, self.y, self.z)

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y, self.z+other.z)
 
    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    def __mul__(self, value):
        return Point(value*self.x, value*self.y, value*self.z)

    __rmul__ = __mul__

    # def __rmul__(self, value):
    #     return Point(value*self.x, value*self.y, value*self.z)
