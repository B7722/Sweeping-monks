import numpy as np

class ObjectInScene():
    def __init__(self, name):
        self.name = name
        self.position = None
        self.orientation = None
        self.color = None
        self.transparency = None
        self.point1 = None
        self.point2 = None
        self.point3 = None
        self.point4 = None
        self.x = None
        self.y = None
        self.z = None
        self.r = None


    def setPosition(self, pos):
        self.position = pos

    def setOrientation(self, orientation):
        self.orientation = orientation

    def setColor(self, color):
        self.color = color

    def setTransparency(self, transparency):
        self.transparency = transparency

    def setPoint1(self, point1):
        self.point1 = point1

    def setPoint2(self, point2):
        self.point2 = point2

    def setPoint3(self, point3):
        self.point3 = point3

    def setPoint4(self, point4):
        self.point4 = point4

    def setX(self, x):
        self.x = x

    def setY(self, y):
        self.y = y

    def setZ(self, z):
        self.z = z

    def setR(self, r):
        self.r = r

    def getName(self):
        return self.name

    def getPosition(self):
        return self.position

    def getOrientation(self):
        return self.orientation

    def getColor(self):
        return self.color

    def getTransparency(self):
        return self.transparency

    def getPoint1(self):
        return self.point1

    def getPoint2(self):
        return self.point2

    def getPoint3(self):
        return self.point3

    def getPoint4(self):
        return self.point4

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getZ(self):
        return self.z

    def getR(self):
        return self.r

def reshape_array(array):
    temp = np.zeros(len(array))
    for i in range(temp.shape[0]):
        temp[i] = float(array[i])
    return temp

def load(filename):
    file = open('test.txt','r')
    models = []
    while True:
        line = file.readline()
        if line == '':
            break
        array = line.split()
        if array[0] == 'Name:':
            models.append(ObjectInScene(array[1]))
        else:
            field = array[0]
            array = reshape_array(array[1:])
            model = models[len(models)-1]
            if field == 'Position:':
                model.setPosition(array)
            elif field == 'Orientation:':
                model.setOrientation(array)
            elif field == 'Color:':
                model.setColor(array)
            elif field == 'Transparency:':
                model.setTransparency(array)
            elif field == 'Point1:':
                model.setPoint1(array)
            elif field == 'Point2:':
                model.setPoint2(array)
            elif field == 'Point3:':
                model.setPoint3(array)
            elif field == 'Point4:':
                model.setPoint4(array)
            elif field == 'X:':
                model.setX(array)
            elif field == 'Y:':
                model.setY(array)
            elif field == 'Z:':
                model.setZ(array)
            elif field == 'Radius:':
                model.setR(array)
            else:
                print('Unknown data is read.')
                break
    file.close()
    return models