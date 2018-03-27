import numpy as np

class ObjectInScene():
    def __init__(self, name):
        self.name = name
        self.position = None
        self.orientation = None
        self.color = None
        self.transparency = None

    def setPosition(self, pos):
        self.position = pos

    def setOrientation(self, orientation):
        self.orientation = orientation

    def setColor(self, color):
        self.color = color

    def setTransparency(self, transparency):
        self.transparency = transparency

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

def reshape_array(array):
    temp = np.zeros(len(array))
    for i in range(temp.shape[0]):
        temp[i] = int(array[i])
    return temp

def load(filename):
    file = open('data.txt','r')
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
            else:
                print 'Unknown data is read.'
                break
    file.close()
    return models
