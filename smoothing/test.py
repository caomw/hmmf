
class myclass:
    def __init__(self, a):
        self.data = a


a = myclass(3)
b = myclass(4)
c = a
c.data = 5

