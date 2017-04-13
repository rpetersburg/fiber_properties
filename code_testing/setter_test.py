class A(object):
    def __init__(self, a=0):
        self._a = a

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a += value

obj = A()
print obj.a

obj = A(3)
obj.a = 4
print obj.a