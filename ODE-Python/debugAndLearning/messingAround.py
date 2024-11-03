class Foo():
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def bar(self):
        print("bar", self.a)

    def zop(self):
        print("zip", self.b)

foo = Foo(1,2)
