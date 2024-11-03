class Foo():
    def __init__(self, a, b):
        self.parameters = {key: value for key, value in locals().items() if not key.startswith('__') and key != 'self'}
        self.a = a
        self.b = b

    def bar(self):
        print("bar", self.a)

    def zop(self):
        print("zip", self.b)

foo = Foo(1,2)
print(foo.parameters)