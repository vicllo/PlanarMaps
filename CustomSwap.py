from MapPermutation import *


class CustomSwap (MapPermutation):
    """
    A custom class that permit to Initializes faster simple swap permutation
    """

    def __init__(self, lst) -> None:
        try:
            assert len(lst) == 1
            assert len(lst[0]) == 2

            self.a = min(lst[0][0], lst[0][1])
            self.b = max(lst[0][0], lst[0][1])

        except:
            raise

    def __eq__(self, other):
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False

    def size(self):
        return self.b

    def apply(self, i):
        if i != self.a and i != self.b:
            return i
        return self.a+self.b-i

    def number_of_fixed_points(self):
        return self.b-(1+self.a != self.b)

    def to_cycles(self):
        if self.a == self.b:
            return [(i,) for i in range(1, self.b+1)]

        return [(i,) for i in range(1, self.a)] + [(self.a, self.b)] + [(j,) for j in range(self.a, self.b)]

    def inverse(self):
        return self
