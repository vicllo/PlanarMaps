from sage.all import Permutation, Permutations
from MapError import InvalidMapPermutationArgument


class MapPermutation:
    def __init__(self, lst) -> None:
        # if isinstance(lst, Permutation):
        #     self._init_from_permutation(lst)
        #     return
        # try:
        #     if lst == int(lst) or lst > 0:
        #         self._init_from_number(lst)
        #         return
        # except:
        #     pass

        # try:
        #     if type(lst[0]) == type((42,)):
        #         self._init_from_cycle_list(lst)
        #         return
        #     self._init_from_list(lst)
        # except:
        #     raise InvalidMapPermutationArgument()

        if isinstance(lst, Permutation):
            self._init_from_permutation(lst)
        elif isinstance(lst, int):
            self._init_from_number(lst)
        elif isinstance(lst, list) and isinstance(lst[0], tuple):
            self._init_from_cycle_list(lst)
        elif isinstance(lst, list) and isinstance(lst[0], int):
            self._init_from_list(lst)
        else:
            raise InvalidMapPermutationArgument()

    def _init_from_cycle_list(self, lst):
        self._perm = Permutation(lst)

    def _init_from_number(self, n):
        self._perm = Permutations(n).identity()

    def _init_from_permutation(self, perm):
        self._perm = perm

    def _init_from_list(self, lst):
        self._perm = Permutation(lst)

    def size(self):
        return self._perm.size()

    def apply(self, i):
        if i > self.size():
            return i
        return self._perm(i)

    def __repr__(self) -> str:
        return str(list(self))

    def pretty_repr(self):
        return f"MapPermutation: {self.to_cycles()}"

    def pretty_print(self):
        """
        Print self in a more pretty form
        """
        print(self.pretty_repr())

    def to_cycles(self):
        return self._perm.to_cycles()

    def inverse(self):
        """
        This function calculate  the inverse of self
        -------
        Returns:
            - The inverse of self
        -------
        O(n)
        where n is the number of element of the permutation
        -------
        """
        return MapPermutation(self._perm.inverse())

    def __call__(self, i):
        return self.apply(i)

    def __getitem__(self, id):
        return self(id)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index < len(self):
            result = self(self.index + 1)
            self.index += 1
            return result
        else:
            raise StopIteration

    def __len__(self):
        return self.size()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return list(other) == list(self)
        return False

    def left_action_product(self, rperm):
        """
        This function calculate self*perm where * is the composition operation between permutation
        -------
        Args:
            -rperm: Another MapPermutation
        Returns:
            -A MapPermutation of size max(rperm.size(),self.size()) representing the composition self*rperm
        -------
        O(max(n,t))
        where n is the size of self and t the size of rperm
        -------
        """
        outSize = max(self.size(), rperm.size())

        outList = [self(rperm(i)) for i in range(1, outSize + 1)]

        return MapPermutation(outList)

    def number_of_fixed_points(self):
        """
        Returns: the number of fixed point ( we only consider i such that i<=self.size())
        """
        return self._perm.number_of_fixed_points()

    def right_action_product(self, lperm):
        """
        This function calculate lperm*self where * is the composition operation between permutation
        -------
        Args:
            -lperm: Another map permutation
        Returns:
            -A MapPermutation of size max(lperm.size(),self.size()) representing the composition lperm*self
        -------
        O(max(n,t))
        where n is the size of self and t the size of lperm
        -------
        """
        return lperm.left_action_product(self)

    def __mul__(self, b):
        return self.left_action_product(b)
