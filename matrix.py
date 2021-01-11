from abc import ABCMeta
from abc import abstractmethod
from random import randint


class MGP(metaclass=ABCMeta):
    @abstractmethod
    def check_primitive(self):
        pass

    @abstractmethod
    def check_nonlin(self):
        pass

    def non(self):
        self.check_nonlin()

    def prim(self):
        self.check_primitive()

    def check(self):
        return self.prim(), self.non()



class Matrix:
    def __init__(self, rows):
        self.rows = rows
        
    @property
    def cols(self):
        return zip(*self.rows)

    def __repr__(self):
        rv = ''
        for row in self.rows:
            rv += ' '.join(map(str, row)) + '\n'
        return rv
    
    @staticmethod
    def scalar(row, col):
        rv = 0
        for row_elem, col_elem in zip(row, col):
            rv += row_elem * col_elem 
        return rv

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            raise Exception('not matrix at all')
        cols = list(other.cols)
        new_rows = []
        for row in self.rows:
            new_row = tuple(self.scalar(row, col) for col in cols)
            new_rows.append(new_row)
        rv = self.__class__(new_rows)
        # print(rv)
        return rv

    def __le__(self, other):
        for row1, row2 in zip(self.rows, other.rows):
            if not all(x <= y for x, y in zip(row1, row2)):
                return False
        return True


class SquareMatrix(Matrix):
    def __init__(self, rows):
        if not all(len(rows) == len(row) for row in rows):
            raise Exception('sorry, not a square matrix')
        self.size = len(rows)
        super().__init__(rows)

    @property
    def unity(self):
        def delta(x, y):
            if x == y:
                return 1
            return 0

        rows = []
        for i in range(self.size):
           row = tuple(delta(i, j) for j in range(self.size))
           rows.append(row)
        return self.__class__(rows)


    def __pow__(self, deg):
        rv = self.unity
        for _ in range(deg):
            rv *= self
        return rv


class BinMatrix(SquareMatrix):
    def __init__(self, rows):
        for row in rows:
            if not all(el in (0, 1) for el in row):
                raise Exception('sorry, not a binary matrix')
        super().__init__(rows)

    @staticmethod
    def scalar(row, col):
        return 0 if sum(x * y for x, y in zip(row, col)) == 0 else 1


class TernaryMatrix(SquareMatrix, MGP):
    def __init__(self, rows):
        for row in rows:
            if not all(el in (0, 1, 2) for el in row):
                print(row)
                raise Exception('sorry, not a ternary matrix')
        super().__init__(rows)

    @property
    def nonlinear_matrix(self):
        rows = []
        for i in range(self.size):
           row = tuple(2 for j in range(self.size))
           rows.append(row)
        return self.__class__(rows)

    @property
    def primitive_matrix(self):
        rows = []
        for i in range(self.size):
           row = tuple(1 for j in range(self.size))
           rows.append(row)
        return self.__class__(rows)

    @staticmethod
    def scalar(row, col):
        def mul(a, b):
            if a not in (0, 1, 2) or b not in (0, 1, 2):
                raise Exception('not ternary')
            if a == 0 or b == 0:
                return 0
            elif a == b == 1:
                return 1
            else:
                return 2

        return max(mul(x, y) for x, y in zip(row, col))

    def check_nonlin(self):
        print('NONLIN CHECK:')
        A = self
        for i in range(100):
            if self.nonlinear_matrix <= A:
                print('<2>-alpha epxonent:', i + 1)
                return i + 1
            A *= self
        print('not after 100')
        return None

    def check_primitive(self):
        print('PRIMITIVE CHECK:')
        A = self
        for i in range(100):
            if self.primitive_matrix <= A:
                print('epxonent:', i + 1)
                return i + 1
            A *= self
        print('not after 100')
        return None

    def test_local(self, power, n, m, limit=40):
        prev_mat = self ** (power  - 1)
        if prev_mat.rows[n][m] == 2:
            print('error: not min')
            return prev_mat
        for i in range(limit):
            new_mat = self ** (power + i)
            if new_mat.rows[n][m] != 2:
                print('error ', str(power + i))
                return new_mat
        print('ok')

    @property
    def last_col(self):
        return list(self.cols)[-1]

    @property
    def nl(self):
        for i, el in enumerate(reversed(self.last_col)):
            if el == 2:
                return self.size - 1 - i
        return None

    @property
    def delt(self):
        return self.size - self.nl

    @property
    def n(self):
        return [el for el in self.last_col if el == 2]
            
    @property
    def d(self):
        return [el for el in self.last_col if el != 0]

    @property
    def D(self):
        Ds = [[] for _ in range(self.m)]
        i = -1
        for el in self.last_col:
            if el != 0:
                i += 1
            Ds[i].append(el)
        return Ds

    @property
    def E(self):
        Es = [[] for _ in range(self.l)]
        i = -1
        for el in self.last_col:
            if el == 2:
                i += 1
            Es[i].append(el)
        return Es

    @property
    def l(self):
        return len(self.n) + 1

    @property
    def m(self):
        return len(self.d) + 1

    @property
    def delta(self):
        return self.size - self.nl

    def tau(self, u):
        for i, el in enumerate(reversed(self.last_col[:u + 1])):
            if el != 0:
                return u - i
        return 0

    def niu(self, u):
        for i, el in enumerate(reversed(self.last_col[:u + 1])):
            if el == 2:
                return u - i
        return self.nl

    @property
    def nu(self):
        for i, el in enumerate(reversed(self.last_col)):
            if el == 2:
                num = self.size - 1 - i
                if (self.size - num) % 2 == 1:
                    return num
        return None

    @property
    def du(self):
        for i, el in enumerate(reversed(self.last_col)):
            if el != 0:
                num = self.size - 1 - i
                if (self.size - num) % 2 == 1:
                    return num
        return None

    @property
    def dla(self):
        for i, el in enumerate(self.last_col):
            if el != 0:
                if i % 2 == 1:
                    return i
        return None

    @property
    def nla(self):
        for i, el in enumerate(self.last_col):
            if el == 2:
                if i % 2 == 1:
                    break
        else:
            return None
        for i, el in enumerate(self.last_col):
            if el == 2:
                if i % 2 == 0:
                    return i
        return None

    def u_n_1(self, v):
        # return 2 * self.size - v - max(self.du - 2, self.nu) - 2
        if self.check_primitive() is None:
            print('not prim')
            raise ValueError('non prim')
        n = self.size
        # return 2 * n - v - max(self.du - self.delt, self.nu) - 2
        if self.nu is None:
            return 2 * self.size - v - self.du
        return 2 * self.size - v - max(self.du - 2, self.nu) - 2

    def o_dla(self, u, v):
        if u >= self.dla:
            print('u ', u, ' is less then dla: ', self.dla)
            raise ValueError(u)
        if self.check_primitive() is None:
            print('not prim')
            raise ValueError('non prim')
        add = [self.tau(u) + self.du - 2]
        if self.nu is not None:
            add.append(self.tau(u) + self.nu)
        if self.niu(u) < self.nl:
            add.append(self.niu(u) + self.du)
        return self.size * 2 - v + u - 1 - max(add)

    def hi(self, u):
        if u < self.dla:
            print('condition dla not met ', self.dla)
            raise ValueError(u)
        els = self.last_col[:u + 1]
        a = None
        b = None
        for i, el in enumerate(reversed(els)):
            if el != 0:
                num = u - i
                if a is None and num % 2 == 0:
                    a = num
                if b is None and num % 2 == 1:
                    b = num
            if a is not None and b is not None:
                return abs(a - b)
        raise ValueError('no odd elem')

    def hi2(self, u):
#        if u < self.nla:
#            print('condition nla not met ', self.nla)
#            raise ValueError(u)
        a = self.tau(u)
        for i, el in enumerate(reversed(self.last_col[:u + 1])):
            if el != 0:
                num = u - i
                if num % 2 != a % 2 and el == 2:
                    return a - num
        return None

    def ksi(self, u):
        return self.tau(u) - self.niu(u)

    @property
    def delt2(self):
        ddu = self.size - self.du + 2
        if self.nu is None:
            return ddu
        return min(ddu, self.size - self.nu)

    def contains(self, u):
        n_to_u = [
            (u - i) % 2 for i, el in
            enumerate(self.last_col[self.niu(u):u + 1])
            if el != 0
        ]
        return 0 not in n_to_u or 1 not in n_to_u

    def d_i(self, u, v):
        base = self.size - v + u - 1 - self.tau(u)
        if self.contains(u):
            add = [
                self.delt2,
                self.hi(u) + 2,
                self.ksi(u) + self.size - self.du,
            ]
            hdl = self.hi2(u)
            if hdl is not None:
                add.append(hdl)
            print('does not contain')
            return base + min(add)
        print('contains')
        return base + min(
            max(self.tau(u) - self.niu(u), 2),
            self.delt2,
            self.hi(u) + 2,
        )


# def ShiftRegisterMatrix(*col):
def SRM(*col):
    size = len(col)

    def delta(x, y):
        return 1 if y == x - 1 else 0
    
    print('  '.join(map(str, col)))
    print(' '.join(map(lambda x: str(x[0]).ljust(2), enumerate(col))))
    rows = []
    for i, last in enumerate(col):
        row = tuple(delta(i, j) for j in range(size - 1)) + (last,)
        rows.append(row)
    matr = TernaryMatrix(rows)
    print(matr.size)
    matr.check()
    # print(matr)
    return matr

NUM = 24
            
#if __name__ == '__main__':
#    def delta(x, y, z):
#        if x == 0: return 1
#        elif x == z: return 1
#        elif x == y != 0: return 2
#        else: return 0
#
#    def undelta(x, y, z):
#        if x == 0: return 2
#        elif x == z: return 1
#        elif x == y != 0: return 1
#        else: return 0
#
#    for i in range(NUM - 5, NUM):
#    # for i in range(5):
#        ran = randint(1, NUM - 1)
#        ShiftRegisterMatrix(*tuple(delta(x, i, ran) for x in range(NUM)))
#        ShiftRegisterMatrix(*tuple(undelta(x, i, ran) for x in range(NUM)))


def h(m, u, v):
    m.test_local(m.d_i(u, v), u, v)

def hdi(m, v):
    for i in range(m.dla, m.nl + 1):
        h(m, i, v)
