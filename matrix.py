from abc import ABCMeta
from abc import abstractmethod
from random import randint
from random import choice


global DEBUG
DEBUG = False

def printf(*s):
    if DEBUG:
        print(*s)


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
        printf('NONLIN CHECK:')
        A = self
        for i in range(100):
            if self.nonlinear_matrix <= A:
                printf('<2>-alpha epxonent:', i + 1)
                return i + 1
            A *= self
        printf('not after 100')
        return None

    def check_primitive(self):
        printf('PRIMITIVE CHECK:')
        A = self
        for i in range(100):
            if self.primitive_matrix <= A:
                printf('epxonent:', i + 1)
                return i + 1
            A *= self
        printf('not after 100')
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
    def reg(self):
        return list(self.cols)[-1]

    def state(self):
        print('  '.join(map(str, self.reg)))
        print(' '.join(map(
            lambda x: str(x[0]).ljust(2),
            enumerate(self.reg),
        )))

    @property
    def nl(self):
        for i, el in enumerate(reversed(self.reg)):
            if el == 2:
                return self.size - 1 - i
        return None

    @property
    def n0(self):
        for i, el in enumerate(self.reg):
            if el == 2:
                return i
        return None


    @property
    def delt(self):
        return self.size - self.nl

    @property
    def n(self):
        return [el for el in self.reg if el == 2]
            
    @property
    def d(self):
        return [el for el in self.reg if el != 0]

    @property
    def Ds(self):
        rv = [[] for _ in range(self.m)]
        i = -1
        for num, el in enumerate(self.reg):
            if el != 0:
                i += 1
            rv[i].append(num)
        return rv

    @property
    def Es(self):
        rv = [[] for _ in range(self.l)]
        i = -1
        for num, el in enumerate(self.reg):
            if el == 2:
                i += 1
            rv[i].append(num)
        return rv

    @property
    def E(self):
        return [i for i, el in enumerate(self.reg) if el == 2]

    @property
    def D(self):
        return [i for i, el in enumerate(self.reg) if el != 0]

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
        for i, el in enumerate(reversed(self.reg[:u + 1])):
            if el != 0:
                return u - i
        return 0

    def niu(self, u):
        for i, el in enumerate(reversed(self.reg[:u + 1])):
            if el == 2:
                return u - i
        return self.nl

    @property
    def nu(self):
        for i, el in enumerate(reversed(self.reg)):
            if el == 2:
                num = self.size - 1 - i
                if (self.size - num) % 2 == 1:
                    return num
        return None

    @property
    def du(self):
        for i, el in enumerate(reversed(self.reg)):
            if el != 0:
                num = self.size - 1 - i
                if (self.size - num) % 2 == 1:
                    return num
        return None

    @property
    def dla(self):
        for i, el in enumerate(self.reg):
            if el != 0:
                if i % 2 == 1:
                    return i
        return None

    @property
    def nla(self):
        for i, el in enumerate(self.reg):
            if el == 2:
                if i % 2 == 1:
                    break
        else:
            return None
        for i, el in enumerate(self.reg):
            if el == 2:
                if i % 2 == 0:
                    return i
        return None

    def u_n_1(self, v):
        if self.check_primitive() is None:
            print('not prim')
            raise ValueError('non prim')
        return self.size - v + self.delt2 - 2

    def o_dla(self, u, v):
        if u >= self.dla:
            print('u ', u, ' is less then dla: ', self.dla)
            raise ValueError(u)
        if self.check_primitive() is None:
            print('not prim')
            raise ValueError('non prim')
        base = self.size - v + u - self.tau(u) - 1
        add = [self.delt2]
        if u >= self.n0:
            add.append(self.ksi(u) + self.size - self.du)
        return base + min(add)

    def hi(self, u):
        if u < self.dla:
            print('condition dla not met ', self.dla)
            raise ValueError(u)
        els = self.reg[:u + 1]
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
        a = self.tau(u)
        for i, el in enumerate(reversed(self.reg[:u + 1])):
            if el != 0:
                num = u - i
                if num % 2 != a % 2 and el == 2:
                    return a - num

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
            enumerate(self.reg[self.niu(u):u + 1])
            if el != 0
        ]
        return 0 in n_to_u and 1 in n_to_u

    def local_exp_2(self, u, v):
        if u < self.dla:
            printf('u < dla')
            return self.o_dla(u, v)
        if u == self.size - 1:
            printf('u = n - 1')
            return self.u_n_1(v)
        return self.d_i(u, v)

    @property
    def exp_2(self):
        return max(self.local_exp_2(u, 0) for u in range(self.size))

    def eta(self, u):
        if self.tau(u) - self.hi(u) not in self.E:
            return 2
        if self.hi(u) == 1 and self.tau(u) not in self.E:
            return 1
        return 0

    def d_i(self, u, v):
        base = self.size - v + u - 1 - self.tau(u)
        add = [self.delt2, self.hi(u) + self.eta(u)]
        if not self.contains(u):
            if u >= self.n0:
                add.append(self.ksi(u) + self.size - self.du)
            printf('does not contain')
            return base + min(add)
        printf('contains')
        return base + min(add)


# def ShiftRegisterMatrix(*col):
def SRM(*col):
    size = len(col)

    def delta(x, y):
        return 1 if y == x - 1 else 0
    
    rows = []
    for i, last in enumerate(col):
        rows.append(tuple(delta(i, j) for j in range(size - 1)) + (last,))
    matr = TernaryMatrix(rows)
    matr.state()
    print(matr.size)
    matr.check()
    return matr

NUM = 24
            

def h(m, u, v):
    m.test_local(m.d_i(u, v), u, v)


def hdi(m, v):
    for i in range(m.dla, m.nl + 1):
        h(m, i, v)


def check(m):
    return m.exp_2 == m.check_nonlin()


def gen_reg_n_2(count=10, limit=30):
    for _ in range(count):
        n = randint(3, limit)
        reg = [randint(0, 2) for _ in range(n)]
        reg[0] = choice((1, 2))
        odd = choice([i for i in range(n - 2) if (n - i) % 2 == 1])
        if odd == n - 2:
            raise ValueError('n - 2')
        reg[odd] = choice((1, 2))
        reg[n - 1] = 0
        reg[n - 2] = 2
        print('  '.join(map(str, reg)))
        if check(SRM(*reg)):
            print('ok')
        else:
            print('error')
            mat = SRM(*reg)
            print('exp_2: ', mat.exp_2)
            print('nonlin check: ', mat.check_nonlin())
            import pdb; pdb.set_trace()
            pass
            return
