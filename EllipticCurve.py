# coding: utf-8

class FiniteField:
    """ 有限体 """
    def __init__(self, p):
        self.p = p

    def xgcd(self, a, b):
        x0, x1, y0, y1 = 1, 0, 0, 1
        while b != 0:
            q, a, b = a // b, b, a % b
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        return  a, x0, y0

    def modInv(self, a):
        g, x, y = self.xgcd(a, self.p)
        if g != 1:
            raise Exception("逆数が存在しません")
        return x % self.p

    def legendreSymbol(self, a, p):
        ls = pow(a, (p-1) / 2, p)
        return -1 if ls == p-1 else ls

    def modSqrt(self, a):
        L = self.legendreSymbol(a, self.p)
        if L == -1:
            raise Exception("根が存在しません")
        elif L == 0:
            return 0
        elif self.p == 2:
            return a # a=1->n=1
        # p-1 = s*2^eなるsとeを計算する
        s = (self.p - 1)>>1
        e = 0
        while s & 1 == 0:
            s >>= 1
            e += 1
        # L(n/p) = -1なるnを見つける
        n = 2
        while self.legendreSymbol(n, self.p) != -1:
            n += 1
        # 根を計算する
        x = pow(a, (s + 1)/2, self.p)
        b = pow(a, s, self.p)
        g = pow(n, s, self.p)
        while True:
            t, m = b, 0
            for m in xrange(e):
                if t == 1:
                    break
                t = pow(t, 2, self.p)
            if m == 0:
                return x
            gs = pow(g, 2 ** (e - m - 1), self.p)
            g = (gs ** 2) % self.p
            x = (x * gs) % self.p
            b = (b * g) % self.p
            e = m
            
class Point:
    """ 点 """
    def __init__(self, EC, x, y):
        self.EC = EC
        self.F = self.EC.F
        self.x = x
        self.y = y
        
    def __add__(self, Q):
        """ 加算の定義 """
        if not isinstance(Q, Point):
            raise("無効な型との加算です")
        # 異なる楕円曲線
        if Q.EC.F.p != self.EC.F.p:
            raise("異なる楕円曲線上の点どうしの加算です")
        # 無限遠点の場合
        if Q.x == -1:
            return self
        if self.x == -1:
            return Q
        # 無限遠点以外
        if self.x == Q.x:
            if self.y == Q.y:
                # R = P + P
                m = (3*Q.x*Q.x + self.EC.A) * self.F.modInv(2*Q.y) % self.F.p
            else:
                # R = P + (-P) = O
                return self.EC.O
        else:
            # R = P + Q
            m = (Q.y - self.y) * self.F.modInv(Q.x - self.x) % self.F.p
        x = (m*m - self.x - Q.x) % self.F.p
        y = (m*(self.x - x) - self.y) % self.F.p
        return Point(self.EC, x, y)

    def __mul__(self, n):
        """ スカラー倍の定義 """
        if not isinstance(n, int):
            raise("無効な型との乗算です")
        Q = self.EC.O
        R = self
        while n > 0:
            if n & 1 == 1:
                Q += R
            R += R
            n >>= 1
        return Q

    def __str__(self):
        return "({0}, {1})".format(self.x, self.y)

class EllipticCurve:
    """ 楕円曲線 """
    def __init__(self, F, tup):
        self.F = F
        self.A, self.B = tup
        self.D = (4 * pow(self.A, 3) + 27 * pow(self.B, 2)) % self.F.p
        if self.D == 0:
            raise("判別式が0になります")
        self.O = Point(self, -1, -1)
    
    def f(self, x):
        """ xからyを求める """
        y2 = (pow(x, 3) + self.A*x + self.B) % self.F.p
        y = self.F.modSqrt(y2)
        return y

    def __str__(self):
        return "EC: y^2 = x^2 + {A}x + {B} mod {p}".format(
            A=self.A, B=self.B, p=self.F.p
        )

if __name__ == '__main__':
    p = 17
    # 楕円曲線
    F = FiniteField(p)
    A, B = 1, 2
    EC = EllipticCurve(F, (A, B))
    # 例
    O = EC.O
    P = Point(EC, 0, EC.f(0)) # P(0, 4)
    Q = Point(EC, 1, EC.f(1)) # Q(1, 16)
    print(EC)
    print("P = " + str(P))
    print("Q = " + str(Q))
    print("P + Q = " + str(P + Q))
    print("P + P = " + str(P + P))
    print("2P = " + str(P*2))
    print("10P = " + str(P*10))
