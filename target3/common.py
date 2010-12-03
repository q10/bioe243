from math import sqrt, log, exp

IMAX = 2.0**31.0

# This code follows exactly the algorithm described in Numerical Recipes
class Ran3:
    m = IMAX-1
    mseed = 161803398

    def __init__(self, init):
        mj = abs(Ran3.mseed-abs(init)) % Ran3.m
        self.rand_table = [0]*55
        self.rand_table[54] = mj
        self.iseed, ii, mk = 0, 0, 1

        # initialize rand_table using init as seed and large number mseed
        for i in range(0, 54):
            ii = (21*i) % 55
            self.rand_table[ii] = mk
            mk = mj - mk
            if mk < 0:
                mk += Ran3.m
            mj = self.rand_table[ii]
        
        # randomize table warming up the generator
        for k in range(0, 4):
            for i in range(0, 55):
                self.rand_table[i] -= self.rand_table[1+((i+30)%54)]
                if (self.rand_table[i] < 0):
                    self.rand_table[i] += Ran3.m

        

    # generate the next random number
    def generate(self):
        self.rand_table[self.iseed] = (self.rand_table[self.iseed] - self.rand_table[(self.iseed+31) % 55])
        if self.rand_table[self.iseed] < 0:
            self.rand_table[self.iseed] += Ran3.m
        self.iseed = (self.iseed+1) % 55

        return self.rand_table[self.iseed] / Ran3.m


class BoxMueller:
    def __init__(self, seed):
        self.ran3 = Ran3(seed)

    def generate(self):
        r = 1.0
        while r >= 1.0:
            x1 = (2.0*self.ran3.generate()) - 1.0
            x2 = (2.0*self.ran3.generate()) - 1.0
            r = (x1**2) + (x2**2)

        r = sqrt(-2.0*log(r)/r)
        return (x1*r, x2*r)
