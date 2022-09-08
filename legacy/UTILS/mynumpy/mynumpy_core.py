import sys
import math

int = type(1)
float = type(0.5)
FLOAT=float
INT=int
int64=int
float64=FLOAT
float32=FLOAT
pi=3.141592653589793

def zeros (n, dtype=FLOAT):
    seq = []
    for i in range (n):
        seq.append (0)
    return array (seq, dtype)    

def concatenate (l):
    if (not isinstance (l, list)) and (not isinstance (l, tuple)):
        print "not list here:", l, type(l), isinstance (l, list)
        return NotImplemented
    
    mtype = INT
    seq = []
    for e in l:
        if not isinstance (e, array):
            try:
                e = array (e)
            except:
                raise TypeError
        if e.type == FLOAT:
            mtype = FLOAT
        for a in e.seq:
            seq.append(a)
    
    return array (seq, mtype)

def deg2rad (a):
    return a * pi / 180.

def sin (a):
    if isinstance (a, int) or isinstance (a, float):
        return math.sin (a)
    elif isinstance (a, array):
        return array([math.sin(x) for x in a.seq], dtype = FLOAT)
    else:
        raise TypeError

def cos (a):
    if isinstance (a, int) or isinstance (a, float):
        return math.cos (a)
    elif isinstance (a, array):
        return array([math.cos(x) for x in a.seq], dtype = FLOAT)
    else:
        raise TypeError

def sqrt (a):
    if isinstance (a, int) or isinstance (a, float):
        return math.sqrt (a)
    elif isinstance (a, array):
        return array([math.sqrt(x) for x in a.seq], dtype = FLOAT)
    else:
        raise TypeError

def arccos (a):
    if isinstance (a, int) or isinstance (a, float):
        return math.acos (a)
    elif isinstance (a, array):
        return array([math.acos(x) for x in a.seq], dtype = FLOAT)
    else:
        raise TypeError
        
def dot (a, b):
    if not (isinstance (a, array) and isinstance (b, array)):
        return NotImplemented

    if a.length != b.length:
        print >> sys.stderr, "Cannot take dot product of vectors of different length. Aborting now\n"
        sys.exit (2)
    
    if a.isnumbers and b.isnumbers:
        ret = 0.
        for i in range (a.length):
            ret = ret + a.seq[i] * b.seq[i]
        return ret
    elif a.isarray and b.isnumbers:
        ret = []
        for i in range (a.length):
            ret.append (dot (a.seq[i], b))
        return array (ret, dtype=FLOAT)
    elif a.isarray and b.isnumbers:
        return dot (b, a)
    else:
        return NotImplemented

def cross (a, b):
    if not (isinstance (a, array) and isinstance (b, array)):
        return NotImplemented

    if a.length != 3 or b.length != 3:
        print >> sys.stderr, "Cannot take cross product of vectors of length != 3. Aborting now\n"
        sys.exit (2)
    
    if a.type == FLOAT or b.type == FLOAT:
        type = FLOAT
    else:
        type = INT
        
    return array([a[1] * b[2] - a[2] * b[1], -a[0] * b[2] + a[2] * b[0], a[0] * b[1] - a[1] * b[0]], type)

def floor (a):
    if isinstance (a, int):
        return a
    elif isinstance (a, float):
        return float(int(a))
    elif not isinstance (a, array):
        raise TypeError
    else:
        pass

    if a.type == INT:
        return a.copy()

    return array([float(int(x)) for x in a.seq], FLOAT)
    
def rint (a):
    if isinstance (a, int):
        return a
    elif isinstance (a, float):
        return float(int(a + 0.5))
    elif not isinstance (a, array):
        raise TypeError
    else:
        pass
    
    if a.type == INT:
        return a.copy()

    return array([float(int(x + 0.5)) for x in a.seq], FLOAT)

class array (object):
    def __init__ (self, seq, dtype=FLOAT):
        if isinstance (seq, array):
            seq = seq.seq
            
        if not isinstance (seq, list) and not isinstance (seq, tuple):
            print >> sys.stderr, "Cannot initialize array from not list or tuple. Aborting", seq, type (seq)
            raise TypeError
            return None
        
        if isinstance (seq, tuple):
            seq = list (seq)
        
        # controlliamo la sanita' della sequenza
        self.isarray = False
        self.isnumbers = False
        self.seq = []
        for x in seq:
            if isinstance (x, array):
                self.seq.append (x)
                self.isarray = True
            elif isinstance (x, list) or isinstance (x, tuple):
                self.seq.append (array(x))
                self.isarray = True
            elif isinstance (x, FLOAT) or isinstance (x, INT):
                self.isnumbers = True
                if dtype == INT:
                    self.seq.append (int (x))
                elif dtype == FLOAT:
                    self.seq.append (float (x))
                else:
                    print >> sys.stderr, "dtype must be either INT or FLOAT"
                    raise TypeError
                    return None
            else:
                print >> sys.stderr, "seq can contain arrays or numbers, not frigging both", seq
                raise TypeError
                return None
        
        if self.isnumbers and self.isarray:
            print >> sys.stderr, "# mynumpy error: cannot mix arrays and numebers"
            raise TypeError
            return None
        
        self.type = dtype
        self.length = len (self.seq)
        self.index = 0
    
    #http://docs.python.org/tutorial/classes.html#iterators
    def __iter__ (self):
        return self
    
    def next (self):
        if self.index == self.length:
            self.index = 0
            raise StopIteration
        self.index += 1
        return self.seq[self.index - 1]

    def __getitem__ (self, index):
        if not isinstance (index, int):
            raise TypeError
        if index >= self.length or index < 0:
            raise IndexError
        return self.seq[index]
    
    def __setitem__ (self, index, value):
        if not isinstance (index, int):
            raise TypeError
        if index >= self.length or index < 0:
            raise IndexError
        if self.type == INT:
            self.seq[index] = int(value)
        elif self.type == FLOAT:
            self.seq[index] = float(value)
        else:
            return NotImplemented
    
    def copy (self):
        return array ([x for x in self.seq], dtype=self.type)
    
    def __len__ (self):
        return self.length
    
    def __str__ (self):
        if self.length == 0:
            return "[]"

        if self.isnumbers:
            ret = "["
            if self.type == INT:
                for s in self.seq[:-1]:
                    ret += "%i " % (s)
                ret += "%i]" % (self.seq[-1])
                return ret
            elif self.type == FLOAT:
                for s in self.seq[:-1]:
                    ret += "%g " % (s)
                ret += "%g]" % (self.seq[-1])
                return ret
            else:
                pass
        elif self.isarray:
            ret = "["
            for s in self.seq:
                ret += str(s)
            ret += ']'
            return ret
        else:
            print >> sys.stderr, "Shouldn't get here..."
            pass

    # emulating numeric types: http://docs.python.org/reference/datamodel.html#object.__imul__
    def __add__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        elif isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        
        if self.length == 0:
            return other.copy ()
        if other.length == 0:
            return self.copy ()

        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch"
            raise ValueError
        
        if self.type == FLOAT or other.type == FLOAT:
            type = FLOAT
        else:
            type = INT
        
        ret = []
        for i in range (self.length):
            ret.append (self.seq[i] + other.seq[i])
            
        return array(ret, type)
        
    def __sub__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        
        if self.length == 0:
            return -other.copy ()
        if other.length == 0:
            return self.copy ()

        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch"
            raise ValueError
        
        if self.type == FLOAT or other.type == FLOAT:
            type = FLOAT
        else:
            type = INT
        
        ret = []
        for i in range (self.length):
            ret.append (self.seq[i] - other.seq[i])
            
        return array(ret, type)
    
    def __mul__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass

        if self.length != other.length:
            print >> sys.stderr, "Cannot multiply vectors of different length. Aborting now\n"
            sys.exit (2)
        
        if self.type == FLOAT or other.type == FLOAT:
            type = FLOAT
        else:
            type = INT
        
        if self.isnumbers and other.isnumbers:
            # moltiplicazione membro a membro come fa numpy
            ret = []
            for i in range (self.length):
                ret.append (self.seq[i] * other.seq[i])
            return array(ret, type)
        elif self.isarray and other.isnumbers:
            # moltiplicazione matrice per vettore; numpy
            # fa come se il vettore fosse una matrice per colonne
            if len (self.seq) != other.length:
                return NotImplemented
            newself = self.transpose ()
            newseq  = []
            for i in range (other.length):
                newseq.append (newself.seq[i] * other[i])
            ret = array (newseq, dtype=type)    
            return ret.transpose() 
        elif self.isnumbers and other.isarray:
            return other.__mul__ (self)    
        else:
            return NotImplemented
    
    def __div__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        
        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch"
            raise ValueError
        
        if self.type == FLOAT or other.type == FLOAT:
            type = FLOAT
        else:
            type = INT
        
        ret = []
        for i in range (self.length):
            ret.append (self.seq[i] / other.seq[i])
        
        return array(ret, type)
        
    def __truediv__ (self, other):
        return self.__div__(other)
    
    def __radd__ (self, other):
        return self.__add__(other)
        
    def __rsub__ (self, other):
        other = array (other)
        return other.__sub__(self)
    
    def __rmul__ (self, other):
        if isinstance (other, float):
            other = array ([other,] * self.length, FLOAT)
        elif isinstance (other, int):
            other = array ([other,] * self.length, INT)
        else:
            other = array (other)
        
        return other.__mul__(self)
        
    def __rdiv__ (self, other):
        if isinstance (other, float):
            other = array ([other,] * self.length, FLOAT)
        elif isinstance (other, int):
            other = array ([other,] * self.length, INT)
        else:
            other = array (other)
        
        return other.__div__(self)
    
    def __iadd__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch", self, other
            raise ValueError
        
        for i in range (self.length):
            self.seq[i] += other.seq[i]
        
        if other.type == FLOAT:
            self.type == FLOAT
     
        return self    
    
    def __isub__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        
        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch"
            raise ValueError
        
        for i in range (self.length):
            self.seq[i] -= other.seq[i]
        
        if other.type == FLOAT:
            self.type == FLOAT
       
        return self 
    
    def __imul__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        
        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch"
            raise ValueError
        
        for i in range (self.length):
            self.seq[i] *= other.seq[i]
        
        if other.type == FLOAT:
            self.type == FLOAT
       
        return self 
    
    def __idiv__ (self, other):
        if isinstance (other, float):
            other = array ([other, ] * self.length, FLOAT)
        if isinstance (other, int):
            other = array ([other, ] * self.length, INT)
        elif isinstance (other, list) or isinstance (other, tuple):
            other = array (other, dtype = FLOAT)
        elif not isinstance (other, array):
            return NotImplemented
        else:
            pass
        
        if self.length != other.length:
            print >> sys.stderr, "Shape mismatch"
            raise ValueError
        
        try:
            for i in range (self.length):
                self.seq[i] /= other.seq[i]
        except ZeroDivisionError:
            print >> sys.stderr, "Warning: invalid value encountered in divide"
            raise ZeroDivisionError
            
        if other.type == FLOAT:
            self.type == FLOAT
       
        return self

    def __neg__ (self):
        if self.type == INT:
            return self * (-1)
        else:
            return self * (-1.)

    def transpose (self):
        # la implementiamo solo per matrici 3x3
        if self.length != 3:
            print >> sys.stderr, "Transpose not implemented for matices other than 3x3"
            raise ValueError
            return NotImplemented
        for s in self.seq:    
            if s.length != 3:
                print >> sys.stderr, "Transpose not implemented for matices other than 3x3"
                raise ValueError
        
        newseq = []
        newseq.append (array ([self.seq[0][0], self.seq[1][0], self.seq[2][0]], dtype=self.type)) 
        newseq.append (array ([self.seq[0][1], self.seq[1][1], self.seq[2][1]], dtype=self.type)) 
        newseq.append (array ([self.seq[0][2], self.seq[1][2], self.seq[2][2]], dtype=self.type)) 
        
        return array(newseq, self.type)
    
    T = property (transpose)        

class matrix (array):
    pass

        
import random as DefRandom
class Random (object):
    def __init__ (seed=None):
        DefRandom.seed (seed)
    
    def random_sample (self, n):
        seq = []
        for i in range (n):
            seq.append (DefRandom.random())
        return array(seq, FLOAT)

    def random (self):
        return DefRandom.random ()
    
    def rand (self, n):
        seq = []
        for i in range (n):
            seq.append (DefRandom.random ())
        return array (seq, FLOAT)
    
    def randint (self, low, high=0, size=1):
        try:
            low = int(low)
            high = int(high)
            size = int(size)
        except:
            raise ValueError
            return None
        if low >= high:
            print >> sys.stderr, "low >= high"
            raise ValueError
            return None
        seq = []
        for i in range (size):
            # randint and numpy.randint behave differently, thus the -1
            seq.append (DefRandom.randint (low, high - 1))

        return array(seq, INT)

# default instance of the random object
random = Random ()

class Linalg (object):
    def __init__ (self):
        self.ciao = True
    
    def det (self, a):
        if not isinstance (a, array):
            raise TypeError
        
        # la implementiamo solo per matrici 3x3
        if a.length != 3:
            print >> sys.stderr, "Determinant not implemented for matices other than 3x3"
            raise ValueError
            return NotImplemented
        for s in a.seq:    
            if s.length != 3:
                print >> sys.stderr, "Determinant not implemented for matices other than 3x3"
                raise ValueError
        res =  a.seq[0][0] * (a.seq[0][1] * a.seq[1][2] - a.seq[0][2] * a.seq[1][1])
        res +=-a.seq[1][0] * (a.seq[0][1] * a.seq[2][2] - a.seq[0][2] * a.seq[2][1])
        res += a.seq[2][0] * (a.seq[1][1] * a.seq[2][2] - a.seq[1][2] * a.seq[2][1])
        return res 
        
linalg = Linalg ()            
