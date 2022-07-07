import base, generators, readers, utils
import doctest

if __name__ == "__main__":
    doctest.testmod(base)
    doctest.testmod(generators)
    doctest.testmod(readers)
    doctest.testmod(utils)
