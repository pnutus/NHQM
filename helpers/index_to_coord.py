from imports import *

def index_to_coord(i, size):
    if not 0 <= i < size**3:
    raise IndexError('Index is out of bounds')
    x = i % size
    tmp = i // size
    y = tmp % size
    z = tmp / size
    return sp.array([x, y, z]) - (size - 1)/2.

#
# Tests
#

import unittest
   
class IndexToCoordinateTests(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def testSizeOne(self):
        size = 1
        test = index_to_coord(0, size)
        correct = sp.array([0, 0, 0])
        self.assertTrue( sp.array_equal(test, correct) )
        
    def testSizeThree(self):
        size = 3
        test = index_to_coord(13, size)
        correct = sp.array([0, 0, 0])
        self.assertTrue( sp.array_equal(test, correct) )
        
    def testOverflow(self):
        with self.assertRaises(IndexError):
            size = 2
            test = index_to_coord(10, size)
        with self.assertRaises(IndexError):
            size = 5
            test = index_to_coord(-1, size)

if __name__ == '__main__':
    unittest.main()