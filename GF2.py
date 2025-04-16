# -*- coding: utf-8 -*-
#       __
#      /  \____
#    (  )      \
#    /  \      |  |
#    \__/      |__|_
#  Have fun!      >^.^<
"""
@author: byt3d1$t1ll3r

This script allows to build a finite field with 2^m elements (for m<25) from
an irreductible polynomial.

If this script is executed, some basic tests are displayed.

Important!:
This script was written only for didactic purposes, it is not optimal neither 
secure so it should not be used in production environments or for serious 
research.
"""
import numpy as np

class GF2():
    
    def __init__(self, m = None):
        """
        Class constructor.
        """
        self.__m = int(m)
        self.__q = int(2**m)
        auxp = np.zeros(24)
        auxp[0] = 0x3
        auxp[1] = 0x7
        auxp[2] = 0xb
        auxp[3] = 0x13
        auxp[4] = 0x25
        auxp[5] = 0x43
        auxp[6] = 0x89
        auxp[7] = 0x11b
        auxp[8] = 0x211
        auxp[9] = 0x409
        auxp[10] = 0x805
        auxp[11] = 0x1009
        auxp[12] = 0x201b
        auxp[13] = 0x4443
        auxp[14] = 0x8003
        auxp[15] = 0x1002b
        auxp[16] = 0x20009
        auxp[17] = 0x40081
        auxp[18] = 0x80027
        auxp[19] = 0x100009
        auxp[20] = 0x200005
        auxp[21] = 0x400003
        auxp[22] = 0x800021
        auxp[23] = 0x1000087
        self.__irrp =  int(auxp[m-1])
    
    def __del__(self):
        """
        Class destructor.        
        """
        del(self.__q)
        del(self.__irrp)
        
    def add(self,a,b):
        """
        This function executes the addition defined over GF(2^m) between a and b.
            
        @param a: finite field element [integer].
        @param b: finite field element [integer].
        
        @return a: a xor b.
        """
        return a^b
    
    def mult(self,a,b):
        """
        This function executes the multiplication defined over GF(2^m) between a and b.
            
        @param a: finite field element [integer].
        @param b: finite field element [integer].
        
        @return a: a . b mod (0x3).
        """
        i=0
        r=0
        while True:
            if (int(b)>>i)&1:
                r=r^(int(a)<<i)
            if (int(b)>>(i+1))==0:
                break
            i+=1
        a=r
        ncp=np.ceil(np.log2(self.__irrp+1))
        while True:
            nca=np.ceil(np.log2(a+1))
            if nca < ncp:
                break
            a = a ^ (self.__irrp<<int(nca-ncp))
        return a
    
    def pot(self,a,n):
        """
        This function executes the multiplication defined over GF(2^m) between a with itself n times.
            
        @param a: finite field element [integer].
        @param n: exponent value [integer].
        
        @return a: a^n.
        """
        if n == 0:
            return 0x1
        b = a
        for i in range(0,n-1):
            b = self.mult(b,a)
        return b

    def inv(self,a):
        """
        This function compute the multiplicative inverse over GF(2^m) of a using the extended Euclidean algorithm.
            
        @param a: finite field element [integer].
        
        @return a: a^{-1}.
        """
        r0, r1 = a, self.__irrp    
        s0, s1 = 1, 0    
        while True:
            d, b = r0, r1
            ncb=np.ceil(np.log2(b+1))
            c=0
            while True:
                ncd=np.ceil(np.log2(d+1))
                if ncd < ncb:
                    break
                c = c^1<<int(ncd-ncb)
                d = d ^ (b<<int(ncd-ncb))
            if not d:
                break
            r1, r0 = d, r1
            i=0
            r=0
            while True:
                if (int(c)>>i)&1:
                    r=r^(int(s1)<<i)
                if (int(c)>>(i+1))==0:
                    break
                i+=1
            s = s0^r
            s1, s0 = s, s1
        return s1

    
    def m(self):
        return self.__m


    def q(self):
        return self.__q

    def irrp(self):
        return self.__irrp

if __name__ == "__main__":

    import random as rnd
    ########################################################
    #Example: Operations over GF(2^16)                                               
    ########################################################   
    m=16
    a=int(rnd.randint(1,2**m))
    b=6
    print("a = " + str(a))
    print("b = " + str(b))
    field1 = GF2(m)
    s = field1.add(a,b)
    print("addition of a and b = " + str(s))
    m = field1.mult(a,b)    
    print("multiplication of a and b = " + str(m))
    i = field1.inv(a)   
    print("inverse of a = " + str(i))
