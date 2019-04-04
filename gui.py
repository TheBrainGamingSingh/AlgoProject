### THE QUADRATIC SIEVE FACTORING ALGORITHM

### MANPREET SINGH JUNEJA, DIVYANSHI GULATI & STUTI SHARMA
### PUNJAB ENGINEERING COLLEGE (DEEMED TO BE UNIV.)

### MANPREET SINGH JUNEJA.......17103022
### DIVYANSHI GULATI............17103029
### STUTI SHARMA................17103122
import math
import sys

sievingArraySize = 1000000

def SieveOfEratosthenes(B):
    prime = [True for i in range(B+1)]
    i = 2
    while (i * i <= B):
        if (prime[i] == True):
            for j in range(i * 2, B+1, i):
                prime[j] = False
        i += 1
    primes = [i for i in range(2, B) if prime[i]]
    return primes

def gcd(x, y):

    if x < 0:
        x *= -1

    if y < 0:
        y *= -1

    if x < y:
        x, y = y, x

    while x % y != 0:
        x = x % y
        if x < y:
            x, y = y, x
    return y

def eulerCriteria(N, primes):
    # return only values for which N is a quadratic residue
    return [prime for prime in primes if N ** ((prime - 1) // 2) % prime == 1]

def tonelliShanksAlgo(n, p):

    assert p % 2 == 1

    if n ** ((p - 1) // 2) % p != 1:
        return 0

    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    for z in range(2, 100):
        eulerCrit = z ** ((p - 1) // 2) % p
        if eulerCrit == p - 1:
            break
    else:

        return 0

    c = z ** Q % p
    R = n ** ((Q + 1) // 2) % p
    t = n ** Q % p
    M = S

    while t != 1:
        for i in range(1, M):
            if t ** (2 ** i) % p == 1:
                break
        else:
            print("Error")   #Remove this check

        b = c ** (2 ** (M - i - 1)) % p
        R = (R * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        M = i
    return R

def isSmooth(n, factorBase):
    if n == 0:
        return False

    for factor in factorBase:
        while n % factor == 0:
            n = n / factor
    return n == 1

def generateSmooth(n,factorBase, maxNum):
    ret = set({})
    startpoint = int(math.sqrt(n)) - sievingArraySize // 2

    endpoint = startpoint + sievingArraySize

    sieve = [x * x - n for x in range(startpoint, endpoint)]

    for factor in factorBase:
        # tonelli shanks algo doesn't work with factor 2
        if factor == 2:
            if n % 2 == 0:
                    residues = [0]
            else:
                residues = [1]
        else:
            R = tonelliShanksAlgo(n, factor)
            assert R != 0
            residues = [R, factor - R]



        for R in residues:
            startFrom = (startpoint - R + (factor - 1)) // factor
            endAt = (startFrom + (endpoint - (R + factor * startFrom) + (factor - 1)) // factor)

            for k in range(startFrom, endAt):
                x = (R + factor * k) - startpoint
                val = x + startpoint
                assert sieve[x] % factor == 0

                sieve[x] //= factor
                while(sieve[x] % factor == 0):
                    sieve[x] //= factor

                if sieve[x] == 1 and x != 0:
                    ret.add(val)
                    number = val * val - n
                    assert isSmooth(number, factorBase)

                    sieve[x] = 0
                    if len(ret) > maxNum:

                        return list(ret)

    return list(ret)


def generateVector(n, factorBase):
    ret = []
    for factor in factorBase:
        times = 0
        while n % factor == 0:
            times += 1
            n /= factor
        if times % 2 == 0:
            ret.append(0)
        else:
            ret.append(1)
    return ret


def getY(x, factorBase):
    #gets a modular square root
    y = 1
    for factor in factorBase:
        while x % (factor ** 2) == 0:
            x = x // (factor ** 2)
            y *= factor
    return y



def identityMatrix(height):
    return [[1 if i == j else 0 for j in range(height)] for i in range(height)]


def findLinearCombination(vectorList):
    height = len(vectorList)
    if height == 0:
        return None

    width = len(vectorList[0])
    if height < width:
        return None


    combinations = identityMatrix(height)

    for offset in range(width):
        if vectorList[offset][offset] == 0:

            for x in range(width):
                if vectorList[offset][x] != 0:
                    break
            else:

                return combinations[offset]

            for y in range(offset + 1, height):
                if vectorList[y][offset] != 0:
                    vectorList[y], vectorList[offset] = vectorList[offset], vectorList[y]
                    combinations[y], combinations[offset] = combinations[offset], combinations[y]
                    break
            else:
                continue

        for y in range(offset + 1, height):
            if vectorList[y][offset] == 0:
                continue
            for x in range(width):
                vectorList[y][x] *= -1
                vectorList[y][x] += vectorList[offset][x]
                vectorList[y][x] %= 2


            for x in range(height):
                combinations[y][x] *= -1
                combinations[y][x] += combinations[offset][x]
                combinations[y][x] %= 2

    return combinations[-1]

def generateDependentSubset(n,U, factorBase):

    vectorList = [generateVector(u * u - n, factorBase) for u in U]

    linearCombination = findLinearCombination(vectorList)

    if not linearCombination:
        return None

    return [u for num, u in enumerate(U) if linearCombination[num] == 1]



def factorize(N):
    #B = math.ceil(math.sqrt(math.exp(0.5 * math.log(N)* math.log(math.log(N)))))
    #print(N)
    B = 10000
    primes = SieveOfEratosthenes(B)

    factorBase = eulerCriteria(N, primes)
    #print("Factor Base used is : %s" % factorBase)


    U = generateSmooth(N,factorBase, len(factorBase) + 20)

    while True:
        print("...\n")
        dependentSet = generateDependentSubset(N,U, factorBase)
        print(dependentSet)
        if not dependentSet:
            print("No factors found - Algorithm doesn't work ")
            return 0, 0
        x = 1
        for u in dependentSet:
            x *= u

        yPrevious = 1
        for u in dependentSet:
            yPrevious *= u * u - N

        y = getY(yPrevious, factorBase)

        if x == y:
            U.remove(dependentSet[0])
            continue

        f1, f2 = gcd(x + y, N), gcd(x - y, N)

        if f1 != N and f1 != 1 and f2 != N and f2 != 1:
            return f1, f2
        U.remove(dependentSet[0])

if __name__ == "__main__":
    from tkinter import *
    root = Tk()
    root.configure(background="black")
    root.geometry("300x150")
    root.title("Quadratic Sieve Factorising Algorithm")
    topframe = Frame(root,bg = "black")
    topframe.pack()
    label1 = Label(topframe,bg="black",fg="white", text = "Quadratic Sieve Factorising Algorithm")
    label1.pack()
    label2 = Label(topframe,bg="red",fg="white", text = "Enter A Number")
    label2.pack()
    bottomframe = Frame(root,bg="black")
    bottomframe.pack(side = BOTTOM)

    number = StringVar()
    entry1 = Entry(root,width=40,textvariable = number).pack()

    factors = StringVar()
    def justDoIt():
        RSA = int(str(number.get()))
        print(RSA)
        factors.set("Code running...")
        if RSA % 2 == 0:
            factors.set("Even numbers not allowed")
        else:
            f1, f2 = factorize(RSA)
            factors.set("Factors are: %d and %d" % (f1, f2))


    button1 = Button(bottomframe,bg="green",fg="white",text = "Factorise",command = justDoIt)
    button1.pack()
    label4 = Label(bottomframe,bg="black",fg="white", text = "A project made by The Brain Gaming Singh")
    label4.pack()
    label5 = Label(bottomframe,bg="black",fg="white",text = "  ")
    label5.pack()

    #print("Factors(p and q): %d and %d" % (f1, f2))
    label3 = Label(root,bg="black",fg="white",textvariable = str(factors)).pack()

    root.mainloop()
