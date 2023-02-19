import concurrent.futures
import math

PRIMES = \
    [
    112272535095293,112582705942171,112272535095293,115280095190773,115797848077099,1099726899285419,
    109972700231,109972700261,109972700269,109972700323,109972700381,109972700389,109972700521,109972700573,109972700599,
    109972700603,109972700653,109972700683,109972700693,109972700743,109972700797,109972700801,109972700833,109972700839,
    109972700843,109972700863
    ]
PRIMES += PRIMES
PRIMES += PRIMES
PRIMES += PRIMES
PRIMES += PRIMES
PRIMES += PRIMES

def is_prime(n):
    if n % 2 == 0:
        return False

    sqrt_n = int(math.floor(math.sqrt(n)))
    for i in range(3, sqrt_n + 1, 2):
        if n % i == 0:
            return False
    return True

def main():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for number, prime in zip(PRIMES, executor.map(is_prime, PRIMES)):
            print('%d is prime: %s' % (number, prime))
            pass

if __name__ == '__main__':
    main()