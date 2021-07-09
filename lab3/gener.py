import random
import string

def generator(size, chars=string.ascii_lowercase + ' '):
	return ''.join(random.choice(chars) for _ in range(size))

n = int(input())
print(n)
print(generator(n))

print(generator(n))