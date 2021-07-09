import time
import random


A = 12345
B = 67890

def elliptic_curve(x, y, p):
    return (y ** 2) % p == (x ** 3 + A * x + B) % p


def print_curve(p):
    print("y^2 = x^3 + {0} * x + {1} (mod {2})".format(A, B, p))


def extended_euclidean_algorithm(a, b):
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = b, a

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    return old_r, old_s, old_t


def inverse_of(n, p):
    gcd, x, y = extended_euclidean_algorithm(n, p)
    assert (n * x + p * y) % p == gcd

    if gcd != 1:
        raise ValueError('inverse error')
    else:
        return x % p


def add_points(P, Q, p):
    if P == (0, 0):
        return Q
    elif Q == (0, 0):
        return P
    elif P[0] == Q[0] and P[1] != Q[1]:
        return (0, 0)
   
    if P == Q:
        m = ((3 * P[0] ** 2 + A) * inverse_of(2 * P[1], p)) % p
    else:
        m = ((P[1] - Q[1]) * inverse_of(P[0] - Q[0], p)) % p

    x = (m ** 2 - P[0] - Q[0]) % p
    y = (P[1] + m * (x - P[0])) % p
    return (x, -y % p)


def order_point(point, p):
    i = 1
    check = point
    while check != (0, 0):
        check = add_points(check, point, p)
        i += 1
    return i


if __name__ == '__main__':
    p = 27017
    print_curve(p)
    points = []
    start = time.time()
    for x in range(0, p):
        for y in range(0, p):
            if elliptic_curve(x, y, p):
                points.append((x, y))

    cnt_points = len(points)+1
    print("Group(curve) oder is {0}".format(cnt_points))
    point = random.choice(points)
    print("Order of P={0} is {1}".format(point, order_point(point, p)))
    print("Time: {} min.".format((time.time() - start)/60))