import gmpy2
import time
from gmpy2 import mpz

def gen_n(k, secure_k):
    prime = gmpy2.next_prime(2 ** (2 * k + secure_k + 2))
    while True:
        n_candidate = prime
        if n_candidate % 4 == 3:
            break
        prime = gmpy2.next_prime(prime)
    return n_candidate

def share(x, N):
    # P1
    random_state = gmpy2.random_state(int(time.time() * 10000))
    x1 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10001))
    x2 = gmpy2.mpz_random(random_state, N)
    x3 = (x - x1 - x2) % N
    return [x1, x2, x3]

def ranp(N):
    # P1
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a11 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10001))
    a12 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10002))
    a13 = gmpy2.mpz_random(random_state, N)
    # P2
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a21 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10001))
    a22 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10002))
    a23 = gmpy2.mpz_random(random_state, N)
    # P3
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a31 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10001))
    a32 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10002))
    a33 = gmpy2.mpz_random(random_state, N)
    # compute sharing
    a1 = (a11 + a21 + a31) % N
    a2 = (a12 + a22 + a32) % N
    a3 = (a13 + a23 + a33) % N
    return (a1, a2, a3)

def mult(share_tuple, N):
    share_sum = (share_tuple[0] + share_tuple[1] + share_tuple[2]) % N
    share_sum_square = (share_sum * share_sum) % N
    random_state = gmpy2.random_state(int(time.time() * 10000))
    share1 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10001))
    share2 = gmpy2.mpz_random(random_state, N)
    share0 = (share_sum_square - share1 - share2) % N
    return [share0, share1, share2]

def reveal(share_tuple, N):
    return (share_tuple[0] + share_tuple[1] + share_tuple[2]) % N

def square_root(x, N):
    square_root_x = gmpy2.powmod(x, gmpy2.t_div(N+1, 4), N)
    return square_root_x

def PRandBit(N):
    while True:
        ap_share = ranp(N)
        a2p_share = mult(ap_share, N)
        a2p = reveal(a2p_share, N)
        if a2p != 0:
            break
    b = square_root(a2p, N)
    c_share = 3 * [0]
    d_share = 3 * [0]
    for i in range(0, 3):
        c_share[i] = gmpy2.invert(b, N) * ap_share[i] % N
    share_1 = share(1, N)
    for i in range(0, 3):
        d_share[i] = (c_share[i] + share_1[i]) * gmpy2.invert(2, N) % N
    return d_share

def PRandFld(modulus):
    # P1
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a1 = gmpy2.mpz_random(random_state, int(modulus/3))
    # P2
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a2 = gmpy2.mpz_random(random_state, int(modulus/3))
    # P3
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a3 = gmpy2.mpz_random(random_state, int(modulus/3))
    return (a1, a2, a3)

def PRandInt(k, N):
    rand_int = PRandFld(gmpy2.powmod(2, k, N))
    return rand_int

def add_share(a, b, N):
    if not isinstance(a, list):
        a_share = share(a, N)
    else:
        a_share = a
    if not isinstance(b, list):
        b_share = share(b, N)
    else:
        b_share = b
    c_share = 3 * [0]
    for i in range(3):
        c_share[i] = (a_share[i] + b_share[i]) % N
    return c_share

def BitLTC(a, bi_list, k, secure_k, N):
    y_share = BitLTMap(a, bi_list, k, N)
    s_share = LSB(y_share, k, secure_k, N)
    return s_share

def dec_2_bin_share(a, k, N):
    a_bin_share = []
    a_bin_list = k * [0]
    a_bin = a.digits(2)
    for i in range(len(a_bin)):
        a_bin_list[i] = int(a_bin[i])
    for i in range(k):
        a_i = share(a_bin_list[i], N)
        a_bin_share.append(a_i)
    return a_bin_share

def reveal(share_tuple, N):
    return (share_tuple[0] + share_tuple[1] + share_tuple[2]) % N

def inner_prod(x_list, y_list, N):
    sum_share = [0, 0, 0]
    for i in range(len(x_list)):
        prod_share = mult_share(x_list[i], y_list[i], N)
        sum_share = add_share(sum_share, prod_share, N)
    return sum_share

def PRandFld(modulus):
    # P1
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a1 = gmpy2.mpz_random(random_state, int(modulus/3))
    # P2
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a2 = gmpy2.mpz_random(random_state, int(modulus/3))
    # P3
    random_state = gmpy2.random_state(int(time.time() * 10000))
    a3 = gmpy2.mpz_random(random_state, int(modulus/3))
    return [a1, a2, a3]

def mult_share(share_tuple_1, share_tuple_2, N):
    if not isinstance(share_tuple_1, list):
        share_1 = share_tuple_1
    else:
        share_1 = (share_tuple_1[0] + share_tuple_1[1] + share_tuple_1[2]) % N
    if not isinstance(share_tuple_2, list):
        share_2 = share_tuple_2
    else:
        share_2 = (share_tuple_2[0] + share_tuple_2[1] + share_tuple_2[2]) % N
    share_mult = (share_1 * share_2) % N
    random_state = gmpy2.random_state(int(time.time() * 10000))
    share1 = gmpy2.mpz_random(random_state, N)
    random_state = gmpy2.random_state(int(time.time() * 10001))
    share2 = gmpy2.mpz_random(random_state, N)
    share0 = (share_mult - share1 - share2) % N
    return [share0, share1, share2]

def PRandInv(modulus):
    while True:
        x = PRandFld(modulus)
        y = PRandFld(modulus)
        u = reveal(mult_share(x, y, modulus), modulus)
        if u != 0:
            break
    u_inv_y = mult_share(gmpy2.invert(u, modulus), y, modulus)
    return [x, u_inv_y]

def PreMulC(ai_list, k, N):
    ri_list = []
    si_list = []
    mi_list = []
    pj_list = []
    # step 1&2
    for i in range(k):
        ri_list.append(PRandInv(N))
    # step 3
    si_list.append(ri_list[0][0])
    # step 4&5
    for i in range(k-1):
        si = mult_share(ri_list[i][1], ri_list[i + 1][0], N)
        si_list.append(si)
    # step 6&7
    for i in range(k):
        mi = mult_share(si_list[i], ai_list[i], N)
        mi_list.append(reveal(mi, N))
    # print("mi_list =", mi_list)
    # step 8&9
    for j in range(k):
        m_prod = 1
        for i in range(j+1):
            m_prod = m_prod * mi_list[i] % N
        pj_list.append(mult_share(ri_list[j][1], m_prod, N))
    return pj_list

def BitLTMap(a, bi_list, k, N):
    d = [[0, 0, 0]]
    ai_list = k * [0]
    for i in range(len(a)):
        ai_list[i] = int(a[i])
    a_bin_share = dec_2_bin_share(a, k, N)
    # print(a_bin_share)
    y = []
    # step 1&2
    for i in range(len(bi_list) - 1):
        sum1 = add_share(ai_list[i + 1], bi_list[k - i - 2], N)
        sum2 = mult_share(-2*ai_list[i + 1], bi_list[k - i - 2], N)
        d.append(add_share(sum1, sum2, N))
    # print("d = ", d, len(d))
    # step 3
    dk_list = []
    for i in range(k-1):
        dk_list.append(add_share(d[k - i - 1], 1, N))
    # print("dk_list =", dk_list)
    x = PreMulC(dk_list, k-1, N)
    # print("x = ", x)
    # step 4&5
    for i in range(k):
        y.append(mult_share(bi_list[k - i - 1], 1-int(ai_list[i]), N))
    # print("y = ", y)
    # step 6
    y_share = add_share(y[k-1], inner_prod(list(reversed(x)), y[0:k-1], N), N)
    return y_share

def LSB(a_share, k, secure_k, N):
    # step 1
    u = PRandBit(N)
    # u = share(0, N)
    # step 2
    r = PRandInt(secure_k + k - 1, N)
    # r = share(3, N)
    # step 3
    c = 2 ** (k - 1) + reveal(a_share, N) + 2 * reveal(r, N) + reveal(u, N)
    # step 4
    c0 = int(c.digits(2)[-1])
    sum1 = add_share(c0, u, N)
    prod1 = mult_share(-2*c0, u, N)
    v = add_share(sum1, prod1, N)
    return v

def delete_LSB(x, m, N):
    x_trunc = (x * gmpy2.invert(gmpy2.powmod(2, m, N), N)) % N
    # print("invert 2^m = ", gmpy2.invert(gmpy2.powmod(2, m, N), N))
    return x_trunc

def bin_2_dec(bin_str, type):
    dec_value = 0
    bit_length = len(bin_str)
    if type == "int":
        for i in range(bit_length):
            dec_value += 2 ** (bit_length - i - 1) * int(bin_str[i])
    if type == "frac":
        for i in range(bit_length):
            dec_value += 2 ** -(i+1) * int(bin_str[i])
    return dec_value

def ff_2_fix(x, m, N):
    singed_x = x if x < int(N / 2) else N - x
    # print("singed_x = ", singed_x, mpz(singed_x).digits(2))
    bin_str = singed_x.digits(2)
    int_str = bin_str[0:-m]
    frac_str = bin_str[-m:]
    # print("int_str =", int_str)
    # print("frac_str =", frac_str)
    frac_value = bin_2_dec(frac_str, "frac")
    int_value = bin_2_dec(int_str, "int")
    fix_value = int_value + frac_value
    return fix_value if x < int(N / 2) else -fix_value

def dec_2_bin(x, m):
    ff_rac_value = 0
    xx = x
    for i in range(m):
        ff_rac = 2 ** (m - 1 - i) * int(xx * 2)
        ff_rac_value += ff_rac
        xx = xx * 2 - int(xx * 2)
    return ff_rac_value

def fix_2_ff(x, m):
    x_int = abs(int(x))
    ff_int = x_int * 2 ** m
    ff_frac = dec_2_bin(abs(x) - x_int, m)
    ff_value = ff_int + ff_frac
    return ff_value if int(x) > 0 else -ff_value

def reveal_bin(x, N):
    singed_x = x if x < int(N / 2) else N - x
    singed_x_bin = singed_x.digits(2)
    return singed_x_bin

def Trunc(a, k, m, secure_k, N):
    b = 3 * [0]
    ri_list = []
    r_p = 3 * [0]
    r = 3 * [0]
    c = 3 * [0]
    u = 3 * [0]
    a_p = 3 * [0]
    a_a_p = 3 * [0]
    d = 3 * [0]
    # step 1
    two_k_1_share = share(gmpy2.powmod(2, k - 1, N), N)
    for i in range(3):
        b[i] = a[i] + two_k_1_share[i]
    # print("b = ", b, reveal(b, N))
    # step 2 & 3
    for i in range(m):
        ri = PRandBit(N)
        ri_list.append(ri)
        # print("ri = ", ri, reveal(ri, N))
    # ri_list.append(share(1, N))
    # ri_list.append(share(0, N))
    # ri_list.append(share(1, N))
    # ri_list.append(share(0, N))
    # step 4
    for i in range(m):
        r_p = add_share(r_p, mult_share(ri_list[i], gmpy2.powmod(2, i, N), N), N)
    # step 5
    r_pp = PRandInt(secure_k + k - m, N)
    # r_pp = share(2, N)
    # print("r_pp = ", r_pp, reveal(r_pp, N))
    # step 6
    for i in range(3):
        r[i] = (gmpy2.powmod(2, m, N) * r_pp[i] + r_p[i]) % N
        # r[i] = r_p[i]
    # print("r = ", r, reveal(r, N))
    # step 7 & 8
    for i in range(3):
        c[i] = (b[i] + r[i]) % N
    # print("c = ", c, reveal(c, N))
    c_p = reveal(c, N) % gmpy2.powmod(2, m, N)
    # print(reveal(c, N), gmpy2.powmod(2, m, N))
    c_p_share = share(c_p, N)
    # BitLt
    u = BitLTC(c_p, list(reversed(ri_list)), m, secure_k, N)
    for i in range(3):
        a_p[i] = (c_p_share[i] - r_p[i] + gmpy2.powmod(2, m, N) * u[i]) % N
        # a_p[i] = (c_p_share[i] - r_p[i]) % N
    # step 10
    for i in range(3):
        a_a_p[i] = (a[i] - a_p[i]) % N
        d[i] = delete_LSB(a_a_p[i], m, N)

    # print("a =", reveal(a, N))
    # print("m =", m)
    # print("c_p =", reveal(c_p_share, N))
    # print("r_p =", reveal(r_p, N))
    # print("u = ", u, reveal(u, N))
    # print("a_p = ", a_p, reveal(a_p, N))
    # print("a_a_p = ", a_a_p, reveal(a_a_p, N))

    return d

def TruncPr(a, k, m, secure_k, N):
    b = 3 * [0]
    ri_list = []
    r_p = 3 * [0]
    r = 3 * [0]
    c = 3 * [0]
    a_p = 3 * [0]
    a_a_p = 3 * [0]
    d = 3 * [0]
    # step 1
    two_k_1_share = share(gmpy2.powmod(2, k - 1, N), N)
    for i in range(3):
        b[i] = a[i] + two_k_1_share[i]
    # step 2 & 3
    for i in range(m):
        ri = PRandBit(N)
        ri_list.append(ri)
    # step 4
    for i in range(m):
        r_p = add_share(r_p, mult_share(ri_list[i], gmpy2.powmod(2, i, N), N), N)
    # step 5
    r_pp = PRandInt(secure_k + k - m, N)
    # r_pp = share(2, N)
    # step 6
    for i in range(3):
        r[i] = (gmpy2.powmod(2, m, N) * r_pp[i] + r_p[i]) % N
        # r[i] = r_p[i]
    # step 7 & 8
    for i in range(3):
        c[i] = (b[i] + r[i]) % N
    # print("c = ", c, reveal(c, N))
    c_p = reveal(c, N) % gmpy2.powmod(2, m, N)
    c_p_share = share(c_p, N)
    # step 9
    for i in range(3):
        a_p[i] = (c_p_share[i] - r_p[i]) % N
    # step 10
    # a_a_p = share(-36, N)
    for i in range(3):
        a_a_p[i] = (a[i] - a_p[i]) % N
        d[i] = delete_LSB(a_a_p[i], m, N)

    # print("a = ", a, reveal(a, N))
    # print("b = ", b, reveal(b, N))
    # print("r_p = ", r_p, reveal(r_p, N))
    # print("r_pp = ", r_pp, reveal(r_pp, N))
    # print("r = ", r, reveal(r, N))
    # print("c = ", c, reveal(c, N))
    # print("c_p = ", c_p)
    # print("c_p_share = ", c_p_share, reveal(c_p_share, N))
    # print("a_p = ", a_p, reveal(a_p, N))
    # print("a_a_p = ", a_a_p, reveal(a_a_p, N), reveal(a_a_p, N).digits(2))

    return d

def LTZ(a, k, secure_k, N):
    s = Trunc(a, k, k-1, secure_k, N)
    signed_s = add_share(s, 1, N)
    return signed_s

def main():

    secure_k = 80
    k = 64
    m = 52
    N = gen_n(k, secure_k)
    print("N = ", N)

    # -2048 < x, y < 2048
    x = -12.5646535
    y = 45.55678178973
    z = x * y

    x_value = fix_2_ff(x, m)
    y_value = fix_2_ff(y, m)
    z_value = x_value * y_value % N

    print("x =", x)
    print("y =", y)
    print("real_z =", z)

    x_share = share(x_value, N)
    y_share = share(y_value, N)
    z_share = share(z_value, N)

    print("x_share =", x_share, ", x_value =", reveal(x_share, N))
    print("y_share =", y_share, ", y_value =", reveal(y_share, N))
    print("z_share =", z_share, ", z_value =", reveal(z_share, N))

    z_trunc_share = TruncPr(z_share, 2 * k, m, secure_k, N)
    z_trunc_value = reveal(z_trunc_share, N)
    z_trunc = ff_2_fix(z_trunc_value, m, N)

    print("z_trunc_share =", z_trunc_share, ", z_trunc_value =", reveal(z_trunc_share, N))
    print("z_share_bin =      ", reveal_bin(z_value, N))
    print("z_trunc_share_bin =", reveal_bin(z_trunc_value, N))
    print("z_trunc =", z_trunc)

    signed_bit_share = LTZ(z_trunc_share, k, secure_k, N)
    signed_bit = int(reveal(signed_bit_share, N))
    print("signed_bit_share =", signed_bit_share, signed_bit)
    if signed_bit == 0:
        print("z_trunc < 0")
    elif signed_bit == 1:
        print("z_trunc >= 0")
    else:
        raise Exception("Invalid signed_bit!", signed_bit)


    # for i in range(20):
    #     signed_bit_share = LTZ(z_trunc_share, k, secure_k, N)
    #     signed_bit = int(reveal(signed_bit_share, N))
    #     print("signed_bit_share =", signed_bit_share, signed_bit)
    #     if signed_bit == 0:
    #         print("z_trunc < 0")
    #     elif signed_bit == 1:
    #         print("z_trunc >= 0")
    #     else:
    #         raise Exception("Invalid signed_bit!", signed_bit)


if __name__ == "__main__":
    main()
