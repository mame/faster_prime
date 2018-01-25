module FasterPrime
  SMALL_PRIMES = []
  SMALL_PRIME_TABLE = []

  def self.setup_table(n)
    n = (n + 59) / 60 * 60
    SMALL_PRIMES.clear << 2
    n /= 2
    SMALL_PRIME_TABLE.replace([true] * n)
    i = 1
    while i < n
      if SMALL_PRIME_TABLE[i]
        SMALL_PRIMES << j = i * 2 + 1
        k = i + j
        while k < n
          SMALL_PRIME_TABLE[k] = false
          k += j
        end
      end
      i += 1
    end
  end
  setup_table(100000)

  module Utils
    module_function

    FLOAT_BIGNUM = Float::RADIX ** (Float::MANT_DIG - 1)
    # Find the largest integer m such that m <= sqrt(n)
    def integer_square_root(n)
      if n < FLOAT_BIGNUM
        Math.sqrt(n).floor
      else
        # newton method
        a, b = n, 1
        a, b = a / 2, b * 2 until a <= b
        a = b + 1
        a, b = b, (b + n / b) / 2 until a <= b
        a
      end
    end

    # Find b such that a * b = 1 (mod n)
    def mod_inv(a, n)
      raise ArgumentError, "not coprime for mod_inv" if a == 0 || a.gcd(n) != 1

      # Extended GCD Algorithm:
      # Find $$(x, y)$$ such that $$ax + by = gcd(a, b), x > 0 and y > 0$$
      x0, x1 = 1, 0
      y0, y1 = 0, 1
      b = n
      while b != 0
        q = a / b
        a, b = b, a % b
        x0, x1 = x1, x0 - x1 * q
        y0, y1 = y1, y0 - y1 * q
      end

      x0 % n
    end

    # Calculate a^b mod n
    def mod_pow(a, b, n)
      r = 1
      return 1 if b == 0
      len = b.bit_length
      len.times do |i|
        if b[i] == 1
          r = (a * r) % n
          return r if i == len - 1
        end
        a = (a * a) % n
      end

      raise "assert not reached"
    end

    # Find x such that x^2 = a (mod prime)
    # See ``Complexity of finding square roots'' section of wikipedia:
    # http://en.wikipedia.org/wiki/Quadratic_residue
    def mod_sqrt(a, prime)
      return 0 if a == 0

      case
      when prime == 2
        a.odd? ? 1 : 0

      when prime % 4 == 3
        mod_pow(a, (prime + 1) / 4, prime)

      when prime % 8 == 5
        # Let v = (2a)^((prime-5)/8) and r = (2av^2 - 1)av.
        #
        # 2^((prime-1)/2) = -1 because (2|prime) = -1.
        # a^((prime-1)/2) =  1 because (a|prime) =  1.
        # So, (2a)^((prime-1)/2) = -1.
        #
        # r = ((2a)^((prime-1)/4) - 1) * a * (2a)^((prime-5)/8)
        # r^2 = ((2a)^((prime-1)/2) - 2(2a)^((prime-1)/4) + 1)
        #        * a^2 * (2a)^((prime-5)/4)
        #     = -2(2a)^((prime-1)/4) * a^2 * (2a)^((prime-5)/4)
        #     = -a * (2a)^((prime-1)/2)
        #     = a
        #
        # So r is the solution of x^2 = a (mod prime).
        v = mod_pow(a * 2, (prime - 5) / 8, prime)
        (2 * a * v * v - 1) * a * v % prime

      else # prime % 8 == 1
        # Use Tonelli–Shanks algorithm, see:
        # http://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm

        # Calculate (q, s) such that prime - 1 = 2^s * q and q is odd
        s, q = 3, prime - 1
        s += 1 until q[s] == 1
        q >>= s

        # Find quadratic nonresidue z.
        z = (2...prime).find {|z_| kronecker(z_, prime) != 1 }

        # The main part of the algorithm
        c = mod_pow(z, q, prime)
        x = mod_pow(a, (q - 1) / 2, prime)
        r = a * x % prime  # a^((q+1)/2)
        t = r * x % prime  # a^q
        until t == 1
          i, d = 0, t
          i, d = i + 1, d * d % prime while d != 1
          b = mod_pow(c, 1 << (s - i - 1), prime)
          c = b * b % prime
          r = r * b % prime
          t = t * c % prime
          s = i
        end
        r
      end
    end

    # p-adic valuation: Find the highest exponent v such that p^v divides n
    def valuation(n, p)
      e = 0
      while n % p == 0
        n /= p
        e += 1
      end
      e
    end

    # Enumerate divisors of n, including n itself
    def each_divisor(n)
      return enum_for(:each_divisor, n) unless block_given?
      t = 1
      while t * t < n
        yield t if n % t == 0
        t += 1
      end
      yield t if n == t * t
      while t > 1
        t -= 1
        yield n / t if n % t == 0
      end
    end

    # Calculate the number of positive integers less than or equal to n that
    # are coprime to n
    def totient(n)
      n = n.abs
      raise "n must not be zero" if n == 0
      r = n
      PrimeFactorization.prime_factorization(n) do |prime, exp|
        r = r * (prime - 1) / prime
      end
      r
    end

    # Finds the least primitive root modulo n.
    # Algorithm 1.4.4 (Primitive Root).
    def primitive_root(n)
      t = Utils.totient(n)
      ps = PrimeFactorization.prime_factorization(t).to_a
      2.upto(t) do |a|
        found = true
        ps.each do |q,|
          if mod_pow(a, t / q, n) == 1
            found = false
            break
          end
        end
        return a if found
      end
      nil
    end

    # Calculate Kronecker symbol (a|b) in O((log a)(log b)), see pp. 29--30,
    # Cohen, Henri (1993), A Course in Computational Algebraic Number Theory
    def kronecker(a, b)
      return a.abs != 1 ? 0 : 1 if b == 0
      return 0 if a.even? && b.even?

      k = 1

      case a & 7
      when 1, 7 then b    = b / 2     while b.even?
      when 3, 5 then b, k = b / 2, -k while b.even?
      end
      if b < 0
        b = -b
        k = -k if a < 0
      end

      while a != 0
        case b & 7
        when 1, 7 then a    = a / 2     while a.even?
        when 3, 5 then a, k = a / 2, -k while a.even?
        end
        k = -k if 2 & a & b == 2
        a = a.abs
        a, b = b % a, a
        a -= b if a > b / 2
      end

      return b > 1 ? 0 : k
    end
  end
end
