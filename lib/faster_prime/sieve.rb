# Sieve of Atkin
# http://cr.yp.to/papers/primesieves.pdf

require "faster_prime/utils"

module FasterPrime
  module Sieve
    module_function

    TBL = []
    [
      [:mark_type_1, ->(f, g) { 4*f*f + g*g }, 15, 30,  1,  4],
      [:mark_type_2, ->(f, g) { 3*f*f + g*g }, 10, 30,  7, 12],
      [:mark_type_3, ->(f, g) { 3*f*f - g*g }, 10, 30, 11, 12],
    ].each do |type, expr, fmax, gmax, n, mod|
      gmax.times do |g|
        fmax.times do |f|
          d = expr[f, g] % 60
          TBL << [type, f, g, d] if d % mod == n && d.gcd(60) == 1
        end
      end
    end

    MODULO60 = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]

    def each(ubound = nil, &blk)
      if ubound && SMALL_PRIMES.last >= ubound
        SMALL_PRIMES.each do |n|
          break if n > ubound
          yield n
        end
        return
      end

      u = SMALL_PRIMES.last / 60 * 60
      SMALL_PRIMES.pop until SMALL_PRIMES.last < u
      SMALL_PRIMES.each(&blk)

      b = Utils.integer_square_root(ubound * 3000) if ubound
      b = 100000 if !b || b > 100000
      l = u / 60

      # Algorithm 3.1, 3.2, 3.3
      flags = [0] * b
      while true
        # step
        flags.fill(0)

        # sieve
        TBL.each do |name, f, g, d|
          send(name, flags, f, g, d, l, b)
        end

        # squarefree
        squarefree(flags, l, b)

        # enumerate primes
        if ubound && ubound < 60 * (l + b)
          enumerate(flags, l, b) do |p|
            return if p >= ubound
            yield p
          end
        else
          enumerate(flags, l, b, &blk)
        end

        l += b
      end
    end

    def enumerate(flags, l, b)
      b.times do |x|
        base = 60 * (l + x)
        flag = flags[x]
        MODULO60.each do |d|
          next if flag[d] == 0
          p = base + d
          SMALL_PRIMES << p
          yield p
        end
      end
    end

    def squarefree(flags, l, b)
      SMALL_PRIMES.each do |q|
        q2 = q * q
        break if q2 > 60 * (l + b)
        next if q < 7
        minus_inv_60 = -Utils.mod_inv(60, q2)
        MODULO60.each do |d|
          # 60(l+x)+d == 0 (mod q2)
          # x = -d * 60^{-1} - l (mod q2)
          x = (d * minus_inv_60 - l) % q2
          m = ~(1 << d)
          while x < b
            flags[x] &= m
            x += q2
          end
        end
      end
    end

    # Algorithm 4.1
    def mark_type_1(flags, f, g, d, l, b)
      x, y0 = f, g
      k0 = (4 * f * f + g * g - d) / 60 - l
      m = 1 << d
      while k0 < l + b
        k0 += 2*x + 15
        x += 15
      end
      while true
        x -= 15
        k0 -= 2*x + 15
        return if x <= 0
        while k0 < 0
          k0 += y0 + 15
          y0 += 30
        end
        k, y = k0, y0
        while k < b
          flags[k] ^= m # mark
          k += y + 15
          y += 30
        end
      end
    end

    # Algorithm 4.2
    def mark_type_2(flags, f, g, d, l, b)
      x, y0 = f, g
      k0 = (3 * f * f + g * g - d) / 60 - l
      m = 1 << d
      while k0 < b
        k0 += x + 5
        x += 10
      end
      while true
        x -= 10
        k0 -= x + 5
        return if x <= 0
        while k0 < 0
          k0 += y0 + 15
          y0 += 30
        end
        k, y = k0, y0
        while k < b
          flags[k] ^= m # mark
          k += y + 15
          y += 30
        end
      end
    end

    # Algorithm 4.3
    def mark_type_3(flags, f, g, d, l, b)
      x, y0 = f, g
      k0 = (3 * f * f - g * g - d) / 60 - l
      m = 1 << d
      while true
        while k0 >= b
          return if x <= y0
          k0 -= y0 + 15
          y0 += 30
        end
        k, y = k0, y0
        while k >= 0 && y < x
          flags[k] ^= m # mark
          k -= y + 15
          y += 30
        end
        k0 += x + 5
        x += 10
      end
    end
  end
end
