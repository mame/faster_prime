require "faster_prime/utils"

module FasterPrime
  module PrimalityTest
    module_function

    DEFAULT_RANDOM = Random.new

    def prime?(n, random = DEFAULT_RANDOM)
      r = probable_prime?(n, random)
      return r if r != nil
      return APRCL.prime?(n)
    end

    # Miller-Rabin primality test
    # Returns true for a prime, false for a composite, nil if unknown.
    def probable_prime?(n, random = Random.new)
      return n >= 2 if n <= 3
      return true  if n == 2 || n == 3
      return true if n == 5
      return false unless 30.gcd(n) == 1
      return SMALL_PRIME_TABLE[n / 2] if n / 2 < SMALL_PRIME_TABLE.size
      m = n % 6
      return false if m != 1 && m != 5

      list = DETERMINISTIC_MR_TEST_TABLE.find {|t,| n < t }
      if list
        ret = true
        list = list.last
      else
        ret = nil
        list = (1..20).map { random.rand(n - 2) + 1 }
      end

      d, s = n - 1, 0
      while d.even?
        d >>= 1
        s += 1
      end
      list.each do |a|
        y = Utils.mod_pow(a, d, n)
        if y != 1
          f = s.times do
            break if y == n - 1
            y = (y * y) % n
          end
          return false if f
        end
      end

      ret
    end

    # These numbers are taken from:
    #   http://code.google.com/p/sympy/
    #   http://primes.utm.edu/prove/prove2_3.html
    DETERMINISTIC_MR_TEST_TABLE = [
      [          1_373_653, [2, 3]],
      [        170_584_961, [350, 3_958_281_543]],
      [      4_759_123_141, [2, 7, 61]],
      [     75_792_980_677, [2, 379_215, 457_083_754]],
      [  1_000_000_000_000, [2, 13, 23, 1_662_803]],
      [  2_152_302_898_747, [2, 3, 5, 7, 11]],
      [  3_474_749_660_383, [2, 3, 5, 7, 11, 13]],
      [341_550_071_728_321, [2, 3, 5, 7, 11, 13, 17]],
    ]

    # APRCL Primarity Test
    #
    # Henry Cohen. A course in computational algebraic number theory.
    # 9.1 The Jacobi Sum Test

    class APRCL
      # usage:
      #
      # (1) simplest
      #
      #   APRCL.prime?(N) #=> true or false
      #
      # (2) reuse-table
      #
      #   aprcl = APRCL.new(B)
      #   aprcl.prime?(N1) #=> true or false
      #   aprcl.prime?(N2) #=> true or false
      #   aprcl.bound #=> upper bound (>= B)
      #
      # (3) manual-t (for test or debug)
      #
      #   aprcl = APRCP.new(B)
      #   aprcl.set_t(t)
      #   aprcl.prime?(N) #=> true or false

      # An internal exception for terminating the algorithm because n is composite
      class Composite < StandardError
      end

      # An exception representing that the test has failed (highly improbable)
      class Failed < StandardError
      end

      # Returns true if n is prime, false otherwise.
      def self.prime?(n)
        new(n).prime?
      end

      # Creates an APRCL instance to test an integer that is less than bound
      def initialize(n, jacobi_sum_table: DefaultJacobiSumTable, t: nil)
        # an integer in question
        @n = n

        # precompute n-1 and (n-1)/2 because n is normally big
        @n1 = n - 1
        @n1_2 = @n1 / 2

        # @eta[pk]: a set of p^k-th roots of unity
        @eta = {}

        # a cache of Jacobi sums
        @js = jacobi_sum_table

        # compute e(t)
        if t
          raise "t must be even" if t.odd?
          @t = t
          @et = compute_e(@t)
          raise "t is too small to test n" if @n >= @et ** 2
        else
          @t = find_t(@n)
          @et = compute_e(@t)
        end
      end

      # Algorithm 9.1.28 (Jacobi Sum Primality test)
      def prime?
        return false if @n <= 1

        # 1. [Check GCD]
        g = @n.gcd(@t * @et)
        return g == @n && PrimalityTest.prime?(g) if g != 1

        # 2. [Initialize]
        lp = {}
        PrimeFactorization.prime_factorization(@t) do |p, |
          lp[p] = false if p == 2 || Utils.mod_pow(@n, p - 1, p ** 2) == 1
        end

        # 3. [Loop on characters]
        Utils.each_divisor(@t) do |d|
          q = d + 1
          next unless PrimalityTest.prime?(q)
          PrimeFactorization.prime_factorization(d) do |p, k|
            # 4. [Check (*beta)]
            lp.delete(p) if step4(q, p, k)
          end
        end

        # 5. [Check conditions Lp]
        lp.keys.each {|p| step5(p) }

        # 6. [Final trial division]
        r = 1
        1.upto(@t - 1) do
          r = r * @n % @et
          return false if @n % r == 0 && 1 < r && r < @n
        end

        return true

      rescue Composite
        return false
      end


      private

      def step4(q, p, k)
        case
        when p >= 3 then step4a(q, p, k)
        when k >= 3 then step4b(q, k)
        when k == 2 then step4c(q)
        when k == 1 then step4d(q)
        end
      end

      # p is odd
      def step4a(q, p, k)
        pk = p ** k

        # s1 = J(p,q)^theta, s3 = J(p,q)^alpha
        s1 = s3 = ZZeta.one(p, k)
        r = @n % pk
        jpq = @js.get_j(p, q)
        1.upto(pk - 1) do |x|
          next if x % p == 0
          j = jpq.substitute(x).reduce_mod!(@n)
          s1 = s1 * j.mod_pow(x, @n)
          s3 = s3 * j.mod_pow(r * x / pk, @n)
        end

        s2 = s1.mod_pow(@n / pk, @n)

        # S(p,q) = s2 J(p,q)^alpha
        s = (s2 * s3).reduce_mod!(@n)

        i = check_root_of_unity(p, k, pk, s)

        i % p != 0 # it is a primitive root of unity
      end

      # p = 2, k >= 3
      def step4b(q, k)
        pk = 2 ** k

        # s1 = J_3(q)^theta, s3 = J_3(q)^alpha
        s1 = s3 = ZZeta.one(2, k)
        r = @n % pk
        j3, j2 = @js.get_j3_j2(q)
        diff = 1
        0.step(pk - 1, 4) do |x|
          x, diff = x + diff, -diff # x = 8m+1 or 8m+3, i.e., 1, 3, 9, 11, ...
          j = j3.substitute(x).reduce_mod!(@n)
          s1 = (s1 * j.mod_pow(x         , @n)).reduce_mod!(@n)
          s3 = (s3 * j.mod_pow(r * x / pk, @n)).reduce_mod!(@n)
        end

        # S(2,q) = s2 J_3(q)^alpha        if r = 8m+1 or 8m+3
        #          s2 J_3(q)^alpha J_2(q) if r = 8m+5 or 8m+7
        s = s3 * s1.mod_pow(@n / pk, @n)
        s = s * j2 if r[2] == 1 # r = 8m+5 or 8m+7
        s.reduce_mod!(@n)

        i = check_root_of_unity(2, k, pk, s)

        i[0] == 1 && Utils.mod_pow(q, @n1_2, @n) == @n1
      end

      # p = 2, k = 2
      def step4c(q)
        s0 = @js.get_j(2, q).sqr
        s1 = s0.scalar(q)
        s2 = s1.mod_pow(@n / 4, @n)
        # S(2,q) = s2          if n = 4m+1
        #          s2 J(2,q)^2 if n = 4m+3
        s = (@n[1] == 0 ? s2 : s2 * s0).reduce_mod!(@n)

        i = check_root_of_unity(2, 2, 4, s)

        i[0] != 1 && Utils.mod_pow(q, @n1_2, @n) == @n1
      end

      # p = 2, k = 1
      def step4d(q)
        # S(2,q) = (-q)^((n-1)/2) mod N
        s = Utils.mod_pow(-q, @n1_2, @n)

        raise Composite if s != 1 && s != @n1

        s == @n && @n[1] == 0
      end

      # check whether there exists a p^k-th root of unity e such that s = e mod n
      def check_root_of_unity(p, k, pk, s)
        unless @eta[pk]
          table = {}
          pk.times do |i|
            table[ZZeta.monomial(p, k, i).reduce_mod!(@n)] = i
          end
          @eta[pk] = table
        end

        @eta[pk][s] || raise(Composite)
      end


      TRIAL = 2000
      def step5(p)
        bound_k = 4

        # Find q and k by brute force.
        # We first try to find a small k because step4 might take long time
        # when k is big (>= 4).  (An example case: prime?(19541) with TRIAL = 128)
        # If we cannot do so, we use any k.
        2.times do
          q = p + 1
          TRIAL.times do
            if @et % q != 0 && PrimalityTest.prime?(q)
              if @n % q == 0
                return if @n == q
                raise Composite
              end
              k = Utils.valuation(q - 1, p)
              return if k <= bound_k && step4(q, p, k)
            end
            q += p
          end
          bound_k = Float::INFINITY
        end

        raise Failed, "APRCL primality test failed (highly improbable)"
      end


      # The pairs of t candidate and floor(log_2(e^2(t))).
      # generated by tool/build-aprcl-t-table.rb
      T_TABLE = [
        [12, 31], [24, 33], [30, 34], [36, 54], [60, 65], [72, 68], [108, 70],
        [120, 78], [180, 102], [240, 104], [360, 127], [420, 136], [540, 153],
        [840, 165], [1008, 169], [1080, 178], [1200, 192], [1260, 206],
        [1620, 211], [1680, 222], [2016, 225], [2160, 244], [2520, 270],
        [3360, 279], [3780, 293], [5040, 346], [6480, 348], [7560, 383],
        [8400, 396], [10080, 426], [12600, 458], [15120, 527], [25200, 595],
        [30240, 636], [42840, 672], [45360, 684], [55440, 708], [60480, 771],
        [75600, 775], [85680, 859], [100800, 893], [110880, 912], [128520, 966],
        [131040, 1009], [166320, 1042], [196560, 1124], [257040, 1251],
        [332640, 1375], [393120, 1431], [514080, 1483], [655200, 1546],
        [665280, 1585], [786240, 1661], [831600, 1667], [917280, 1677],
        [982800, 1728], [1081080, 1747], [1179360, 1773], [1285200, 1810],
        [1310400, 1924], [1441440, 2001], [1663200, 2096], [1965600, 2166],
        [2162160, 2321], [2751840, 2368], [2827440, 2377], [3326400, 2514],
        [3341520, 2588], [3603600, 2636], [3931200, 2667], [4324320, 3028],
        [5654880, 3045], [6652800, 3080], [6683040, 3121], [7207200, 3283],
        [8648640, 3514], [10810800, 3725], [12972960, 3817], [14414400, 3976],
        [18378360, 3980], [21621600, 4761], [36756720, 5067], [43243200, 5657],
        [64864800, 5959], [73513440, 6423], [86486400, 6497],
      ]

      # Find t such that e(t) > n^(1/2).
      def find_t(n)
        n_log2 = n.bit_length

        # e^2(t) >= 2^(et2_log2) >= 2^n_log2 > n
        t, = T_TABLE.find {|_t, et2_log2| et2_log2 >= n_log2 }
        raise "too large" unless t

        t
      end

      # Compute e(t) (the definition comes from Theorem 9.1.12)
      def compute_e(t)
        r = 2
        Utils.each_divisor(t) do |d|
          q = d + 1
          r *= q ** (Utils.valuation(t, q) + 1) if PrimalityTest.prime?(q)
        end
        r
      end


      class JacobiSumTableClass
        # Algorithm 9.1.27 (Precomputations)
        # Note that this computes the table lazily, i.e., on demand

        def initialize
          @f = {} # a table of the function f(x) in 2. (1)
          @j = {} # J(p,q)
          @j3_j2 = {} # J3(p,q) and J2(p,q)
        end

        def clear
          @f.clear
          @j.clear
          @j3_j2.clear
        end

        # compute a table of the function f(x) such that 1 - g^x = g^f(x),
        # where g is a primitive root modulo q
        def calc_f(q)
          return @f[q] if @f[q]

          g = Utils.primitive_root(q)
          f, h = [], {}
          1.upto(q - 2) {|x| h[(1 - Utils.mod_pow(g, x, q)) % q] = x }
          1.upto(q - 2) {|fx| f[h[Utils.mod_pow(g, fx, q)]] = fx }
          @f[q] = f
        end

        # compute J(p,q)
        def get_j(p, q)
          # assumes that p >= 3 or (p == 2 && k >= 2)
          return @j[[p, q]] if @j[[p, q]]

          f = calc_f(q)
          k = Utils.valuation(q - 1, p)
          pk = p ** k

          # J(p,q) = \Sum_{1<=x<=q-2} (zeta_{p^k}^{x+f(x)})
          coefs = [0] * pk
          1.upto(q - 2) {|x| coefs[(x + f[x]) % pk] += 1 }
          @j[[p, q]] = ZZeta.new(coefs, p, k).reduce!
        end

        # compute J_3(q) and J_2(q)
        def get_j3_j2(q)
          # assumes that p == 2 && k >= 3
          return @j3_j2[q] if @j3_j2[q]

          f = calc_f(q)
          k = Utils.valuation(q - 1, 2)
          pk = 2 ** k

          # coefs3: \Sum_{1<=x<=q-2} (zeta_{2^k}^{2x+f(x)})
          # coefs2: \Sum_{1<=x<=q-2} (zeta_{2^3}^{3x+f(x)})
          coefs3 = [0] * pk
          coefs2 = [0] * pk
          1.upto(q - 2) do |x|
            coefs3[(2 * x + f[x]) % pk] += 1
            coefs2[(3 * x + f[x]) % 8 * (pk / 8)] += 1
          end
          # J_3(q) = J(2,q) coefs3, J_2(q) = ( coefs2 )^2
          j3 = (get_j(2, q) * ZZeta.new(coefs3, 2, k)).reduce!
          j2 = ZZeta.new(coefs2, 2, k).reduce!.sqr.reduce!
          @j3_j2[q] = [j3, j2]
        end
      end
      DefaultJacobiSumTable = JacobiSumTableClass.new


      # A group algebra Z[zeta_{p^k}]
      # (linear sum of p^k-th roots of unity)
      class ZZeta
        def initialize(coefs, p, k)
          @coefs = coefs
          @p, @k = p, k
          @pk1 = p ** (k - 1) # p^(k-1)
          @pk = @pk1 * p      # p^k
          @p1pk1 = @pk - @pk1 # (p-1) p^(k-1)
          if @pk != @coefs.size
            raise "coefs.size (#{ coefs.size }) must be equal to #{ @p }**#{ @k }"
          end
        end

        attr_reader :coefs
        protected :coefs

        # generates 1
        def self.one(p, k)
          @@one[p ** k] ||= monomial(p, k, 0)
        end
        @@one = {}

        # generates zeta_{p^k}^i
        def self.monomial(p, k, i)
          coefs = [0] * (p ** k)
          coefs[i] = 1
          ZZeta.new(coefs, p, k)
        end

        # for hashing
        def ==(other); @coefs == other.coefs; end
        def hash; @coefs.hash; end
        alias eql? ==

        # scalar multiplication
        def scalar(n)
          ZZeta.new(@coefs.map {|x| x * n }, @p, @k)
        end

        # polynomial multiplication
        def *(other)
          coefs = [0] * @pk
          other.coefs.each_with_index do |oc, j|
            next if oc == 0
            @pk.times do |i|
              coefs[(i + j) % @pk] += @coefs[i] * oc
            end
          end
          ZZeta.new(coefs, @p, @k)
        end

        # polynomial power modulo n
        def mod_pow(exp, n)
          return ZZeta.one(@p, @k) if exp == 0

          r = nil
          z = reduce_mod!(n)
          len = exp.bit_length
          len.times do |i|
            if exp[i] == 1
              r = r ? (z * r).reduce_mod!(n) : z
              return r if i == len - 1
            end
            z = z.sqr.reduce_mod!(n)
          end

          raise "assert not reached"
        end

        # polynomial square (self * self)
        # assumes self is already reduced
        def sqr
          coefs = [0] * @pk
          @p1pk1.times do |j|
            oc = @coefs[j]
            next if oc == 0
            @p1pk1.times do |i|
              coefs[(i + j) % @pk] += @coefs[i] * oc
            end
          end
          ZZeta.new(coefs, @p, @k)
        end

        # generate a reminder divided by CyclotomicPolynomial(p^k)
        def reduce!
          # CyclotomicPolynomial(n) = \Prod_{d|n} (x^{n/d} - 1)^\moeb{d}
          #
          # CyclotomicPolynomial(p^k) =
          #   (x^{p^ k   } - 1)^\moeb{p^0}  (case d = p^0)
          # * (x^{p^{k-1}} - 1)^\moeb{p^1}  (case d = p^1)
          # ...
          # * (x^{p^ 0   } - 1)^\moeb{p^k}  (case d = p^k)
          # = (x^{p^k} - 1) / (x^{p^{k-1}} - 1)
          # = \Sum_{0 <= i < p} (x^{i p^{k-1}})
          @p1pk1.times do |i|
            @coefs[i] -= @coefs[@p1pk1 + (i % @pk1)]
          end
          @pk1.times do |i|
            @coefs[@p1pk1 + i] = 0
          end
          self
        end

        # reduce modulo n
        def reduce_mod!(n)
          @pk.times do |i|
            # I would use centermod if it was built-in...
            @coefs[i] = (@coefs[i] - @coefs[@p1pk1 + (i % @pk1)]) % n
          end
          self
        end

        # replace zeta with zeta^x
        def substitute(x)
          coefs = []
          @pk.times {|i| coefs << @coefs[(x * i) % @pk] }
          ZZeta.new(coefs, @p, @k)
        end

        def inspect
          s = []
          (@pk - 1).downto(0) do |i|
            c = @coefs[i]
            next if c == 0
            s << (c > 0 ? ?+ : ?-) +
              case
              when i == 0 then c.abs.to_s
              when c.abs == 1 then i == 1 ? "z" : "z^#{ i }"
              else "#{ c.abs }#{ i == 1 ? "z" : "z^#{ i }" }"
              end
          end
          return "0" if s.empty?
          s.first[0, 1] = "" if s.first[0] == ?+
          s.join
        end
      end
    end
  end
end
