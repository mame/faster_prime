require "faster_prime/utils"

module FasterPrime
  module PrimeFactorization
    module_function

    # Factorize an integer
    def prime_factorization(n)
      return enum_for(:prime_factorization, n) unless block_given?

      return if n == 0
      if n < 0
        yield [-1, 1]
        n = -n
      end

      SMALL_PRIMES.each do |prime|
        if n % prime == 0
          c = 0
          begin
            n /= prime
            c += 1
          end while n % prime == 0
          yield [prime, c]
        end
        if prime * prime > n
          yield [n, 1] if n > 1
          return
        end
        return if n == 1
      end

      if PrimalityTest.prime?(n)
        yield [n, 1]
      else
        d = nil
        until d
          [PollardRho, MPQS].each do |algo|
            begin
              d = algo.try_find_factor(n)
            rescue Failed
            else
              break
            end
          end
        end

        pe = Hash.new(0)
        prime_factorization(n / d) {|p, e| pe[p] += e }
        prime_factorization(d) {|p, e| pe[p] += e }
        pe.keys.sort.each do |p|
          yield [p, pe[p]]
        end
      end
    end
  end

  class Failed < StandardError
  end

  # An implementation of Pollard's rho algorithm (Brent's variant)
  #
  # https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants
  module PollardRho
    module_function

    def self.try_find_factor(n)
      x0 = rand(n)
      y = x0
      m = Math.log(n, 2).ceil
      c = 2
      r = q = 1
      while true
        x = y
        r.times { y = (y * y + c) % n }
        k = 0
        while true
          ys = y
          [m, r - k].min.times do
            y = (y * y + c) % n
            q = (q * (x - y)) % n
          end
          g = q.gcd(n)
          k += m
          break if k >= r || g > 1
        end
        r *= 2
        break if g > 1
      end
      if g == n
        while true
          ys = (ys * ys + c) % n
          g = (x - ys).gcd(n)
          break if g > 1
        end
      end
      if g == n
        raise Failed, "failed to find a factor: #{ n }"
      else
        g
      end
    end
  end

  # An implementation of MPQS factoring algorithm.
  #
  # R. D. Silverman,
  # "The Multiple Polynomial Quadratic Sieve." (1987)
  # http://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866119-8/S0025-5718-1987-0866119-8.pdf

  class MPQS
    # MPQS parameters
    # (N digits, factor base size, sieve interval, and large prime tolerance)
    PARAMETERS = [
      [24,  100,    250, 1.5 ],
      [30,  200,  1_250, 1.5 ],
      [36,  400,  1_250, 1.75],
      [42,  900,  2_500, 2.0 ],
      [48, 1200,  7_000, 2.0 ],
      [54, 2000, 12_500, 2.2 ],
      [60, 3000, 17_500, 2.4 ],
      [66, 4500, 25_000, 2.6 ],
    ]

    SMALL_PRIMES = [2, 3, 5, 7, 11]

    # Generates a table for finding k such that kN = 1 mod 4 and that kN is
    # rich in small quadratic residues, i.e., (kN | p) = 1 for all p in
    # SMALL_PRIMES.
    def self.k_table
      return @@k_table if defined?(@@k_table)

      k_table = { 1 => {}, 3 => {} }
      Sieve.each do |p|
        next if p <= SMALL_PRIMES.last
        state = SMALL_PRIMES.map {|q| Utils.kronecker(p, q) }
        k_table[p % 4][state] ||= p
        break if k_table.all? {|r, t| t.size == 2 ** SMALL_PRIMES.size }
      end

      @@k_table = k_table
    end

    def self.try_find_factor(n)
      new(n).try_find_factor
    end

    S = 65536

    def initialize(n)
      @n = n

      if @n <= 1 || PrimalityTest.prime?(@n)
        @trivial_factor = 1
        return
      end
      SMALL_PRIMES.each do |p|
        if @n % p == 0
          @trivial_factor = p
          return
        end
      end

      m = Utils.integer_square_root(n)
      if m * m == n
        @trivial_factor = m
        return
      end

      # select parameters
      # @k: a multiplier
      @k = MPQS.k_table[@n % 4][SMALL_PRIMES.map {|q| Utils.kronecker(@n, q) }]
      @kn = @k * @n

      # @f: the size of the factor base
      # @m: the length of the sieve interval (2M+1)
      # @t: the large prime tolerance
      _digit, @f, @m, @t =
        PARAMETERS.find {|d,| Math.log(@n, 10) < d } || PARAMETERS.last

      # compute factor base
      # @ps: the factor base
      # @log_ps: log of each factor base
      # @sqrts: sqrt(kN) mod p
      @ps, @log_ps, @sqrts = [], [], []
      Sieve.each do |p|
        if @n % p == 0
          @trivial_factor = p
          return
        end
        if p >= 3 && Utils.kronecker(@kn, p) == 1
          @ps << p
          @log_ps << (Math.log(p) * S).floor
          @sqrts << Utils.mod_sqrt(@kn, p)
          break if @ps.size >= @f
        end
      end
      @pmax = @ps.last # max factor base

      # for enumerating coefficients
      @base_d = Math.sqrt(Math.sqrt(@kn) / (Math.sqrt(2) * @m)).floor | 3
      @diff_d = 0 # 0, 4, -4, 8, -8, 12, -12, ...

      # sieve table
      @ws = [0] * (2 * @m + 1)

      # a test value for sieve
      @test_value = (Math.log(@m * Math.sqrt(@kn) / @pmax ** @t) * S).round

      # relations found
      @incomplete_relations = {}

      # for gaussian elimination
      @elim = GaussianElimination.new(@ps.size + 2)
    end

    # try to find any factor (not necessarily a prime)
    def try_find_factor
      return @trivial_factor if defined?(@trivial_factor)

      catch(:factor_found) do
        while true
          p :foo
          select_q
          solve_roots
          sieve
          @ws.each_with_index do |w, i|
            check_relation(i - @m) if w > @test_value
          end
        end
      end
    end

    # Selects Q(x) = Ax^2 + Bx + C for quadratic sieve
    def select_q
      # find a prime D such that (kN | D) = 1 and D = 3 mod 4
      begin
        @d = @base_d + @diff_d
        @diff_d = @diff_d > 0 ? -@diff_d : 4 - @diff_d
      end until Utils.kronecker(@kn, @d) == 1 && PrimalityTest.probable_prime?(@d) != false && !@ps.include?(@d)

      # select coefficients (see the paper)
      @a = @d * @d

      h0 = Utils.mod_pow(@kn, @d / 4, @d)
      h1 = (@kn * h0) % @d
      # assert: h1*h1 = kN mod D
      # assert: h0*h1 = 1 mod D
      h2 = (Utils.mod_inv(2, @d) * h0 * ((@kn - h1 * h1) / @d)) % @d
      # assert: kN - h1^2 = 0 mod D

      @b = (h1 + h2 * @d) % @a
      @b -= @a if @b[0] == 0
      # assert: B^2 = kN mod 4A

      #@c = (@b * @b - @kn) / (4 * @a) # not needed

      # compute real roots of Q(x) = 0
      s = Utils.integer_square_root(@kn)
      @real_root_1 = (-@b - s) / (2 * @a)
      @real_root_2 = (-@b + s) / (2 * @a)
      # x is between real roots if @real_root_1 <= x && x <= @real_root_2 
    end

    # Computes the roots of Q(x) = 0 mod p (p in the factor base)
    def solve_roots
      @roots = []
      @ps.zip(@sqrts) do |p, s|
        a_inv = Utils.mod_inv(2 * @a, p)
        x1 = (-@b + s) * a_inv % p
        x2 = (-@b - s) * a_inv % p
        @roots << [x1, x2]
      end
    end

    # Performs the quadratic sieve
    def sieve
      @ws.fill(0)
      @ps.zip(@log_ps, @roots) do |p, lp, (s1, s2)|
        ((@m + s1) % p).step(2 * @m, p) {|i| @ws[i] += lp }
        ((@m + s2) % p).step(2 * @m, p) {|j| @ws[j] += lp }
      end
    end

    # Tests whether Q(x) is actually a product of factor bases
    def check_relation(x)
      h, qx = compute_q(x)
      return if qx == 0
      es, l = exponent_bitvector(qx)

      # discard this x if the residue L is too big
      return if l > @pmax ** @t

      if l == 1
        # complete relation found:
        # Q(x) = p0^e0 * p1^e1 * ... * pk^ek (pi in the factor base)
        qx_vec = Hash.new(0)
        PrimeFactorization.prime_factorization(qx) {|p, e| qx_vec[p] += e }
        collect_relation(es, h, qx_vec)

      elsif @incomplete_relations[l]
        # large prime procedure:
        # make a complete relation by multiplying two incomplete relations
        es2, h2, qx2 = @incomplete_relations[l]

        # XXX: use FactoredInteger
        qx_vec = Hash.new(0)
        PrimeFactorization.prime_factorization(qx) {|p, e| qx_vec[p] += e }
        PrimeFactorization.prime_factorization(qx2) {|p, e| qx_vec[p] += e }

        collect_relation(es ^ es2, h * h2 % @kn, qx_vec)

      else
        @incomplete_relations[l] = [es, h, qx]
      end
    end

    # Adds a complete relation found to gaussian elimination list
    def collect_relation(es, h, q)
      #p "%0#{ @ps.size + 1 }b" % es
      @elim[es] = [h, q]

      if @elim.size > 0.9 * @f && @elim.size % [1, @f / 100].max == 0
        @elim.eliminate do |relations|
          # factor candidate found
          check_factor_candidate(relations)
        end
      end

      if @elim.size > @f * 1.5
        raise Failed, "failed to find a factor: #{ @n }"
      end
    end

    # Computes and checks a factor candidate
    def check_factor_candidate(relations)
      # computes the factor candidate
      p1 = 1
      p2_vec = Hash.new(0) # XXX: use FactoredInteger
      relations.each do |h, qx_vec|
        p1 = p1 * h % @kn
        qx_vec.each {|p, e| p2_vec[p] += e }
      end
      p2 = 1
      p2_vec.each {|p, e| p2 = (p2 * Utils.mod_pow(p, e / 2, @kn)) % @kn }

      return if p1 == p2 || (p1 + p2) % @kn == 0

      factor = (p1 + p2).gcd(@n)

      # check the factor candidate
      if factor != 1 && factor != @n
        throw :factor_found, factor
      end
    end

    # Computes Q(x)
    def compute_q(x)
      h = (2 * @a * x + @b) * Utils.mod_inv(2 * @d, @kn) % @kn

      q = h * h % @kn
      q -= @kn if @real_root_1 <= x && x <= @real_root_2
      # assert: q == @a*x*x + @b*x + @c

      return h, q
    end

    # Factorizes n in factor base and returns an exponent bitvector
    def exponent_bitvector(n)
      vec = 0
      if n < 0
        vec = 1
        n = -n
      end
      bit = 2
      ([2] + @ps).each do |p|
        e = false
        while n % p == 0
          n /= p
          e = !e
        end
        vec |= bit if e
        bit <<= 1
      end
      return vec, n
    end


    # Incremental gaussian elimination.
    #
    # A. K. Lenstra, and M. S. Manasse,
    # "Compact incremental Gaussian Elimination over Z/2Z." (1988)
    # http://dl.acm.org/citation.cfm?id=896266
    class GaussianElimination
      def initialize(m)
        @n = 0
        @m = m
        @d = []
        @e = []
        @u = []
        @v = []
      end

      def []=(es, val)
        @d << es
        @e << nil
        @v << val
      end

      def size
        @d.size
      end

      def eliminate
        n2 = @d.size

        @n.times do |i|
          if @u[i]
            @e[@u[i]] = i
            @n.upto(n2 - 1) do |j|
              @d[j] ^= @d[i] if @d[j][@u[i]] == 1
            end
          end
        end

        @n.upto(n2 - 1) do |j|
          found = false

          @m.times do |i|
            if @d[j][i] == 1 && !@e[i]
              @u[j] = i
              @e[@u[j]] = j
              @d[j] &= ~(1 << @u[j])
              (j + 1).upto(n2 - 1) do |k|
                @d[k] ^= @d[j] if @d[k][@u[j]] == 1
              end
              found = true
              break
            end
          end

          if !found
            # linear dependency found
            @u[j] = nil
            c = {}
            c[j] = true
            @m.times do |k|
              c[@e[k]] = true if @d[j][k] == 1
            end
            yield c.keys.map {|k| @v[k] }
          end
        end

        @n = n2
      end
    end
  end
end
