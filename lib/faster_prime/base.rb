require "faster_prime/utils"
require "faster_prime/primality_test"
require "faster_prime/prime_factorization"
require "faster_prime/sieve"
require "faster_prime/core_ext"

module FasterPrime
  VERSION = "1.0.1"

  module Core
    include Enumerable

    class SuccEnumerator < Enumerator
      alias succ next
    end

    def each(ubound = nil, &block)
      return SuccEnumerator.new(Float::INFINITY) do |y|
        Sieve.each(ubound) {|p| y << p }
      end unless block_given?

      Sieve.each(ubound, &block)
    end

    def prime?(n)
      raise ArgumentError, "Expected an integer, got #{ n }" unless n.respond_to?(:integer?) && n.integer?
      PrimalityTest.prime?(n)
    end

    def prime_division(n)
      PrimeFactorization.prime_factorization(n).to_a
    end

    def int_from_prime_division(pd)
      pd.inject(1) {|v, (p, e)| v * (p ** e) }
    end
  end

  def self.instance
    self
  end

  include Core
  extend Core
end
