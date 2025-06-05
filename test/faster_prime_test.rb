require "simplecov"
SimpleCov.start

require "test_helper"

class FasterPrimeTest < Minitest::Test
  def test_each
    primes = [
      2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
      71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113
    ]
    assert_equal(primes, FasterPrime.take(30))

    primes = [
      104743, 104759, 104761, 104773, 104779, 104789, 104801, 104803, 104827,
      104831, 104849, 104851, 104869, 104879, 104891, 104911, 104917, 104933,
      104947, 104953, 104959, 104971, 104987, 104999, 105019, 105023, 105031,
      105037, 105071, 105097
    ]
    assert_equal(primes, FasterPrime.lazy.drop(10000).take(30).to_a)
  end

  def test_each_ubound
    assert_equal(997, FasterPrime.each(1000).to_a.last)
  end

  def test_prime?
    mersenne_primes = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607]
    200.times do |n|
      assert_equal(mersenne_primes.include?(n), FasterPrime.prime?(2 ** n - 1))
    end
  end

  def test_prime_division
    100.times do |n|
      next if n == 0
      pe = FasterPrime.prime_division(2 ** n - 1)
      assert_equal(2 ** n - 1, FasterPrime.int_from_prime_division(pe))
      pe.each do |p,|
        assert(FasterPrime.prime?(p))
      end
    end
  end
end
