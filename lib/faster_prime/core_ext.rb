class Integer
  def Integer.from_prime_division(pd)
    FasterPrime.int_from_prime_division(pd)
  end

  def prime_division(generator = nil)
    FasterPrime.prime_division(self)
  end

  def prime?
    FasterPrime.prime?(self)
  end

  def Integer.each_prime(ubound, &block) # :yields: prime
    Faster.Prime.each(ubound, &block)
  end
end
