# FasterPrime

This is an alternative implementation to the standard `lib/prime.rb`.  It is (almost) comapatible and faster than `lib/prime.rb`.

## Usage

* Install `faster_prime` gem.
* Replace `require "prime"` with `require "faster_prime"`.
* Replace `Prime` with `FasterPrime`.

If you want to use this without the core extension, you can use `require "faster_prime/base"` and `FasterPrime.prime?(num)` instead of `Integer#prime?`.

## Benchmark

### `Integer#prime?`

|test                  |     Prime|FasterPrime|
|:---------------------|---------:|----------:|
|`(2**31-1).prime?`    | 0.00355 s| 0.000035 s|
|`(10**100+267).prime?`| *n/a*    |   0.7003 s|

### `Integer#prime_division`

|test                       |     Prime|FasterPrime|
|:--------------------------|---------:|----------:|
|`(2**59-1).prime_division` |  0.1293 s|  0.01259 s|
|`(2**256+1).prime_division`| *n/a*    |    71.71 s|

### `Prime.each`

|test                       |     Prime|FasterPrime|
|:--------------------------|---------:|----------:|
|`Prime.each until 10**7`   |  0.6039 s|   0.4855 s|
|`Prime.each until 10**9`   |   57.21 s|    49.75 s|

## Internal

* `Integer#prime?` uses [APR-CL primality test](https://en.wikipedia.org/wiki/Adleman%E2%80%93Pomerance%E2%80%93Rumely_primality_test).
* `Integer#prime_division` uses [Pollard's rho algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants) and [Multiple Polynomial Quadratic Sieve](https://en.wikipedia.org/wiki/Quadratic_sieve#Multiple_polynomials).
* `Prime.each` uses [the sieve of Atkin](https://en.wikipedia.org/wiki/Sieve_of_Atkin).

## Caveat

* This library is not intended for practical use because it is written only in Ruby.  If you are serious about enjoying number theory, use [PARI/GP](https://pari.math.u-bordeaux.fr/).
* I have never evaluated the library throughout.  It may be slower than `lib/prime.rb` in some range.
* The author does not understand the algorithms.  I just read and implemented some papers.  If you find a bug, I'd appreciate if you could give me a patch :-)

## License

The MIT License (MIT)

Copyright (c) 2017 Yusuke Endoh

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
