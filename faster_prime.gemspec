lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "faster_prime"

Gem::Specification.new do |spec|
  spec.name          = "faster_prime"
  spec.version       = FasterPrime::VERSION
  spec.authors       = ["Yusuke Endoh"]
  spec.email         = ["mame@ruby-lang.org"]

  spec.summary       = %q{A faster substitute for `lib/prime.rb`.}
  spec.description   = %q{This provides `Integer#prime?`, `Integer#prime_division`, and `Prime#each` that are almost compatible to and faster than `lib/prime.rb`.}
  spec.homepage      = "https://github.com/mame/faster_prime"

  spec.files         = `git ls-files -z`.split("\x0").reject do |f|
    f.match(%r{^(test|spec|features)/})
  end
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.16"
  spec.add_development_dependency "rake", "~> 10.0"
end
