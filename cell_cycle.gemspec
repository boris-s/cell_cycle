# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'cell_cycle/version'

Gem::Specification.new do |spec|
  spec.name          = "cell_cycle"
  spec.version       = CellCycle::VERSION
  spec.authors       = ["boris"]
  spec.email         = ["\"boris@iis.sinica.edu.tw\""]
  spec.summary       = %q{A model of eukaryotic cell cycle.}
  spec.description   = %q{Eukaryotic cell cycle modelled at different levels of precision. Has two levels at the moment: Simple and Virginia Tech.}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.6"
  spec.add_development_dependency "rake"

  spec.add_dependency "y_nelson"
  spec.add_dependency "sy"
end
