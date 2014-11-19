# CellCycle

Cell cycle model.

## Installation

You can install it by

    $ gem install cell_cycle

## Usage

This gem contains two versions of cell cycle. Both of them contain an engine,
a Petri net whose execution results in cycling, and 3 output places: A_phase,
S_phase and Cdc20A.

Simple version of the cell cycle is called by:

    require 'cell_cycle/simple'

You can then try

    s = CELL_CYCLE.default_simulation
    s.run!
    s.recording.print
    s.recording.plot except: [ Timer ]        # relies on Gnuplot

And a more complex one, developed at Virginia Institute of Technology, and called by:

    require 'cell_cycle/virginia_tech'

You can then try

    run!

The current version of the Virginia Tech cell cycle takes a long time to execute, but
once the waiting is over, you can display the results

## Contributing

1. Fork it ( https://github.com/[my-github-username]/cell_cycle/fork )
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create a new Pull Request
