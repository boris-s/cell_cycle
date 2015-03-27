#! /usr/bin/ruby
# encoding: utf-8

require 'minitest/autorun'
require 'sy'
require 'y_nelson' and include YNelson

describe :simple_cell_cycle do
  before do
    require './../lib/cell_cycle/simple'
  end

  it "should work" do
    CELL_CYCLE.must_be_kind_of YNelson::Net
    CELL_CYCLE.places.must_include A_phase
    CELL_CYCLE.places.must_include S_phase
    CELL_CYCLE.places.must_include Cdc20A
    S_phase_duration.must_equal 12.h
    S_phase_start.must_equal 5.h
    S_phase_end.must_equal 17.h
    A_phase_start.must_equal 3.h
    A_phase_end.must_equal 17.h
    Cdc20A_start.must_equal 22.h
    Cdc20A_end.must_equal 1.h
    CELL_CYCLE.A_tt.names.must_equal [ :A_phase_f, :S_phase_f, :Cdc20A_f ]
    # the simulation should look good
    sim = CELL_CYCLE.simulation time: 0..36.h.in( :s ),
                                step: 1.min.in( :s ),
                                sampling: 20.min.in( :s )
    sim.simulation_method.must_equal :pseudo_euler
    sim.run!
    sim.recording.marking( S_phase, A_phase, Cdc20A ).plot
    sleep 15
  end
end
