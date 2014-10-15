# encoding: utf-8

# A simplistic cell cycle. It has a single input (Timer place) and 3 outputs
# (A_phase, S_phase, Cdc20A). A_phase is the phase when the cell cycle enzyme
# machinery is synthesized. S_phase has the standard meaning: DNA synthesis
# phase. Cdc20A represents the APC (Anaphase Promoting Complex). When present,
# it degrades the cell cycle enzyme machinery.

require 'y_nelson' and include YNelson
require 'sy'
require './../cell_cycle'

# Constants that control the cell cycle settings.
S_phase_duration = 12.h
S_phase_start = 5.h
S_phase_end = S_phase_start + S_phase_duration
A_phase_start = 3.h
A_phase_end = S_phase_end
Cdc20A_start = 22.h
Cdc20A_end = 1.h

# Alternative setting for shorter cycle (4.h).
=begin
# Constants that control the cell cycle settings.
S_phase_duration = 4.h
S_phase_start = 1.h + 40.min
S_phase_end = S_phase_start + S_phase_duration
A_phase_start = 1.h
A_phase_end = S_phase_end
Cdc20A_start = 7.h + 30.min
Cdc20A_end = 20.min
=end

# Figure them out as numbers in seconds.
Sα = S_phase_start.in :s
Sω = S_phase_end.in :s
Aα = A_phase_start.in :s
Aω = A_phase_end.in :s
Cdc20Aα = Cdc20A_start.in :s
Cdc20Aω = Cdc20A_end.in :s

# Timer place
Timer = Place m!: 0

# The clock transition
Clock = Transition stoichiometry: { Timer: 1 }, rate: 1

# Empirical places (in arbitrary units); output of the cell cycle.
A_phase = Place m!: 0
S_phase = Place m!: 0
Cdc20A = Place m!: 1

# Include them in the CELL_CYCLE net.
CELL_CYCLE = << Timer << Clock << A_phase << S_phase << Cdc20A

# Assignment transitions that control the state of the places A_phase, S_phase
# and Cdc20A.
# 
A_phase_ϝ = Transition assignment: -> t { t > Aα && t < Aω ? 1 : 0 },
                       domain: Timer,
                       codomain: A_phase
S_phase_ϝ = Transition assignment: -> t { t > Sα && t < Sω ? 1 : 0 },
                       domain: Timer,
                       codomain: S_phase
Cdc20A_ϝ = Transition assignment: -> t { t < Cdc20Aω || t > Cdc20Aα ? 1 : 0 },
                      domain: Timer,
                      codomain: Cdc20A

# Include the A transitions in the CELL_CYCLE net.
CELL_CYCLE << A_phase_ϝ << S_phase_ϝ << Cdc20A_ϝ

def CELL_CYCLE.default_simulation
  simulation time: 0..36.h.in( :s ),
             step: 1.min.in( :s ),
             sampling: 20.min.in( :s )
end
