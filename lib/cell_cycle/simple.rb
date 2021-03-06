# encoding: utf-8

# A simplistic cell cycle. It has a single input (Timer place) and 3 outputs
# (A_phase, S_phase, Cdc20A). A_phase is the phase when the cell cycle enzyme
# machinery is synthesized. S_phase has the standard meaning: DNA synthesis
# phase. Cdc20A represents the APC (Anaphase Promoting Complex). When present,
# it degrades the cell cycle enzyme machinery.

require_relative '../cell_cycle'

# Constants that control the cell cycle settings.
Cell_cycle_duration = 24.h
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
Cc_duration_s = Cell_cycle_duration.in :s
S_start_s = S_phase_start.in :s
S_end_s = S_phase_end.in :s
A_start_s = A_phase_start.in :s
A_end_s = A_phase_end.in :s
Cdc20A_s = Cdc20A_start.in :s
Cdc20A_e = Cdc20A_end.in :s

# Timer place
Timer = Place m!: 0

# The clock transition
Clock = Transition stoichiometry: { Timer: 1 }, rate: 1

# Empirical places (in arbitrary units); output of the cell cycle.
A_phase = Place m!: 0
S_phase = Place m!: 0
Cdc20A = Place m!: 1

# Include them in the CELL_CYCLE net.
CELL_CYCLE << Timer << Clock << A_phase << S_phase << Cdc20A

# Assignment transitions that control the state of the places A_phase, S_phase
# and Cdc20A.
# 
A_phase_f = Transition domain: Timer, codomain: A_phase,
                       assignment: -> timer {
                         cc_progress = timer % Cc_duration_s
                         if cc_progress > A_start_s and cc_progress < A_end_s
                           1
                         else
                           0
                         end
                       }

S_phase_f = Transition domain: Timer, codomain: A_phase,
                       assignment: -> timer {
                         cc_progress = timer % Cc_duration_s
                         if cc_progress > S_start_s and cc_progress < S_end_s
                           1
                         else
                           0
                         end
                       }
Cdc20A_f = Transition domain: Timer, codomain: Cdc20A,
                      assignment: -> timer {
                        cc_progress = timer % Cc_duration_s
                        if cc_progress < Cdc20A_e or cc_progress > Cdc20A_s
                          1
                        else
                          0
                        end
                      }

# Include the A transitions in the CELL_CYCLE net.
CELL_CYCLE << A_phase_f << S_phase_f << Cdc20A_f

def CELL_CYCLE.default_simulation
  simulation time: 0..48.h.in( :s ),
             step: 1.min.in( :s ),
             sampling: 20.min.in( :s )
end
