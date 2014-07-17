# encoding: utf-8

# Generic cell cycle published in Csikasznagy2006agm, alternative syntax of
# the transitions.

# == TRANSITIONS ================================================================

# This creates a timed stoichiometric (TS) transition representing cell growth.
# Cell growth changes Mass (domain), and its stoichiometry is { Mass: 1 }, that
# is, Mass simply increases at the rate indicated by the transition's function.
# The function (rate) is given by the formula m * CELL_GROWTH_RATE, where m is
# the current mass of the cell, and CELL_GROWTH_RATE a constant given by the model
# authors.
# 
Cell_growth = Transition domain: Mass,
                         stoichiometry: { Mass: 1 },
                         rate: -> m { m * CELL_GROWTH_RATE }

# This creates an assignment (A) transition representing cytokinesis. It changes
# 2 places (codomain): Mass, and Ck_license. Mass is the cell mass, Ck_license is
# a special place that was added to the Petri net to prevent Cytokinesis transition
# from accidentally firing twice in a row in the same cycle. Firing depends on,
# in order, Mass, Ck_license, and ActCycB (domain). The function thus takes 3
# arguments (mass, license, b). If ActCycB is under the threshold specified by
# the model authors (CycB_DIVISION_THRESHOLD), and Ck_license is cocked
# (equal to 1), then the mass is halved and the license is consumed (set to 0).
# Otherwise, the [ mass, license ] pair is returned unchanged. The :pseudo_euler
# simulation method fires the transition once after each simulation step.
# 
Cytokinesis = Transition codomain: [ Mass, Ck_license ],
                         domain: [ Mass, Ck_license, ActCycB ],
                         assignment: -> mass, license, b do
                           if license == 1 and b < CycB_DIVISION_THRESHOLD then
                             [ mass / 2, 0 ] # mass is halved, license is set to 0
                           else
                             [ mass, license ] # nothing happens
                           end
                         end

# An assignment transition that controls the Ck_license cocking (codomain), and
# whose firin depends on Ck_license and the level of ActCycB (domain). The
# assignment function cocks the license (sets it to 1) if the level of ActCycB
# is above CycB_DIVISION_THRESHOLD plus 10% margin. Otherwise, the license is
# unchanged. Again, the :pseudo_euler simulation method fires the transition once
# after each simulation step.
# 
License_cocking = Transition codomain: Ck_license,
                             domain: [ Ck_license, ActCycB ],
                             assignment: -> license, b do
                               if b > CycB_DIVISION_THRESHOLD * 1.1 then
                                 1
                               else
                                 license
                               end
                             end

# === Module 1

# Cdc20T synthesis and degradation functions defined as lambda expressions
# according to the definitions in Csikasznagy2006agm.
# 
Cdc20T_synthesis = -> b { x = b ** N; ( Ks20p + Ks20pp * x ) / ( J20 ** N + x ) }
Cdc20T_degradation = -> cdc20T { cdc20T * Kd20 }

# A timed stoichiometric transition representing the change of anaphase-promoting
# factor (Cdc20T). Its stoichiometry is thus { Cdc20T: 1 }. Its rate depends on
# the level of activated cyclin B (ActCycB) and Cdc20T). The function is modified
# to enable larger execution step with :pseudo_euler method.
# 
Cdc20T_change = Transition domain: [ ActCycB, Cdc20T ],
                           stoichiometry: { Cdc20T: 1 },
                           rate: -> b, t {
                             step = world.simulation.step.to_f
                             fine_step = step / 50.0
                             orig = t
                             # Fine-step the function.
                             50.times do
                               t += ( Cdc20T_synthesis.( b ) - Cdc20T_degradation.( t ) ) * fine_step
                             end
                             ( t - orig ) / step # get the positive change rate
                           }

# # Without the modification for higher speed, there are 2 transitions
# # Cdc20T_synthesis and Cdc20T_degradation as follows:
# # 
# Cdc20T_synthesis = Transition domain: [ ActCycB, Cdc20T ],
#                               stoichiometry: { Cdc20T: 1 },
#                               rate: -> b {
#                                 x = b ** N; ( Ks20p + Ks20pp * x ) / ( J20 ** N + x )
#                               }
# Cdc20T_degradation = Transition stoichiometry: { Cdc20T: -1 },
#                                 rate: Kd20

# Cdc20 activation, inactivation and degradation functions defined as lambdas as
# per Csikasznagy2006agm.
# 
Cdc20_activation = -> t, a, apcp { x = t - a; Ka20 * apcp * x / ( Ja20 + x ) }
Cdc20A_inactivation = -> a { a * Ki20 / ( Ji20 + a ) }
Cdc20A_degradation = -> a { a * Kd20 }

Cdc20A_change = TS Cdc20T, Cdc20A, APCP, Cdc20A: 1 do |t, a, apcp|
  step = world.simulation.step.to_f
  fine_step = step / 50.0
  orig = a
  50.times do
    a = a +
      ( Cdc20_activation.( t, a, apcp ) -
        Cdc20A_inactivation.( a ) -
        Cdc20A_degradation.( a )
        ) * fine_step
  end
  ( a - orig ) / step # return the positive change rate
end

# REMARK: This section has to be band-aided for speed:
# Cdc20_activation = TS Cdc20T, Cdc20A, APCP, Cdc20A: 1 do |t, a, apcp|
#   x = t - a; Ka20 * apcp * x / ( Ja20 + x )
# end
# Cdc20A_inactivation = TS Cdc20A: -1 do |a| a * Ki20 / ( Ji20 + a ) end
# Cdc20A_degradation = TS Cdc20A: -1, rate: Kd20

# REMARK: Just like in CI, this section has to be band-aided.

APC_activation = -> b, apcp { x = 1 - apcp; KaAPC * b * x / ( JaAPC + x ) }
APC_inactivation = -> apcp { KiAPC * apcp / ( JiAPC + apcp ) }

APC_change = TS ActCycB, APCP, APCP: 1 do |b, apcp|
  step = world.simulation.step.to_f
  fine_step = step / 50.0
  orig = apcp
  50.times do
    apcp += ( APC_activation.( b, apcp ) - APC_inactivation.( apcp ) ) * fine_step
  end
  ( apcp - orig ) / step
end

# APC_activation = TS ActCycB, APCP, APCP: 1 do |b, apcp|
#   x = 1 - apcp; KaAPC * b * x / ( JaAPC + x )
# APC_inactivation = TS APCP: -1 do |apcp| KiAPC * apcp / ( JiAPC + apcp ) end

# === Module 2

# REMARK: Cdh1 activation and inactivation joined into 1 and band-aided in CI,
# and same has to be done here.

Cdh1_activation = -> cdc14, cdh1 {
  x = 1 - cdh1; ( Kah1p + Kah1pp * cdc14 ) * x / ( Jah1 + x ) # orig. formula
}

Cdh1_inactivation = -> a, b, d, e, cdh1 {
  ( Kih1p + Kih1pp * a + Kih1ppp * b + Kih1pppp * e + Kih1ppppp * d ) *
    cdh1 / ( Jih1 + cdh1 )
}

Cdh1_change =
  TS ActCycA, ActCycB, CycD, ActCycE, Cdc14, Cdh1, Cdh1: 1 do
  |a, b, d, e, cdc14, cdh1|
  step = world.simulation.step.to_f
  fine_step = step / 500
  orig = cdh1
  500.times do # fine-stepped formula
    cdh1 += ( Cdh1_activation.( cdc14, cdh1 ) -
              Cdh1_inactivation.( a, b, d, e, cdh1 ) ) * fine_step
  end
  ( cdh1 - orig ) / step
end

# Cdh1_activation = TS Cdc14, Cdh1, Cdh1: 1 do |cdc14, cdh1|
#   x = 1 - cdh1; ( Kah1p + Kah1pp * cdc14 ) * x / ( Jah1 + x ) # orig. formula
# end
# Cdh1_inactivation = TS ActCycA, ActCycB, CycD, ActCycE, Cdh1, Cdh1: -1 do
#   |a, b, d, e, cdh1|
#   ( Kih1p + Kih1pp * a + Kih1ppp * b + Kih1pppp * e + Kih1ppppp * d ) *
#     cdh1 / ( Jih1 + cdh1 )
# end

# === Module 4

CycB_synthesis = TS Vsb, Mass, CycB: 1 do |v, m| v * m end
CycB_degradation = TS Vdb, CycB, CycB: -1 do |v, b|
  fod { v }.fine_step( world.simulation.step, 50 ).( b )
  # v * b
end
ActCycB_synthesis = TS Vsb, Mass, ActCycB: 1 do |v, m| v * m end

ActCycB_freeing_due_to_degradation_of_CKI =
  TS Vdi, CycB, PreMPF, ActCycB, ActCycB: 1 do |v, b, preMPF, actCycB|
  v * ( b - preMPF - actCycB )
end

ActCycB_freeing_due_dissoociation_from_CKI =
  TS CycB, PreMPF, ActCycB, ActCycB: 1 do |b, preMPF, actCycB|
  Kdib * ( b - preMPF + actCycB )
end

ActCycB_creation_by_dephosphorylation_of_CycB =
  TS V25, CycB, TriB, ActCycB, ActCycB: 1 do |v, b, triB, actCycB|
  v * ( b - triB - actCycB )
end

ActCycB_phosphorylation_by_Wee1 =
  TS Vwee, ActCycB, ActCycB: -1 do |v, b| v * b end
ActCycB_asociation_with_CKI =
  TS FreeCKI, ActCycB, ActCycB: -1 do |freeCKI, b| freeCKI * b * Kasb end
ActCycB_degradation = TS Vdb, ActCycB, ActCycB: -1 do |v, b|
  fod { v }.fine_step( world.simulation.step, 50 ).( b )
  # v * b
end

# === Module 5

MPF_phosphorylation = TS Vwee, CycB, PreMPF, PreMPF: 1 do |v, b, preMPF|
  v * ( b - preMPF )
end

PreMPF_dephosphorylation =
  TS V25, PreMPF, PreMPF: -1 do |v, preMPF| v * preMPF end
PreMPF_degradation =
  TS Vdb, PreMPF, PreMPF: -1 do |v, preMPF| v * preMPF end

# === Module 6

TriB_assembly = TS CycB, TriB, FreeCKI, TriB: 1 do |b, triB, freeCKI|
  Kasb * ( b - triB ) * freeCKI
end

TriB_dissociation =
  TS TriB: -1, rate: Kdib
TriB_decrease_due_to_CycB_degradation =
  TS Vdb, TriB, TriB: -1 do |vdb, triB| vdb * triB end
TriB_decrease_due_to_CKI_degradation =
  TS Vdi, TriB, TriB: -1 do |vdi, triB| vdi * triB end

# === Module 8

CKI_synthesis = TS Vsi, CKI: 1 do |v| v end

CKI_degradation = TS Vdi, CKI, CKI: -1 do |v, cki|
  fod { v }.fine_step( world.simulation.step, 50 ).( cki )
 # v * cki
end

# === Module 10

CycE_synthesis =
  TS TFE, Mass, CycE: 1 do |f, m| ( Ksep + Ksepp * f ) * m end

CycE_degradation = TS Vde, CycE, CycE: -1 do |v, e|
  # v * e
  fod { v }.fine_step( world.simulation.step, 50 ).( e )
end

ActCycE_synthesis =
  TS TFE, Mass, ActCycE: 1 do |f, m| ( Ksep + Ksepp * f ) * m end

ActCycE_freeing_due_to_degradation_of_CKI =
  TS Vdi, TriE, ActCycE: 1 do |v, triE| v * triE end

ActCycE_freeing_due_to_dissociation_from_CKI =
  TS TriE, ActCycE: 1 do |triE| Kdie * triE end

# REMARK: band-aided in CI not to go under 0, fine-stepped here.
ActCycE_degradation = TS Vde, FreeCKI, ActCycE, ActCycE: -1 do |v, freeCKI, e|
  # ( v + Kase * freeCKI ) * e
  fod { v + Kase * freeCKI }.fine_step( world.simulation.step, 50 ).( e )
end

# === Module 13

CycA_synthesis =
  TS TFE, Mass, CycA: 1 do |f, m| ( Ksap + Ksapp * f ) * m end

CycA_degradation = TS Vda, CycA, CycA: -1 do |v, a|
  # v * a
  fod { v }.fine_step( world.simulation.step, 50 ).( a )
end

ActCycA_synthesis =
  TS TFE, Mass, ActCycA: 1 do |f, m| ( Ksap + Ksapp * f ) * m end

ActCycA_freeing_due_to_degradation_of_CKI =
  TS Vdi, TriA, ActCycE: 1 do |v, triA| v * triA end

ActCycA_freeing_due_to_dissociation_from_CKI =
  TS TriA, ActCycA: 1 do |triA| Kdia * triA end

# band-aided
ActCycA_degradation = TS Vda, FreeCKI, ActCycA, ActCycA: -1 do |v, freeCKI, a|
  # ( v + Kasa * freeCKI ) * a
  fod { v + Kasa * freeCKI }.fine_step( world.simulation.step, 50 ).( a )
end

# === finalize

finalize # YNelson command that finishes the net from the prescribed commands.

# ==============================================================================
# THIS CONCLUDES THE MODEL. COMMANDS BELOW ARE JUST SAVED COMMANDS TO SERVE AS
# HELP WHEN PERUSING THE MODEL THROUGH THE COMMAND LINE
# ==============================================================================


















=begin

# === Whole CELL_CYCLE Nelson net ===

CELL_CYCLE = Net() << A_phase << S_phase << Cdc20A

# === Cell growh net ===

CELL_GROWTH_NET = Net() <<
  Mass << Cell_growth << CycD << transition( :CycD_ϝ ) << # cell growth
  CycB << ActCycB <<                                      # cyclin B
  Ck_license << Cytokinesis << License_cocking            # cytokinesis

set_step 5                                               # simulation settings
set_sampling 300
set_target_time 3600 * 24 * 4

new_simulation


pm                            # prints marking of the newly created simulation
simulation.step!              # peforms a single simulation step
pm                            # prints marking again

# simulation.run! upto: 3600    # we can run 1 hour of this simulation

recording.plot                # and plot the recording

# and do other things with the recording, such as resampling,
# feature extraction, saving to a file, reconstruction of a new
# simulation object at a given time etc.

=end
