# encoding: utf-8

# Generic cell cycle published in Csikasznagy2006agm. Mammalian constants used.
# Initial state was obtained as described in Csikasznagy2006agm.

require 'cell_cycle'

# == TRICKS =====================================================================

# The following definitions are used to handle Float::NAN.
class Float; def nan_to_zero; nan? ? 0.0 : self end end
class Numeric; def to_f_nan_0; to_f.nan_to_zero end end

# The following definitions are used to fine-step degradation.
# 
class FirstOrderDegradation < Proc
  def fine_step( step, n )
    fine = step.to_f / n
    -> v { orig = v; n.times { v -= call * v * fine }; ( orig - v ) / step }
  end
end

def fod &nullary_block; FirstOrderDegradation.new &nullary_block end

# ==============================================================================
# THE MODEL ITSELF
# ==============================================================================

# Empirical places (in arbitrary units); they are the output of the cell cycle.
# 
A_phase = Place m!: 0                    # in situ
S_phase = Place m!: 0                    # in situ
# TODO: Cdc20A is a part of the Virginia cycle, but not as 0/1 value, user must be adaptable
# Cdc20A = Place m!: 1                     # in situ p

# Csikasznagy2006agm distinguishes 2 cases and has optional G2 module.
CASE = 1; G2_MODULE = true


# ==============================================================================
# THERE IS ANOTHER FILE WITH MODEL PARAMETERS THAT GETS LOADED HERE:
# cell_cycle/virginia/mammalian_constants.rb
# ==============================================================================

begin; require_relative 'virginia_tech/parameters'; rescue LoadError
  require './parameters'
end

# begin
#   require_relative 'mammalian_constants'
# rescue LoadError
#   begin
#     require './mammalian_constants'
#   rescue LoadError
#     require 'cell_cycle/mammalian_constants.rb'
#   end
# end

# == PLACES =====================================================================

# Cell mass
Mass = Place m!: { 1 => 1.098, 2 => 1.568 }[ CASE ]

# License that ensures that the cytokinesis occurs not more than once per cycle.
Ck_license = Place m!: 0

# === Module 1
Cdc20T = Place m!: { 1 => 2.66, 2 => 2.7 }[ CASE ] # Cdc20 total
Cdc20A = Place m!: { 1 => 0.708, 2 => 0.76 }[ CASE ]
APCP = Place m!: { 1 => 0.717, 2 => 0.78 }[ CASE ]

# === Module 2
Cdh1 = Place m!: { 1 => 0.999, 2 => 0.999 }[ CASE ]

# === Module 4
CycB = Place m!: { 1 => 0.289, 2 => 0.5 }[ CASE ] # Cdc2/CycB
ActCycB = Place m!: { 1 => 0.289, 2 => 0.22 }[ CASE ]

# === Module 5
PreMPF = Place m!: { 1 => 0, 2 => 0.29 }[ CASE ]

# === Module 6
TriB = Place m!: { 1 => 0, 2 => 0 }[ CASE ] # trimer CKI / CycB / Cdk

# === Module 8
CKI = Place m!: { 1 => 0.343, 2 => 0.26 }[ CASE ]

# === Module 10
CycE = Place m!: { 1 => 0.414, 2 => 0.73 }[ CASE ]
ActCycE = Place m!: { 1 => 0.181, 2 => 0.53 }[ CASE ]

# === Module 13
CycA = Place m!: { 1 => 0.0280, 2 => 0.062 }[ CASE ]
ActCycA = Place m!: { 1 => 0.0124, 2 => 0.045 }[ CASE ]

# == PLACES WITH ASSIGNMENT TRANSITIONS =========================================

CycD = ϝ Mass do |m| m * CycD⁰ end

# === Module 1
Cdc14 = ϝ :Cdc20A do |a| a end

# === Module 3
TFB = ϝ :ActCycB do |b| GK.( Kafb * b, Kifb, Jafb, Jifb ) end

# === Module 4
Vsb = ϝ :TFB do |tfb| Ksbp + Ksbpp * tfb end
Vdb = ϝ :Cdh1, :Cdc20A do |cdh1, a| Kdbp + Kdbpp * cdh1 + Kdbppp * a end
Cdk1P_CycB = ϝ :CycB, :ActCycB, :TriB do |b, ab, triB| b - ab - triB end
Cdk1_CycB_CKI = ϝ :CycB, :ActCycB, :PreMPF do |b, ab, f| b - ab - f end

# === Module 5 (G2 module) -- not included in mammalian cycle
Cdc25 = ϝ :ActCycB, :Cdc14 do |b, cdc14|
  # Rescue 0 makes the closure return 0 rather than raising error over Ki25p etc.
  # constants being equal to nil (not present in the mammalian cycle).
  GK.( Ka25 * b, Ki25p + Ki25pp * cdc14, Ja25, Ji25 ).to_f_nan_0 rescue 0
end
V25 = ϝ :Cdc25 do |cdc25| ( K25p + K25pp * cdc25 ).to_f_nan_0 rescue 0 end
Vwee = ϝ :Wee1 do |wee1| ( Kweep + Kweepp * wee1 ).to_f_nan_0 rescue 0 end
Wee1 = ϝ :Cdc14, :ActCycB do |cdc14, b|
  GK.( Kaweep + Kaweepp * cdc14, Kiwee * b, Jawee, Jiwee ).to_f_nan_0 rescue 0
end

# === Module 7 -- Not included in the mammalian cycle
TFI = ϝ :Cdc14, :ActCycB do |cdc14, actCycB|
  GK.( Kafi * cdc14, Kifip + Kifipp * actCycB, Jafi, Jifi ).to_f_nan_0 rescue 0
end

# === Module 12
TriA = ϝ :CycA, :ActCycA do |a, actCycA| [ a - actCycA, 0 ].max rescue 0 end

# === Module 8
Vsi = ϝ :TFI do |tfi| ( Ksip + ( Ksipp * tfi ) ).to_f_nan_0 rescue 0 end

Vdi = case PARAMETER_SET # The function used depends on the parameter set.
      when :BY then
        ϝ :ActCycA, :ActCycB, :ActCycE, :CycD, :Cdc14 do |a, b, e, d, cdc14|
          ( Kdip + Kdipp * a + Kdippp * b + Kdipppp * e + Kdippppp * d ) /
          ( 1 + cdc14 / J14di ) # <---- problem, J14di is 0
        end
      else
        ϝ :ActCycA, :ActCycB, :ActCycE, :CycD, :Cdc14 do |a, b, e, d, _|
          Kdip + Kdipp * a + Kdippp * b + Kdipppp * e + Kdippppp * d
        end
      end

# REMARK: Band-aided not to go under 0
FreeCKI = ϝ :CKI, :TriA, :TriB, :TriE do |cki, a, b, e|
  [ cki - a - b - e, 0.0 ].max
end
# FreeCKI = ϝ :CKI, :TriA, :TriB, :TriE do |cki, a, b, e| cki - a - b - e end

# === Module 9
# REMARK: Band-aided not to go under 0.
TriE = ϝ :CycE, :ActCycE do |e, actCycE|
  [ e - actCycE, 0 ].max
end

# === Module 10
Vde = ϝ :ActCycA, :ActCycB, :ActCycE do |a, b, e|
  Kdep + Kdepp * e + Kdeppp * a + Kdepppp * b
end

# === Module 11
Vatf = ϝ :ActCycA, :ActCycE, :CycD do |actCycA, actCycE, d|
  Katfp + Katfpp * actCycA + Katfppp * actCycE + Katfpppp * d
end

TFE = ϝ :Vatf, :ActCycA, :ActCycB do |v, a, b|
  GK.( v, Kitfp + Kitfpp * b + Kitfppp * a, Jatf, Jitf )
end

# === Module 13
Vda = ϝ :Cdc20A, :Cdc20T do |a, t| Kdap + Kdapp * a + Kdappp * t end

# == TRANSITIONS ================================================================

Cell_growth = TS Mass, Mass: 1 do |m| m * CELL_GROWTH_RATE end

Cytokinesis = AT Mass, Ck_license, domain: [ Mass, Ck_license, ActCycB ] do
  |mass, license, b|
  # Cytokinesis happens when ActCycB drops below the threshold
  if license == 1 and b < CycB_DIVISION_THRESHOLD then
    [ mass / 2, 0 ] # mass is halved, license consumed (set to value 0)
  else
    [ mass, license ]
  end
end

License_cocking = AT Ck_license, domain: [ Ck_license, ActCycB ] do |license, b|
  # License is cocked (value 1) if ActCycB goes above the threshold * 1.1
  if b > CycB_DIVISION_THRESHOLD * 1.1 then 1 else license end
end

# === Module 1

Cdc20T_synthesis = -> b { x = b ** N; ( Ks20p + Ks20pp * x ) / ( J20 ** N + x ) }
Cdc20T_degradation = -> cdc20T { cdc20T * Kd20 }

Cdc20T_change = TS ActCycB, Cdc20T, Cdc20T: 1 do |b, t|
  step = world.simulation.step.to_f
  fine_step = step / 50.0
  orig = t
  50.times do
    t += ( Cdc20T_synthesis.( b ) - Cdc20T_degradation.( t ) ) * fine_step
  end
  ( t - orig ) / step # return the positive change rate
end

# REMARK: This section has to be band-aided for speed:
# Cdc20T_synthesis = TS ActCycB, Cdc20T: 1 do |b|
#   x = b ** N; ( Ks20p + Ks20pp * x ) / ( J20 ** N + x )
# end
# Cdc20T_degradation = TS Cdc20T: -1, rate: Kd20

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


CELL_CYCLE.merge! net   # This is quite fragile, since it assumes that no places
                        # and/or transitions have been created before requiring
                        # the cell cycle gem.


# Net that defines cell growth, to be used eg. for testing or separate use.
CELL_GROWTH_NET = Net()
CELL_GROWTH_NET << Mass << Cell_growth << CycD <<
  transition( :CycD_ϝ ) << CycB << ActCycB <<         # cyclin B
  Ck_license << Cytokinesis << License_cocking        # cytokinesis

# cell growth
def CELL_CYCLE.default_simulation
  simulation time: 0..96.h.in( :s ),
             step: 5,
             sampling: 300
  # FIXME: This method #default_simulation doesn't work, and I don't know why,
  # while standard way of simulation works. I can just use the standard way
  # of simulation by simply calling run! command, and fix this later.
end

set_step 5
set_target_time 96.h.in( :s )
set_sampling 300

# I commented this out, I don't really want to run this lengthy simulation upon
# requiring the gem.
# run!

# Note: Possible things to do with the recording, such as resampling, feature
# extraction, saving to a file, reconstruction of a new simulation object at
# a given time etc.
