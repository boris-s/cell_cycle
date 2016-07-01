# encoding: utf-8

# *******************************************************************************
# YPetri DSL for pathway modelling, DSL source code of the mammalian cell cycle
# model from Virginia Tech made into 1 file for reference purposes.
# *******************************************************************************
# 
# 
# This file contains source code that encodes the mammalian version of the cell
# cycle model published by the team from Virginia Polytechnic Institute:
# 
# Attila Csikász-Nagy, Dorjsuren Battogtokh, Katherine C Chen,
# Béla Novák, John J Tyson, Analysis of a generic model of eukaryotic cell-cycle
# regulation., Biophys J, 2006
#
# In their paper, the Virginia tech team present a generic cell cycle model,
# where generic means that the model is applicable to a wide range of eukaryotic
# species rather than to a single species. To support their claim, the authors
# provide 4 parameter sets: Budding yeast (Saccharomyces cerevisiae, BY), fission
# yeast (Saccharomyces pombe, BY), mammalian cell (MA) and oocyte of Xenopus
# laevis (XE). From these, we choose ammalian dataset (MA). Initial state for
# the simulation was obtained as described in the original paper.
#
# Note that this source code loads YPetri gem via YNelson gem (require
# 'y_nelson'). Moreover, it uses our SY gem (require 'sy').
#
# *******************************************************************************
# For this code to work, you may have to install YNelson and SY gems by
# 'gem install y_nelson' and 'gem install sy'!!!
# *******************************************************************************

# Require the dependencies.
# 
require 'sy'
require 'y_nelson' and include YNelson
# The above line also loads YPetri

# The constant that contains the whole cell cycle net (YPetri::Net class).
# 
CELL_CYCLE = Net()

# "module CellCycle" syntax is the conventional method of specifying version
# in Ruby libraries.
# 
module CellCycle
  VERSION = "0.1.1"
end


# == TRICKS =====================================================================
# 
# The generic model from Virginia Tech reproduced here involves certain drawbacks
# (it seems fairly stiff to the author of this reproduction). These drawbacks
# necessitated defining of certain tricks, otherwise the simulation would be too
# slow. There are also other minor compatibility problems between the published
# model (which is tuned for Jigcell) and YPetri / YNelson.

# The following definitions are used to handle Float::NAN.
# 
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

# A funcion that defines first order degradation, frequently encountered
# in the model.
# 
def fod &nullary_block
  FirstOrderDegradation.new &nullary_block
end


# == ASSETS SPECIFIC TO BIOLOGICAL CONCERNS =====================================

# Golbeter-Koshland function used in Csikasznagy2006agm.
B = -> a1, a2, a3, a4 do a2 - a1 + a3 * a2 + a4 * a1 end
GK = -> a1, a2, a3, a4 do
  b = B.( a1, a2, a3, a4 )
  2 * a4 * a1 / ( b + ( b**2 - 4 * ( a2 - a1 ) * a4 * a1 )**0.5 )
end

# Csikasznagy2006agm distinguishes 2 cases and has optional G2 module.
CASE = 1
G2_MODULE = true

# Constants as per Csikasznagy2006agm
CELL_MASS_DOUBLING_TIME = { 1 => 24.h, 2 => 14.h }[ CASE ].in :s
CELL_GROWTH_RATE = Math.log( 2 ) / CELL_MASS_DOUBLING_TIME
CycB_DIVISION_THRESHOLD = 0.3

# ==============================================================================
# THE DATASETS FROM THE ORIGINAL PUBLICATION
# ==============================================================================
#
# MA -- mammalian,
# FY -- fission yeast,
# G2 -- G2 module,
# XE -- xenopus embryo
#
DATA = {
  J20:       { BY: 100,    MA: 100,      FY: 0.05,  G2: nil,  XE: nil   },
  Ja20:      { BY: 10,     MA: 0.005,    FY: 0.001, G2: nil,  XE: 0.1   },
  Ja25:      { BY: nil,    MA: nil,      FY: 0.01,  G2: 0.1,  XE: 0.1   },
  JaAPC:     { BY: 0.1,    MA: 0.01,     FY: 0.001, G2: nil,  XE: 0.01  },
  Jafb:      { BY: 0.1,    MA: 0.1,      FY: nil,   G2: nil,  XE: nil   },
  Jafi:      { BY: 10,     MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Jah1:      { BY: 0.03,   MA: 0.01,     FY: 0.01,  G2: nil,  XE: nil   },
  Jatf:      { BY: 0.01,   MA: 0.01,     FY: 0.01,  G2: nil,  XE: nil   },
  Jawee:     { BY: nil,    MA: nil,      FY: 0.01,  G2: 0.05, XE: 0.3   },
  Ji20:      { BY: 10,     MA: 0.005,    FY: 0.001, G2: nil,  XE: 0.1   },
  Ji25:      { BY: nil,    MA: nil,      FY: 0.01,  G2: 0.1,  XE: 0.1   },
  JiAPC:     { BY: 0.1,    MA: 0.01,     FY: 0.001, G2: nil,  XE: 0.01  },
  Jifb:      { BY: 0.1,    MA: 0.1,      FY: nil,   G2: nil,  XE: nil   },
  Jifi:      { BY: 10,     MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Jih1:      { BY: 0.03,   MA: 0.01,     FY: 0.01,  G2: nil,  XE: nil   },
  Jitf:      { BY: 0.01,   MA: 0.01,     FY: 0.01,  G2: nil,  XE: nil   },
  Jiwee:     { BY: nil,    MA: nil,      FY: 0.01,  G2: 0.05, XE: 0.3   },
  J14di:     { BY: 0.0833, MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  K25p:      { BY: nil,    MA: nil,      FY: 0.001, G2: 0.05, XE: 0.1   },
  K25pp:     { BY: nil,    MA: nil,      FY: 1,     G2: 05,   XE: 1.9   },

  Ka20:      { BY: 1,      MA: 0.0833,   FY: 0.2,   G2: nil,  XE: 0.1   },
  Ka25:      { BY: nil,    MA: nil,      FY: 1,     G2: 1,    XE: 1     },
  KaAPC:     { BY: 0.1,    MA: 0.0117,   FY: 0.2,   G2: nil,  XE: 2     },
  Kafb:      { BY: 1,      MA: 0.167,    FY: nil,   G2: nil,  XE: nil   },
  Kafi:      { BY: 6,      MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Kah1p:     { BY: 0.01,   MA: 0.175,    FY: 5,     G2: nil,  XE: nil   },
  Kah1pp:    { BY: 0.8,    MA: 2.33,     FY: 50,    G2: nil,  XE: nil   },
  Kasa:      { BY: 50,     MA: 16.7,     FY: 500,   G2: nil,  XE: nil   },
  Kasb:      { BY: 65,     MA: nil,      FY: 1000,  G2: nil,  XE: nil   },
  Kase:      { BY: nil,    MA: 16.7,     FY: nil,   G2: nil,  XE: nil   },
  Katfp:     { BY: nil,    MA: 0.0,      FY: 1.5,   G2: nil,  XE: nil   },
  Katfpp:    { BY: 0.76,   MA: 0.05,     FY: nil,   G2: nil,  XE: nil   },
  Katfppp:   { BY: 0.76,   MA: 0.0833,   FY: nil,   G2: nil,  XE: nil   },
  Katfpppp:  { BY: 3.8,    MA: 0.055,    FY: nil,   G2: nil,  XE: nil   },
  Kaweep:    { BY: nil,    MA: nil,      FY: 0.25,  G2: 0.3,  XE: 0.1   },
  Kaweepp:   { BY: nil,    MA: nil,      FY: 0.25,  G2: nil,  XE: nil   },
  Kd20:      { BY: 0.05,   MA: 0.025,    FY: 0.1,   G2: nil,  XE: nil   },
  Kdap:      { BY: 0.01,   MA: 0.000333, FY: 0.01,  G2: nil,  XE: nil   },
  Kdapp:     { BY: 0.16,   MA: 0.333,    FY: 2,     G2: nil,  XE: nil   },
  Kdappp:    { BY: nil,    MA: nil,      FY: 0.02,  G2: nil,  XE: nil   },

  Kdbp:      { BY: 0.003,  MA: 0.000833, FY: 0.02,  G2: nil,  XE: 0.015 },
  Kdbpp:     { BY: 0.4,    MA: 0.333,    FY: 0.75,  G2: nil,  XE: nil   },
  Kdbppp:    { BY: 0.15,   MA: 0.0167,   FY: 1.5,   G2: nil,  XE: nil   },
  Kdep:      { BY: 0.12,   MA: 0.00167,  FY: nil,   G2: nil,  XE: nil   },
  Kdepp:     { BY: nil,    MA: 0.0167,   FY: nil,   G2: nil,  XE: nil   },
  Kdeppp:    { BY: nil,    MA: 0.167,    FY: nil,   G2: nil,  XE: nil   },
  Kdepppp:   { BY: nil,    MA: 0.167,    FY: nil,   G2: nil,  XE: nil   },
  Kdia:      { BY: 0.06,   MA: 0.167,    FY: 1,     G2: nil,  XE: nil   },
  Kdib:      { BY: 0.05,   MA: nil,      FY: 1,     G2: nil,  XE: nil   },
  Kdie:      { BY: nil,    MA: 0.167,    FY: nil,   G2: nil,  XE: nil   },
  Kdip:      { BY: 0.02,   MA: 0.167,    FY: 0.1,   G2: nil,  XE: nil   },
  Kdipp:     { BY: 0.2,    MA: 0.833,    FY: 2,     G2: nil,  XE: nil   },
  Kdippp:    { BY: 0.9,    MA: 1.67,     FY: 100,   G2: nil,  XE: nil   },
  Kdipppp:   { BY: 0.12,   MA: 0.833,    FY: nil,   G2: nil,  XE: nil   },
  Kdippppp:  { BY: 0.66,   MA: nil,      FY: 1,     G2: nil,  XE: nil   },
  Ki20:      { BY: 0.05,   MA: 0.0417,   FY: 0.05,  G2: nil,  XE: 0.095 },
  Ki25p:     { BY: nil,    MA: nil,      FY: 0.25,  G2: 0.3,  XE: 0.125 },
  Ki25pp:    { BY: nil,    MA: nil,      FY: 0.25,  G2: nil,  XE: nil   },
  KiAPC:     { BY: 0.15,   MA: 0.03,     FY: 0.08,  G2: nil,  XE: 0.15  },
  Kifb:      { BY: 0.15,   MA: 0.0167,   FY: nil,   G2: nil,  XE: nil   },

  Kifip:     { BY: 0.008,  MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Kifipp:    { BY: 0.05,   MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Kih1p:     { BY: 0.001,  MA: nil,      FY: 1,     G2: nil,  XE: nil   },
  Kih1pp:    { BY: 0.64,   MA: 0.2,      FY: 40,    G2: nil,  XE: nil   },
  Kih1ppp:   { BY: 0.1,    MA: 0.667,    FY: 40,    G2: nil,  XE: nil   },
  Kih1pppp:  { BY: 0.032,  MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Kih1ppppp: { BY: 0.01,   MA: nil,      FY: 40,    G2: nil,  XE: nil   },
  Kitfp:     { BY: 0.6,    MA: 0.0417,   FY: 1,     G2: nil,  XE: nil   },
  Kitfpp:    { BY: 8,      MA: 0.0167,   FY: nil,   G2: nil,  XE: nil   },
  Kitfppp:   { BY: nil,    MA: 0.0167,   FY: 10,    G2: nil,  XE: nil   },
  Kiwee:     { BY: nil,    MA: nil,      FY: 1,     G2: 1,    XE: 3     },
  Ks20p:     { BY: 0.001,  MA: nil,      FY: 0.005, G2: nil,  XE: 1     },
  Ks20pp:    { BY: 10,     MA: 2.5,      FY: 0.1,   G2: nil,  XE: nil   },
  Ksap:      { BY: 0.0008, MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Ksapp:     { BY: 0.005,  MA: 0.00417,  FY: 0.02,  G2: nil,  XE: nil   },
  Ksbp:      { BY: 0.004,  MA: 0.00167,  FY: 0.02,  G2: nil,  XE: 0.1   },
  Ksbpp:     { BY: 0.04,   MA: 0.005,    FY: nil,   G2: nil,  XE: nil   },
  Ksep:      { BY: nil,    MA: 0.00133,  FY: nil,   G2: nil,  XE: nil   },
  Ksepp:     { BY: 0.15,   MA: 0.05,     FY: nil,   G2: nil,  XE: nil   },
  Ksip:      { BY: 0.036,  MA: 0.333,    FY: 0.3,   G2: nil,  XE: nil   },

  Ksipp:     { BY: 0.24,   MA: nil,      FY: nil,   G2: nil,  XE: nil   },
  Kweep:     { BY: nil,    MA: nil,      FY: 0.05,  G2: 0.2,  XE: 0.1   },
  Kweepp:    { BY: nil,    MA: nil,      FY: 0.5,   G2: 2,    XE: 0.9   },
  N:         { BY: 1,      MA: 1,        FY: 4,     G2: nil,  XE: nil   },
  CycD⁰:     { BY: 0.108,  MA: 0.5,      FY: 0.05,  G2: nil,  XE: nil   }
}

# Singleton method that chooses a given set of constants. Also deals with the
# quirks of the original dataset. (Constants starting with "K" are considered
# rate constants in min⁻¹. All values are converted to floats. Using eval, the
# keys are also defined as constants.)
# 
def DATA.choose which
  keys.each { |key|
    a = [ "DATA[:#{key}][:#{which}]" ]
    # Symbols started with "K" are converted from min⁻¹ to num. values in s⁻¹
    x = if "#{key}"[0] == "K" then ".min⁻¹.in :s⁻¹" end
    # This version replaces nils with zeroes:
    eval "#{key} = %s.to_f#{x}" % a
    # This version keeps nils:
    # eval "#{key} = if %s.is_a? Numeric then %s.to_f#{x} else %s end" % ( a * 3 )
  }
end

# Mammalian parameter set is chosen for this model.
PARAMETER_SET = :MA
DATA.choose PARAMETER_SET


# ==============================================================================
# THE MODEL ITSELF
# ==============================================================================

# Empirical places (in arbitrary units); they are the output of the cell cycle.
# 
A_phase = Place m!: 0                    # in situ
S_phase = Place m!: 0                    # in situ
# TODO: Cdc20A is a part of the Virginia cycle, but not as 0/1 value, user must be adaptable
# Cdc20A = Place m!: 1                     # in situ p

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

# Cdh1_activation = TS Cdc14, Cdh1, Cdh1: 1 do |cdc14, cdh1|
#   x = 1 - cdh1; ( Kah1p + Kah1pp * cdc14 ) * x / ( Jah1 + x ) # orig. formula
# end
# 
# Cdh1_inactivation = TS ActCycA, ActCycB, CycD, ActCycE, Cdh1, Cdh1: -1 do
#   |a, b, d, e, cdh1|
#   ( Kih1p + Kih1pp * a + Kih1ppp * b + Kih1pppp * e + Kih1ppppp * d ) *
#     cdh1 / ( Jih1 + cdh1 )
# end

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

# This transition includes an ad-hoc modification, but the character of the
# original published model is preserved.
# 
ActCycA_degradation = TS Vda, FreeCKI, ActCycA, ActCycA: -1 do |v, freeCKI, a|
  # ( v + Kasa * freeCKI ) * a
  fod { v + Kasa * freeCKI }.fine_step( world.simulation.step, 50 ).( a )
end

# === Finalization

finalize # YNelson command that finalizes the net from the prescribed commands.
CELL_CYCLE.merge! net # The net defined thus far is put into CELL_CYCLE.

# This line defines a new Net object that contains places and transitions related
# to cell growth. The purpose of defining this as a separate net is eg. due to
# testing or separate use of the cell growth net.
# 
CELL_GROWTH_NET = Net()
CELL_GROWTH_NET << Mass << Cell_growth << CycD <<
  transition( :CycD_ϝ ) << CycB << ActCycB <<         # cyclin B
  Ck_license << Cytokinesis << License_cocking        # cytokinesis


# ==============================================================================
# 
# THIS CONCLUDES THE MODEL.
#
# By running the statements in this file up to this point, the model
# specification task is accomplished. The resulting net object (of YPetri::Net
# class) can be simulated using methods available in YPetri
# 
# ==============================================================================


# There are multiple ways of running the simulation and interacting with the
# results using YPetri / YNelson. One of them is eg.

set_step 5
set_target_time 96.h.in( :s ) # 96 hour simulation
set_sampling 300

# Beware! Simulation aspect of YPetri / YNelson is currently not optimized!
# Running 96 hour simulation takes a very long time!

run!

# Interaction with the results (resampling, feature extraction, saving to a file,
# reconstruction of a new simulation object at a given time etc. is possible).
