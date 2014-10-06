# encoding: utf-8

# ==============================================================================
# MODEL PARAMETERS AND OTHER COMMON ASSETS OF VIRGINIA GENERIC CELL CYCLE MODEL
# PUBLICATION ID: Csikasznagy2006agm
# ==============================================================================

require 'sy'; require 'y_nelson' and include YNelson

# Golbeter-Koshland function used in Csikasznagy2006agm.
B = -> a1, a2, a3, a4 do a2 - a1 + a3 * a2 + a4 * a1 end
GK = -> a1, a2, a3, a4 do
  b = B.( a1, a2, a3, a4 )
  2 * a4 * a1 / ( b + ( b**2 - 4 * ( a2 - a1 ) * a4 * a1 )**0.5 )
end

# Constants as per Csikasznagy2006agm
CELL_MASS_DOUBLING_TIME = { 1 => 24.h, 2 => 14.h }[ CASE ].in :s
CELL_GROWTH_RATE = Math.log( 2 ) / CELL_MASS_DOUBLING_TIME
CycB_DIVISION_THRESHOLD = 0.3

# Following are the data from Czikasznagy2006agm supplementary materials.
# Abbreviations meanings:
# 
#    BY -- budding yeast
#    MA -- mammalian,
#    FY -- fission yeast,
#    G2 -- G2 module,
#    XE -- xenopus embryo
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

# We choose the mammalian parameter set.
PARAMETER_SET = :MA

# Singleton method that chooses a given set of constants. Constants starting
# with "K" are considered rate constants in min⁻¹. All values are converted to
# floats (Using #to_f method). Using eval, the keys are defined as constants
# in the top namespace.
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

# Choose the parameter set.
DATA.choose PARAMETER_SET
