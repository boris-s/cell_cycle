#! /usr/bin/ruby
# encoding: utf-8

require 'minitest/autorun'
require 'sy'
require 'y_nelson' and include YNelson

bp = './../../lib/'
  .tap { |s| s.singleton_class.class_exec do attr_accessor :loaded end }
  .tap { |s| s.loaded = {} }

require_once = -> path do
  bp.loaded[ path ] or
    require( bp + path ) && bp.loaded.update( path => true )
end

describe :mammalian_cell_cycle do
  before do
    require_once.( 'ttp_pathway/version' )
    require_once.( 'michaelis_menten' )
    require_once.( 'general_assumptions' )
    require_once.( 'cell_cycle/virginia' )
  end

  describe "sanity of this test itself" do
    it "should have version loaded" do
      TtpPathway::VERSION.must_be_kind_of String
    end

    it "should have michaelis_menten loaded" do
      [ Vmax, Km_reduced, Occupancy, MMi ].all? &[ :kind_of?, Proc ]
    end

    it "should have general_assumptions loaded" do
      assert defined? Cell_diameter
      assert defined? Cytoplasm_volume
      assert defined? Pieces_per_µM
    end
  end

  describe "basic elements" do
    it "should have basic interface places" do
      A_phase.must_be_kind_of Place
      S_phase.must_be_kind_of Place
      Cdc20A.must_be_kind_of Place
    end

    it "should have parametrizing constants CASE and G2_MODULE" do
      assert defined? CASE
      assert defined? G2_MODULE
    end

    it "should have Godlbeter-Koshland function defined" do
      assert B.is_a? Proc     # B function as per Csikasznagy2006agm
      assert GK.is_a? Proc    # Golbeter-Koshland function
    end

    it "should have growth and cell division related constants" do
      assert defined? CELL_MASS_DOUBLING_TIME
      assert defined? CELL_GROWTH_RATE
      assert defined? CycB_DIVISION_THRESHOLD
    end

    it "should have DATA hash holding the generic cycle parameters for " +
      "budding yeast (BY), mammalian cell (MA), fission yeast (FY), and " +
      "Xenopus embryo (XE)" do
      key, val = DATA.first
      key.must_be_kind_of Symbol
      val.must_be_kind_of Hash
      val.keys.must_equal [ :BY, :MA, :FY, :G2, :XE ]
      # Now we will make an example of a single parameter that must be present
      # both as a key in the DATA hash, and as a constant in the TOP namespace.
      DATA.keys.must_include :J20
      assert defined? J20                # Must be also defined as a constant.
      J20.must_equal DATA[ :J20 ][ :MA ] # It must be set to the mammalian parameter.
    end
  end

  describe "Module 0 -- cell growth and division" do
    describe "places" do
      it "should have Mass (cell mass) place" do Mass.must_be_kind_of Place end
      it "should have CycD (cyclin D) place" do CycD.must_be_kind_of Place end
    end

    describe "assignment transitions" do
      describe "A transition controlling CycD" do
        it "must exist" do
          assert transition( :CycD_ϝ ).A? # it must exist and be an A transition
        end

        describe "its dynamic behavior" do
          before do
            @net = Net() << Mass << CycD << transition( :CycD_ϝ )
            @sim = @net.simulation initial_marking: { Mass: 1.0 }
            @feature = @net.State.Feature.Assignment( CycD )
          end

          it "must assign to CycD a value proportional to the cell mass" do
            v1 = @feature.extract_from( @sim )
            # construct a second simulation with double mass:
            sim2 = @net.simulation initial_marking: { Mass: 2.0 }
            v2 = @feature.extract_from( sim2 )
            # with double mass, the action of CycD_ϝ should be double
            v2.must_be_within_epsilon v1 * 2
          end
        end
      end
    end

    describe "transitions" do
      describe "Cell_growth" do
        it "should be of correct type" do
          assert Cell_growth.TS?
        end

        describe "its dynamic behavior" do
          before do
            @net = Net() << Mass << Cell_growth
            @sim = @net.simulation initial_marking: { Mass: 1.0 }, step: 0.1
          end

          it "should define exponential growth" do
            # Starting mass is 1. After 10 time units, the mas will be:
            @sim.run! upto: 10
            mass_after_10 = @sim.marking( :Mass ).first
            # After 20 time units, the mass should be the square of it:
            @sim.run! upto: 20
            mass_after_20 = @sim.marking( :Mass ).first
            mass_after_20.must_be_within_epsilon( mass_after_10 ** 2 )
          end
        end
      end

      describe "Cytokinesis && License_cocking" do
        it "should be of correct type" do
          assert Cytokinesis.A?
          assert License_cocking.A?
        end

        describe "dynamic behavior" do
          before do
            @net = Net() << Mass << ActCycB << Ck_license << Cytokinesis << License_cocking
          end

          describe "low cyclin B, cytokinesis license present" do
            before do
              @s1 = @net.simulation initial_marking: { Ck_license: 1, Mass: 1.0, ActCycB: 0.01 }, step: 1.0
              @s1.step!
            end

            it "cytokinesis should happen" do
              @s1.m( Ck_license ).first.must_equal 0     # consumed
              assert @s1.m( Mass ).first < 1             # decreased
            end
          end

          describe "high cyclin B, cytokinesis license absent" do
            before do
              @s2 = @net.simulation initial_marking: { Ck_license: 0, Mass: 1.0, ActCycB: 100 }, step: 1.0
              @s2.step!
            end

            it "license should cock" do
              @s2.m( Ck_license ).first.must_equal 1     # cocked
              assert @s2.m( Mass ).first >= 1            # not decreased
            end
          end

          describe "high cyclin B, cytokinesis license present" do
            before do
              @s3 = @net.simulation initial_marking: { Ck_license: 1, Mass: 1.0, ActCycB: 100 }, step: 1.0
              @s3.step!
            end

            it "cytokinesis should not happen" do
              @s3.m( Ck_license ).first.must_equal 1     # unchanged
              assert @s3.m( Mass ).first >= 1            # not decreased
            end
          end

          describe "low cyclin B, cytokinesis license absent" do
            before do
              @s4 = @net.simulation initial_marking: { Ck_license: 0, Mass: 1.0, ActCycB: 0.1 }, step: 1.0
              @s4.step!
            end

            it "nothing should happen" do
              @s4.m( Ck_license ).first.must_equal 0     # cocked
              assert @s4.m( Mass ).first >= 1            # not decreased
            end
          end
        end
      end
    end
  end

  describe "Modules 4, 10, and 13 -- synthesis and degradation of cyclins B, E, and A" do
    # Cyclin E is active primarily at the G1-S boundary, cyclin A is active from S phase
    # to early M phase, and cyclin B is essential for mitosis.

    it "must have certain places" do
      CycB.must_be_kind_of Place      # Mitotic Cdk/cyclin complex
      ActCycB.must_be_kind_of Place   # activated CycB
      CycE.must_be_kind_of Place      # G1/S transition inducer Cdk/cyclin
      ActCycE.must_be_kind_of Place   # activated CycE
      CycA.must_be_kind_of Place      # S-phase Cdk/cyclin complex
      ActCycA.must_be_kind_of Place   # S-phase Cdk/cyclin complex
    end

    describe "assigment transitions" do
      
    end

    describe "transitions" do
      
    end
  end

  describe "Modules 1 and 2 -- regulation of the anaphase promoting complex (APC)" do
    # The APC works in conjunction with Cdc20 and Cdh1 to ubiquitinylate
    # cyclin B, thereby labeling it for degradation by proteasomes. The APC
    # must be phosphorylated by the mitotic CycB kinase before it will associate
    # readily with Cdc20, but not so with Cdh1. On the other hand, Cdh1 can be
    # inactivated by phosphorylation by cyclin-dependent kinases. Cdc14 is a
    # phosphatase that opposes Cdk by dephosphorylating and activating Cdh1.

    describe "places" do
      
    end

    describe "transitions" do
      
    end
  end

  describe "Module 8 -- synthesis and degradation of CKI (cyclin-dependent kinase inhibitor)" do
    # Degradation of CKI is promoted by phoshporylation by cyclin-dependent kinases and
    # inhibited by Cdc14 phosphatase.

    describe "places" do
      
    end
  end

  describe "Modules 6, 9 and 12 -- reversible binding of CKI to cyclin/Cdk dimers to produce " +
    "catalytically inactive trimers (stoichiometric inhibition)." do
    
  end

  describe "Module 3, 7, and 11 -- regulation of the transcription factors that drive " +
    "expression of cyclins and CKI" do
    # TFB is activated by cyclin B-dependent kinase. TFE is activated by some cyclin-dependent
    # kinases and inhibited by others. TFI is inhibited by cyclin B-dependent kinase and activated
    # by Cdc14 phosphatase.
  end

  describe "Module 5 -- regulation of cyclin B-dependent kinase by tyrosine phosphorylation " +
    "and dephosphorylation (by Wee1 kinase and Cdc25 phosphatase, respectively)." do
    # The tyrosine-phosphorylated form is less active than the unphosphorylated form.
    # Cyclin B-dependent kinase phosphorylates both Wee1 (inactivating it) and Cdc25
    # (activating it), and these phosphorylations are reversed by Cdc14 phosphatase.
  end

  # The model is replete with positive feedback loops (CycB activates TFB, which drives synthesis
  # of CycB; CKI inhibits CycB, which inhibits Cdh1), and negative feedback loops (CycB activates
  # APC, which activates Cdc20, which degrades CycB; CycB activates Cdc20, which activates Cdc14,
  # which opposes CycB; TFE drives synthesis of CycA, which inhibits TFE). These complex, inter-
  # woven feedback loops create the interesting dynamical properties of the control system, which
  # account for the characteristic features of cell cycle regulation, as we intend to show.
end
