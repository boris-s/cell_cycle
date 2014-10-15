# encoding: utf-8

# Mammalian cell cycle model. 2 versions available at the moment,
# 'cell_cycle/simple' and 'cell_cycle/virginia_tech'. This is the
# common 'cell_cycle.rb' file.

# Require dependencies.
require 'sy'
require 'y_nelson' and include YNelson

# Define the empty CELL_CYCLE net (its contents is to be added by
# the particular implemented cell cycle types).
CELL_CYCLE = Net()

# Cell cycle version.
require_relative "cell_cycle/version"
