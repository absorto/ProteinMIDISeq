# coding: utf8
import curses
from time import sleep
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

import pyglet


from loop import *
import scales

# a_loop = loop_from_dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

a_loop = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
]

# create a loop object
l = Loop( a_loop,
          scales.persian,
          height=200, width=300 )


# main window
window = pyglet.window.Window(height=800, width=1280)
pyglet.clock.schedule_interval(l.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(l.play_at_head, 60/80.0) # update at 60Hz
pyglet.app.run()

