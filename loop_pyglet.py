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
          bpm = 60.0,
          x = 500, y=600,
          height=200, width=240)


# main window
window = pyglet.window.Window(height=800, width=1280)
pyglet.clock.schedule_interval(l.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(l.play_at_head, l.bpm)
pyglet.app.run()

