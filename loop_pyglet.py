# coding: utf8
import curses
from time import sleep
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

import pyglet

from loop import *


# C3 = 48
scale = range(38,98,2)
name2midi = NoteConverter()
c_major_names = ['A', 'B','C0', 'D0', 'E0', 'F0', 'G0', 'A0', 'B0', 'C1', 'D1', 'E1', 'F1', 'G1', 'A1', 'B1', 'C2', 'D2', 'E2', 'F2', 'G2', 'A2', 'B2', ]
persian = [64+(1*0), 64+(1*1), 64+(1* 4), 64+(1*  5), 64+(1* 6), 64+(1* 8), 64+(1* 11), 64+(2*0), 64+(2*1), 64+(2* 4), 64+(2*  5), 64+(2* 6), 64+(2* 8), 64+(2* 11), 64+(3*0), 64+(3*1), 64+(3* 4), 64+(3*  5), 64+(3* 6), 64+(3* 8), 64+(3* 11)]
scale = [name2midi(n) for n in c_major_names]

a_loop = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
]

# a_loop = loop_from_dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")





# create a loop object

l = Loop( a_loop, scale, height=200 )

window = pyglet.window.Window(height=800, width=1280)

pyglet.clock.schedule_interval(l.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(l.play_at_head, 60/80.0) # update at 60Hz
pyglet.app.run()

