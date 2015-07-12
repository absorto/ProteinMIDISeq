# coding: utf8

import pyglet
from loop import *
from Bio import SeqIO
import scales

aboseq = SeqIO.read('doc/Coding.txt','fasta').seq.tostring()
abo_loop = loop_from_dna(aboseq)

# # create a loop object
l = Loop( abo_loop,
          scales.c_major,
          bpm = 60.0,
          x = 30, y=30,
          height=200, width=1200,
          midi_chan=14) # tubular bells

# create a loop object
ll = Loop( abo_loop,
           scales.c_major,
           bpm = 30.0,
           x = 30, y=500,
           height=200, width=1200,
           midi_port=u'TiMidity port 1', midi_chan=88) # new age

# main window
window = pyglet.window.Window(height=800, width=1280)
pyglet.clock.schedule_interval(l.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(l.play_at_head, l.bpm)

pyglet.clock.schedule_interval(ll.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(ll.play_at_head, ll.bpm)

pyglet.app.run()

