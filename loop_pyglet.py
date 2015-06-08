# coding: utf8

import pyglet


from loop import *

import scales




b_loop = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
]

a_loop = dna2loop("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# create a loop object
l = PianoRoll( a_loop, scale=scales.persian,
               x = 0, y=400,
               height=100, width=1140, midi_port=u'Hydrogen Midi-In', bpm = 120.0)

# create a loop object
ll = PianoRoll( b_loop, scale=scales.persian,
                x = 0, y=0,
                height=100, width=1140, midi_port=u'ZynAddSubFX', bpm = 120.0)

l.render_pianoroll_svg()




# main window
window = pyglet.window.Window(height=800, width=1280)
pyglet.clock.schedule_interval(l.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(l.play_at_head, l.bpm)


pyglet.clock.schedule_interval(ll.update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(ll.play_at_head, ll.bpm)


pyglet.app.run()

