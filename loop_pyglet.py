# coding: utf8

import pyglet
from loop import Loop, loop_from_dna
from Bio import SeqIO
import mingus.core.scales as scales
from mingus.containers import Note



aboseq = SeqIO.read('doc/Coding.txt','fasta').seq.tostring()
abo_loop = loop_from_dna(aboseq)



A = []
for octave in range(3,6):
    for tone in scales.Aeolian("A").ascending():
        A.append( int(Note(tone, octave) ) )
A.sort()
        
a_loop = Loop( abo_loop,
               A,
               bpm = 60.0,
               x = 30, y=30,
               height=200, width=1200,
               midi_port=u'TiMidity port 0', midi_chan=88) # new age


C = []
for octave in range(1,4):
    for tone in scales.Major("C").ascending():
        C.append( int(Note(tone, octave) ) )
C.sort()
        
c_loop = Loop( abo_loop,
               C,
               bpm = 240.0,
               x = 30, y=300,
               height=200, width=1200,
               midi_port=u'TiMidity port 1', midi_chan=38) # synth bass


Em = []
for octave in range(4,7):
    for tone in scales.NaturalMinor("E").ascending():
        Em.append( int(Note(tone, octave) ) )
Em.sort()
        
e_loop = Loop( abo_loop,
               Em,
               bpm = 120.0,
               x = 30, y=530,
               height=200, width=1200,
               midi_port=u'TiMidity port 2', midi_chan=14) # tubular bells





# main window
window = pyglet.window.Window(height=800, width=1280)

pyglet.clock.schedule_interval(a_loop.update_playhead_display, 1/30.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(a_loop.play_at_head, a_loop.bpm)

pyglet.clock.schedule_interval(c_loop.update_playhead_display, 1/30.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(c_loop.play_at_head, c_loop.bpm)

pyglet.clock.schedule_interval(e_loop.update_playhead_display, 1/30.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(e_loop.play_at_head, e_loop.bpm)

pyglet.app.run()

