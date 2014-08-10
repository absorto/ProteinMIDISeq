import curses
from time import sleep
import mido
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

import pyglet

window = pyglet.window.Window(height=900)


class NoteConverter:
    def __init__(self, base=5):
        self.base = base
        names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
        self.notes = {}
        for note in range(128):
            name = names[note % 12]
            octave = int(note / 12) - base
            key = '%s%i' % (name, octave)
            self.notes[key] = note
            if octave == 0:
                self.notes[name] = note

    def __call__(self, note):
        return self.notes.get(str(note).upper(), note)


def loop_from_dna(sequence):
    loop = []
    # amino list from sequence
#    amino_list = list(  Seq(sequence, IUPAC.unambiguous_dna).translate(to_stop=True).tostring() )
    amino_list = list(  Seq(sequence, IUPAC.unambiguous_dna).translate().tostring() )
    
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

    for amino in amino_list:
        column = [0 for n in range(0,21)]
        if amino != '*':
            column[standard_table.back_table.keys().index(amino)] = 1

        loop.append(column)
    return loop



class loop:

    t = 0

    def __init__(self, midi_output, loop):
        self.loop = loop

        self.output = midi_output
        self.output.send(mido.Message('program_change', program=19))

        self.on = pyglet.image.load('resources/on.png') 
        self.on_sprite = pyglet.sprite.Sprite( self.on )

        self.off = pyglet.image.load('resources/off.png') 
        self.off_sprite = pyglet.sprite.Sprite( self.off )

        self.width = len(loop)*30


    def render_playhead(self):
        pyglet.graphics.draw(2, pyglet.gl.GL_LINES,
                             ('v2f', (self.t, 0, self.t, 1000)))


    def render_pianoroll(self):
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                if self.loop[x][y]:
                    self.on.blit(x*30, y*30)
                else:
                    self.off.blit(x*30, y*30)


    def midi_messages(self, t, scale):

        for i in range(0,len(self.loop[t])):
            if self.loop[t][i]:
                if not self.loop[t-1][i]:
                    self.output.send(mido.Message('note_on', note=scale[i], velocity=64))
            else:
                if self.loop[t-1][i]:
                    self.output.send(mido.Message('note_off', note=scale[i]))







# C3 = 48
scale = range(38,98,2)
name2midi = NoteConverter()
c_major_names = ['A', 'B','C0', 'D0', 'E0', 'F0', 'G0', 'A0', 'B0', 'C1', 'D1', 'E1', 'F1', 'G1', 'A1', 'B1', 'C2', 'D2', 'E2', 'F2', 'G2', 'A2', 'B2', ]
persian = [64+(1*0), 64+(1*1), 64+(1* 4), 64+(1*  5), 64+(1* 6), 64+(1* 8), 64+(1* 11), 64+(2*0), 64+(2*1), 64+(2* 4), 64+(2*  5), 64+(2* 6), 64+(2* 8), 64+(2* 11), 64+(3*0), 64+(3*1), 64+(3* 4), 64+(3*  5), 64+(3* 6), 64+(3* 8), 64+(3* 11)]
scale = [name2midi(n) for n in c_major_names]

a_loop = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 1],
]

a_loop = loop_from_dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# create a loop object
l = loop( mido.open_output( u'TiMidity port 0'), 
          a_loop)



def update(dt):

    if l.t > l.width-1:
        l.t = 0
    else:
        l.t += dt * 50

    l.render_pianoroll()
    l.render_playhead()


pyglet.clock.schedule_interval(update, 1/60.0) # update at 60Hz
pyglet.app.run()

