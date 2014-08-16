# coding: utf8
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


# un objeto loop seguramente tendrá dimensiones visuales y lógicas, y siempre va pegado a un sink midi
# al llamarse tendrian que darle sus dimensiones, y el tendria que redimensionarse apropiadamente
# al loop habría que darle: (width, height), (x,y), roll, midi_sink
class loop:

    playhead_x = 0
    tt         = 0
    start      = True

    def __init__(self, midi_output, loop):
        self.loop = loop

        self.output = midi_output
        self.output.send(mido.Message('program_change', program=1))

        self.on = pyglet.image.load('resources/on.png') 
        self.on_sprite = pyglet.sprite.Sprite( self.on )

        self.off = pyglet.image.load('resources/off.png') 
        self.off_sprite = pyglet.sprite.Sprite( self.off )

        self.width = len(loop)*30


    def render_playhead(self):
        pyglet.graphics.draw(2, pyglet.gl.GL_LINES,
                             ('v2f', (self.playhead_x, 0, self.playhead_x, 1000)))


    def render_pianoroll(self):
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                if self.loop[x][y]:
                    self.on.blit(x*30, y*30)
                else:
                    self.off.blit(x*30, y*30)


    def midi_messages(self, scale):

        for i in range(0,len(self.loop[self.tt])):
            if self.loop[self.tt][i]:
                if not self.loop[self.tt-1][i]:
                    self.output.send(mido.Message('note_on', note=scale[i], velocity=64))
            else:
                if self.loop[self.tt-1][i]:
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
    [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
]

# a_loop = loop_from_dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

# create a loop object
l = loop( mido.open_output( u'TiMidity port 0'), 
          a_loop)



def update_playhead_display(dt):
    if l.start:
        l.playhead_x = 1
        l.start = False
    else:
        l.playhead_x += dt * 1.33 * 30
    l.render_pianoroll()
    l.render_playhead()

def play_at_head(dt):
    if l.tt == len(l.loop)-1:
        l.tt = 0
        l.start = True
    else:
        l.tt+=1
    l.midi_messages(scale)

pyglet.clock.schedule_interval(update_playhead_display, 1/60.0) # Update at 60Hz = 1/60
pyglet.clock.schedule_interval(play_at_head, 60/80.0) # update at 60Hz
pyglet.app.run()

