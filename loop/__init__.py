# coding: utf8
import pprint

import mido
import pyglet


# un objeto loop seguramente tendrá dimensiones visuales y lógicas, y siempre va pegado a un sink midi
# al llamarse tendrian que darle sus dimensiones, y el tendria que redimensionarse apropiadamente
# al loop habría que darle: (width, height), (x,y), roll, midi_sink
class Loop:


    def __init__(self, loop, scale, bpm=120, x=0, y=0, width=800, height=600, midi_port=u'TiMidity port 0', midi_chan=1, foreground=(0,0,255), background=(255,0,0)):
        self.loop       = loop
        self.width      = width
        self.height     = height
        self.scale      = scale
        self.x          = x 
        self.playhead_x = x
        self.y          = y
        self.foreground = foreground
        self.background = background
        self.bpm        = 60/bpm
        self.beat_width = float(width)/len(self.loop)
        self.beat_height = float(height)/len(self.loop[0])
        self.tt         = 0
        self.start      = True
        self.vertex_index = []
        self.vertex_coords = []
        self.vertex_colors = []
        self.quads = 0

        self.vertexes_from_loop()

        # initialize midi output
        self.midi_output = mido.open_output( midi_port )
        self.midi_output.send(mido.Message('program_change', program=midi_chan))


    def render_playhead(self):
        if self.playhead_x <= self.width + self.x:
            pyglet.graphics.draw(2, pyglet.gl.GL_LINES,
                                 ('v2f', (self.playhead_x, self.y, self.playhead_x, self.y+self.height)))


    def render_beat(self, x, y):
        return (x, y,
                x + self.beat_width, y,
                x + self.beat_width, y + self.beat_height,
                x, y+self.beat_height)


    def vertexes_from_loop(self):
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                self.vertex_coords += self.render_beat(self.x + (x*self.beat_width),
                                                         self.y + (y*self.beat_height))

                self.vertex_index += [self.quads, self.quads+1, self.quads+2,
                                      self.quads, self.quads+2, self.quads+3]
                self.quads+=4

                if self.loop[x][y]:
                    self.vertex_colors += self.foreground * 4
                else:
                    self.vertex_colors += self.background * 4

        
    def blit_pianoroll(self):
        # here do the img
        pass
    

    def render_pianoroll(self):
        pyglet.graphics.draw_indexed(self.quads, pyglet.gl.GL_TRIANGLES,
                                     self.vertex_index,
                                     ('v2f', self.vertex_coords),
                                     ('c3B', self.vertex_colors) )




    def midi_messages(self):

        for i in range(0,len(self.loop[self.tt])):
            if self.loop[self.tt][i]:
                if not self.loop[self.tt-1][i]:
                    self.midi_output.send(mido.Message('note_on', note=self.scale[i], velocity=64))
            else:
                if self.loop[self.tt-1][i]:
                    self.midi_output.send(mido.Message('note_off', note=self.scale[i]))


    def update_playhead_display(self, dt):

        if self.start:
            self.playhead_x = self.x + 1
            self.start = False
        else:
            self.playhead_x += dt * (self.beat_width/self.bpm)

        self.render_pianoroll()
        self.render_playhead()

    def play_at_head(self, dt):
        if self.tt == len(self.loop)-1:
            self.tt = 0
            self.start = True
        else:
            self.tt+=1
        self.midi_messages()






from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


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
