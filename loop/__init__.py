# coding: utf8
import pprint

from uuid import uuid4 as random_uuid
import svgwrite
from sh import convert, rm

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
        self.beat_height = float(height)/(len(self.loop[0]) - 1)
        self.tt         = 0
        self.start      = True
        self.vertex_index = []
        self.vertex_coords = []
        self.vertex_colors = []
        self.quads = 0

        # initialize midi output
        self.midi_output = mido.open_output( midi_port )
        self.midi_output.send(mido.Message('program_change', program=midi_chan))

        # render SVG
        loop_hash = random_uuid()
        self.dwg = svgwrite.Drawing(filename="%s.svg" % loop_hash, size=(width, height), debug=True)
        self.render_pianoroll_svg()

        # svg to png
        convert("%s.svg" % loop_hash,
                "%s.png" % loop_hash)
        
        # load png
        self.sprite = pyglet.sprite.Sprite( pyglet.image.load("%s.png" % loop_hash,
                                                              file=open("%s.png" % loop_hash)))

        
        rm("%s.svg" % loop_hash)
        rm("%s.png" % loop_hash)

        
        self.sprite.x = self.x
        self.sprite.y = self.y
        
    def render_playhead(self):
        if self.playhead_x <= self.width + self.x:
            pyglet.graphics.draw(2, pyglet.gl.GL_LINES,
                                 ('v2f', (self.playhead_x, self.y, self.playhead_x, self.y+self.height)),
                                 ('c3B', (255,0,0, 255,0,0)))


    def render_beat(self, x, y):
        return (x, y,
                x + self.beat_width, y,
                x + self.beat_width, y + self.beat_height,
                x, y+self.beat_height)

        
    def render_pianoroll_sprite(self):
        self.sprite.draw()

    def render_pianoroll_svg(self):
        self.dwg.add(self.dwg.rect(insert=(0, 0),
                                   size  =(self.width, self.height),
                                   fill  = 'lightgoldenrodyellow', stroke_width=0))
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                if self.loop[x][y]:
                    self.dwg.add(self.dwg.rect(insert=(x*self.beat_width,
                                                       (self.height - (self.beat_height + (y*self.beat_height)))),
                                               size=(self.beat_width,
                                                     self.beat_height),
                                               fill='lightsteelblue', stroke_width=0))
        self.dwg.save()
    


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

        #self.render_pianoroll()
        self.render_pianoroll_sprite()
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
