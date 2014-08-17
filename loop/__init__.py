# coding: utf8


import mido
import pyglet


# un objeto loop seguramente tendrá dimensiones visuales y lógicas, y siempre va pegado a un sink midi
# al llamarse tendrian que darle sus dimensiones, y el tendria que redimensionarse apropiadamente
# al loop habría que darle: (width, height), (x,y), roll, midi_sink
class Loop:


    tt         = 0
    start      = True

    def __init__(self, loop, scale, bpm=120, x=0, y=0, width=800, height=600, midi_port=u'TiMidity port 0', midi_chan=1):
        self.loop       = loop
        self.beat_width = 30
        self.width      = width
        self.height     = height
        self.scale      = scale
        self.x          = x 
        self.playhead_x = x
        self.y          = y
        self.bpm        = 60/bpm
        self.beat_width = float(width)/len(self.loop)
        self.beat_height = float(height)/len(self.loop[0])

        print self.beat_width, self.beat_height

        # initialize midi output
        self.midi_output = mido.open_output( midi_port )
        self.midi_output.send(mido.Message('program_change', program=midi_chan))

        # load sprites
        self.on = pyglet.image.load('resources/on.png') 
        self.on_sprite = pyglet.sprite.Sprite( self.on )

        self.off = pyglet.image.load('resources/off.png') 
        self.off_sprite = pyglet.sprite.Sprite( self.off )



    def render_playhead(self):
        pyglet.graphics.draw(2, pyglet.gl.GL_LINES,
                             ('v2f', (self.playhead_x, self.y, self.playhead_x, self.y+self.height)))


    def render_pianoroll(self):
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                if self.loop[x][y]:
                    self.on.blit(self.x + (x*self.beat_width), self.y + (y*self.beat_height))
                else:
                    self.off.blit(self.x + (x*self.beat_width), self.y + (y*self.beat_height))


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
            self.playhead_x += dt * 1/self.bpm * self.beat_width

        self.render_pianoroll()
        self.render_playhead()

    def play_at_head(self, dt):
        if self.tt == len(self.loop)-1:
            self.tt = 0
            self.start = True
        else:
            self.tt+=1
        self.midi_messages()







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
