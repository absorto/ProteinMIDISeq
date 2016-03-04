# coding: utf8

import mido
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


class Loop:
    def __init__(self, loop, scale, midi_output):
        self.loop       = loop
        self.scale      = scale
        self.tt         = 0
        self.start      = True
        self.midi_output = midi_output


    def midi_messages(self):

        for i in range(0,len(self.loop[self.tt])):
            if self.loop[self.tt][i]:
                if not self.loop[self.tt-1][i]:
                    self.midi_output.send(mido.Message('note_on', note=self.scale[i], velocity=120))
            else:
                if self.loop[self.tt-1][i]:
                    self.midi_output.send(mido.Message('note_off', note=self.scale[i]))


    def print_head(self):
        print( self.loop[self.tt])
        
    def play_at_head(self):
        if self.tt == len(self.loop)-1:
            self.tt = 0
        else:
            self.tt+=1
        self.midi_messages()







def loop_from_dna(sequence):
    loop = []
    # amino list from sequence
    amino_list = Seq(sequence, IUPAC.unambiguous_dna).translate()
    
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

    for amino in amino_list:
        column = [0 for n in range(0,21)]
        if amino != '*':
            column[standard_table.back_table.keys().index(amino)] = amino

        loop.append(column)
    return loop
