from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


import mingus.core.scales as scales
from mingus.containers import Note, Track, Bar
from mingus.containers.instrument import Instrument
from mingus.midi.midi_file_out import write_Track


import argparse
from pprint import pprint


# class Loop:
#     def __init__(self, loop, scale):
#         self.loop       = loop
#         self.scale      = scale

#     def midi_messages(self):

#         for i in range(0,len(self.loop[self.tt])):
#             if self.loop[self.tt][i]:
#                 if not self.loop[self.tt-1][i]:
#                     self.midi_output.send(mido.Message('note_on', note=self.scale[i], velocity=120))
#             else:
#                 if self.loop[self.tt-1][i]:
#                     self.midi_output.send(mido.Message('note_off', note=self.scale[i]))


#     def print_head(self):
#         print( self.loop[self.tt])
        
#     def play_at_head(self):
#         if self.tt == len(self.loop)-1:
#             self.tt = 0
#         else:
#             self.tt+=1
#         self.midi_messages()

            




def loop_from_dna(sequence):
    loop = []
    # amino list from sequence
    amino_list = Seq(sequence, IUPAC.unambiguous_dna).translate()
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    aminos = list(standard_table.back_table.keys())
    aminos.remove(None)
    aminos.sort()
    for amino in amino_list:
        column = [0 for n in range(0,21)]
        if amino != "*":
            column[aminos.index(amino)] = amino

        loop.append(column)
    return loop



def midi_from_loop(loop, path):
    scale = []
    for octave in range(3,6):
        for note in scales.Major("C").ascending():
            scale.append( Note(note, octave) )

    track = Track(Instrument())
    b = Bar()
    b.set_meter((6,8))
    track.add_bar(b)

    for t in range(0,len(loop)):
        for g in range(len(loop[t])):
            if loop[t][g]:
                track.add_notes( scale[g] )

    write_Track(path, track, bpm=240)


parser = argparse.ArgumentParser(description='Convert multirecord fasta file into many midis')
parser.add_argument('fasta', type=argparse.FileType('r'), help='a fasta file')
args = parser.parse_args()

for record in SeqIO.parse(args.fasta, 'fasta'):
    path = "%s.midi" % record.id.split("|")[1]
    midi_from_loop(loop_from_dna(str(record.seq)), path)
    


