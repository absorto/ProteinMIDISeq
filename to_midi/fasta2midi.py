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

def track_from_dna(sequence, scale):

    # create standard table of aminos
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    standard_aminos = list(standard_table.back_table.keys())
    standard_aminos.remove(None)
    standard_aminos.sort()


    # initialize midi track
    track = Track(Instrument())


    # translate sequence argument to amino sequence
    amino_sequence = Seq(sequence, IUPAC.unambiguous_dna).translate()

    for amino_note in amino_sequence:
        if amino_note == "*":
            pass
        else:
            track.add_notes( scale[standard_aminos.index(amino_note)] )

    return track




    

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


parser = argparse.ArgumentParser(description='Convert multirecord fasta file into many MIDI files')
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True, help='a fasta file')
parser.add_argument('--bpm', default=120, help='Beats Per Minute')
args = parser.parse_args()



scale = []
for octave in range(3,6):
    for note in scales.Major("C").ascending():
        scale.append( Note(note, octave) )



for record in SeqIO.parse(args.fasta, 'fasta'):
    write_Track("%s.midi" % record.id.split("|")[1],
                track_from_dna(str(record.seq), scale),
                bpm=args.bpm)

