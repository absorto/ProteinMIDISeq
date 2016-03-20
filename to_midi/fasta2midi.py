from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

from itertools import permutations

import mingus.core.scales as scales
from mingus.containers import Note, Track, Bar
from mingus.containers.instrument import Instrument
from mingus.midi.midi_file_out import write_Track

import argparse

from pprint import pprint

ph = { 'A' : 'non polar',
       'I' : 'non polar',
       'F' : 'non polar',
       'L' : 'non polar',
       'M' : 'non polar',
       'P' : 'non polar',
       'V' : 'non polar',
       'W' : 'non polar',
       
       'D' : 'acidic',
       'E' : 'acidic',

       'H' : 'basic',
       'K' : 'basic',
       'R' : 'basic',
       
       'G' : 'polar',
       'C' : 'polar',
       'N' : 'polar',
       'Q' : 'polar',
       'S' : 'polar',
       'T' : 'polar',
       'Y' : 'polar',
       '*' : 'stop'
}


durations = {transition:None for transition in permutations(['polar', 'non polar', 'acidic', 'basic', 'stop'], 2)}

durations = {
    ('acidic', 'acidic'): 1,    
    ('acidic', 'basic'): 1,
    ('basic', 'acidic'): 1,    
    ('stop', 'non polar'): 1,
    ('stop', 'polar'): 1,
    ('stop', 'stop'): 1,
    ('stop', 'acidic'): 1,
    ('stop', 'basic'): 1,
    
    ('non polar', 'polar'): 16,
    ('polar', 'non polar'): 16,
    ('non polar', 'non polar'): 8,
    ('polar', 'polar'): 8,             

    
    ('acidic', 'non polar'): 8,
    ('basic', 'non polar'): 8,
    ('acidic', 'polar'): 4,
    ('polar', 'acidic'): 4,    
    ('polar', 'basic'): 4,
    ('basic', 'basic'): 4,
    ('basic', 'polar'): 4,
    ('non polar', 'acidic'): 4,
    ('non polar', 'basic'): 4,

    
    ('acidic', 'stop'): 2,
    ('basic', 'stop'): 2,
    ('non polar', 'stop'): 2,
    ('polar', 'stop'): 2,
}

count = {('acidic', 'basic'): 0,
             ('acidic', 'non polar'): 0,
             ('acidic', 'polar'): 0,
             ('acidic', 'stop'): 0,
             ('basic', 'acidic'): 0,
             ('basic', 'non polar'): 0,
             ('basic', 'polar'): 0,
             ('basic', 'stop'): 0,
             ('non polar', 'acidic'): 0,
             ('non polar', 'basic'): 0,
             ('non polar', 'polar'): 0,
             ('non polar', 'stop'): 0,
             ('polar', 'acidic'): 0,
             ('polar', 'basic'): 0,
             ('polar', 'non polar'): 0,
             ('polar', 'stop'): 0,
             ('stop', 'acidic'): 0,
             ('stop', 'basic'): 0,
             ('stop', 'non polar'): 0,
             ('stop', 'polar'): 0,
             ('acidic', 'acidic'): 0,
             ('basic', 'basic'): 0,
             ('stop', 'stop'): 0,
             ('non polar', 'non polar'): 0,
             ('polar', 'polar'): 0,
}


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

    # populate track from amino-notes in the amino-sequence!
    for i in range(len(amino_sequence)):
        amino_note = amino_sequence[i]
        duration = durations[(ph[amino_sequence[i-1]],ph[amino_note])]
        count[(ph[amino_sequence[i-1]],ph[amino_note])] += 1
        if amino_note == "*":
            b = Bar()
            b.place_rest(duration)
            track.add_bar(b)
        else:
            track.add_notes( scale[standard_aminos.index(amino_note)], duration )

    pprint(count)
    return track





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

