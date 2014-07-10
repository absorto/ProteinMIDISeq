import curses
from time import sleep
import mido
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


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

    def __init__(self, stdscr, midi_output, loop):
        self.loop = loop
        self.stdscr = stdscr
        self.output = midi_output
        self.output.send(mido.Message('program_change', program=19))


    def curses_render(self, t):
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                if x == t:
                    color = curses.color_pair(1)
                else:
                    color = curses.color_pair(2)

                if self.loop[x][y]:
                    self.stdscr.addch(21-y,x,curses.ACS_DIAMOND, color)
                else:
                    self.stdscr.addch(21-y,x,' ', color)
        self.stdscr.refresh()


    def midi_messages(self, t, scale):

        for i in range(0,len(self.loop[t])):
            if self.loop[t][i]:
                if not self.loop[t-1][i]:
                    self.output.send(mido.Message('note_on', note=scale[i], velocity=64))
            else:
                if self.loop[t-1][i]:
                    self.output.send(mido.Message('note_off', note=scale[i]))




def main(stdscr):
    # initialize curses environment
    curses.curs_set(0)
    curses.init_pair(1, curses.COLOR_RED, curses.COLOR_WHITE)
    curses.init_pair(2, curses.COLOR_WHITE, curses.COLOR_BLUE)
    stdscr.nodelay(1)

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
        [0, 0, 0, 0, 0, 0, 0, 1],
    ]

    a_loop = loop_from_dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    a_loop = loop_from_dna("GCTGGTAACCTTTTTTACTAGGTAGAGAAGCTGGACCAACTGGGGTTCTTCCAGGGGGAGAATGAGAAAGAGAAACTGTTTTGCAAGTCCGTAGCTATTTCTCTAGGGCCCTGTTAGCTGACATTGACATGCCTTGCATTGCTCTGCAGATCCCCTCGCAGCCCTCTGTCCCTTGTTCATTTCTGGCCTTAGAGAAAGCAAAGCAGGGTCTGTAACAGGGGAGGCTGCCTCTAAACTCAGGGTTTGGTTACAGCTGTTTTCACTTACATCACTGGCCCTGGTTTTTTTTTTTTTTCTGGCATTAAAAAAAAAAATTGGAAGCAGGTGATGTTCCCATTGCTGATGTGGTGGAAACTCTCCAAGTGAACAATATACGTTTTTCTTGGCAGCTGTTTCTTGTGCCCTGCTTGCTCCTGGTCCAGGACAAGCAAGGACCATCTGCCTCTTTCAATAGAACACCTCCAGATCCCTTTGATCAAAAGTTACTCATTGTCTGACTTGCTATTTCTGTGAGATAAATGGGAGAAGATCAATAAATGCACTTGTTTGTCCA")

    # create a loop object
    l = loop(stdscr, 
             mido.open_output( u'TiMidity port 0'), 
             a_loop)

    # initialize main loop
    bpm = 80
    delay = (60.0 / bpm) / 4
    t = 0
    while True:
        l.midi_messages(t, scale)
        l.curses_render(t) 


        # as time goes by
        sleep(delay)
        t+=1
        if t == len(l.loop):
            t = 0




curses.wrapper(main)

