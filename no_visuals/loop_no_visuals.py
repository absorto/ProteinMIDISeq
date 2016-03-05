# coding: utf8
import mido
from loop import Loop, loop_from_dna
from Bio import SeqIO
from Bio.Data import CodonTable
#import mingus.core.scales as scales
#from mingus.containers import Note
from pprint import pprint
from time import sleep

from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilename



def read_seq(*args):
    fname = askopenfilename()
    seq.set(fname)

def play(midi_output, path):

    loop = loop_from_dna(SeqIO.read(path,'fasta').seq.tostring())

    pprint(loop)
    a_loop = Loop( loop,
                   A,
                   midi_output=midi_output)

    bpm = 120
    delay = 60.0 / bpm

    while True:
        sleep(delay)
        a_loop.print_head()
        a_loop.play_at_head()


def play_callback(*args):

    # initialize midi output
    midi_output = mido.open_output( port.get() )
    #midi_output.send(mido.Message('program_change', program=1))

    play(midi_output,
         seq.get())

    midi_output.close()

    
root = Tk()
root.title("DNA sequence to midi events")

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

port = StringVar()
meters = StringVar()
seq = StringVar()

# MIDI PORT
ttk.Label(mainframe, text="MIDI port").grid(column=1, row=1, sticky=E)
mido_out = ttk.Combobox(mainframe, width=30, textvariable=port, values=mido.get_output_names(), state='readonly')
mido_out.set(mido.get_output_names()[0])
mido_out.grid(column=2, row=1, sticky=(W, E))

# Sequence
ttk.Button(mainframe, text="load DNA sequence", command=read_seq).grid(column=1, row=2,  sticky=E)
ttk.Label(mainframe, text="", textvariable=seq).grid(column=2, row=2, sticky=W)

# scale
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
row = 3
for amino in list(standard_table.back_table.keys()):
    ttk.Label(mainframe, text=amino).grid(column=1, row=row, sticky=E)
    row+=1

# play!
ttk.Button(mainframe, text="Play", command=play_callback).grid(column=2, row=4,  sticky=E)


for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

mido_out.focus()
root.bind('<Return>', play_callback)

root.mainloop()



# A = []
# for octave in range(3,6):
#     for tone in scales.NaturalMinor("C").ascending():
#         A.append( int(Note(tone, octave) ) )
# A.sort()


