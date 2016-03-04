# coding: utf8
import mido
from loop import Loop, loop_from_dna
from Bio import SeqIO
#import mingus.core.scales as scales
#from mingus.containers import Note
from pprint import pprint
from time import sleep

from tkinter import *
from tkinter import ttk



def calculate(*args):
    try:
        value = float(feet.get())
        meters.set((0.3048 * value * 10000.0 + 0.5)/10000.0)
    except ValueError:
        pass


def play(*args):

    # initialize midi output
    midi_output = mido.open_output( port.get() )

    #midi_output.send(mido.Message('program_change', program=1))

    pprint(seq.get(1.0, END))
    #aboseq = SeqIO.basestring(seq.get())
    #abo_loop = loop_from_dna(aboseq)

    #pprint(abo_loop)
    # a_loop = Loop( abo_loop,
    #                A,
    #                midi_output=midi_output)

    # bpm = 120
    # delay = 60.0 / bpm

    # while True:
    #     sleep(delay)
    #     a_loop.print_head()
    #     a_loop.play_at_head()

    

    
root = Tk()
root.title("DNA sequence to midi events")

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

port = StringVar()
meters = StringVar()

ttk.Label(mainframe, text="MIDI port").grid(column=1, row=1, sticky=E)
mido_out = ttk.Combobox(mainframe, width=30, textvariable=port, values=mido.get_output_names(), state='readonly')
mido_out.set(mido.get_output_names()[0])
mido_out.grid(column=2, row=1, sticky=(W, E))

seq = Text(mainframe, width=150,height=40).grid(column=1, row=2, columnspan=2,sticky=(W, E))
ttk.Button(mainframe, text="Play", command=play).grid(column=2, row=3,  sticky=E)

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

mido_out.focus()
root.bind('<Return>', calculate)

root.mainloop()



# A = []
# for octave in range(3,6):
#     for tone in scales.NaturalMinor("C").ascending():
#         A.append( int(Note(tone, octave) ) )
# A.sort()


