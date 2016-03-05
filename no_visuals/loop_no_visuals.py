# coding: utf8
import mido
from loop import Loop, loop_from_dna
from Bio import SeqIO
from Bio.Data import CodonTable

from pprint import pprint
from time import sleep

from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilename



def read_seq(*args):
    fname = askopenfilename()
    seq.set(fname)

def play(midi_output, path, scale):

    loop = loop_from_dna(SeqIO.read(path,'fasta').seq.tostring())

    pprint(loop)
    a_loop = Loop( loop,
                   scale,
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

    scale = [int(g.get()) for g in wscale]
    
    play(midi_output,
         seq.get(),
         scale)

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
C_major = [36,38,40,41,43,45,47,
           48,50,52,53,55,57,59,
           60,62,64,65,67,69,71]
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
aminoacids = []
for a in list(standard_table.back_table.keys()):
    if a:
        aminoacids.append(a)
aminoacids.sort()
ttk.Label(mainframe, text='amino').grid(column=1, row=3, sticky=E)
ttk.Label(mainframe, text='mini note').grid(column=2, row=3, sticky=W)
row = 4
wscale = []
for amino in aminoacids:
    ttk.Label(mainframe, text=amino).grid(column=1, row=row, sticky=E)
    wscale.append(StringVar())
    tmp = ttk.Entry(mainframe, width=4, textvariable=wscale[len(wscale)-1])
    tmp.grid(column=2, row=row, sticky=W)
    tmp.insert(0,string=C_major.pop(0))
    row+=1


# play!
ttk.Button(mainframe, text="Play", command=play_callback).grid(column=3, row=2,  sticky=E)


for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

mido_out.focus()
root.bind('<Return>', play_callback)

root.mainloop()



# A = []
# for octave in range(3,6):
#     for tone in scales.NaturalMinor("C").ascending():
#         A.append( int(Note(tone, octave) ) )
# A.sort()


