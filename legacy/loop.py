import pprint
from time import sleep
import mido
import random

output = mido.open_output( u'TiMidity port 0')
#woodblock
#output.send(mido.Message('program_change', program=115))

output.send(mido.Message('program_change', program=80))


loop = [
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 1, 0, 1, 0],
    [1, 0, 1, 1, 0, 0, 0, 0],
    [0, 1, 0, 0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 1, 0, 1],
    [1, 0, 0, 0, 0, 1, 0, 1],
    [0, 1, 0, 1, 0, 0, 0, 1],
    [0, 0, 1, 0, 1, 0, 0, 1],
    [0, 1, 0, 1, 1, 0, 0, 0],
    [1, 0, 0, 0, 1, 0, 0, 0],
    ]

# C3 = 48
notas = range(48,80,2)

bpm = 111
delay = 60.0 / bpm

try:
    t = 0
    while True:
        sleep(delay)
        
        for i in range(0,len(loop[t])):
            if loop[t][i]:
                if not loop[t-1][i]:
                    output.send(mido.Message('note_on', note=notas[i], velocity=64))
            else:
                if loop[t-1][i]:
                    output.send(mido.Message('note_off', note=notas[i], velocity=64))

        print loop[t]

        # as time goes by
        t+=1
        if t == len(loop):
            t = 0
        

except KeyboardInterrupt:
    for i in range(0,len(loop[0])):
        output.send(mido.Message('note_off', note=notas[i], velocity=64))
