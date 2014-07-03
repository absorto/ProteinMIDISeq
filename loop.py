import pprint
from time import sleep
import mido
import random

output = mido.open_output( u'TiMidity port 0')
#woodblock
#output.send(mido.Message('program_change', program=115))

output.send(mido.Message('program_change', program=80))


loop = [
    [1, 0, 1, 0, 1, 0, 0, 0],
    [0, 1, 0, 1, 0, 1, 0, 1],
    [0, 1, 0, 1, 0, 1, 0, 1],
    [0, 0, 1, 0, 1, 0, 1, 0],
    [1, 0, 1, 0, 1, 0, 0, 0],
    [1, 0, 1, 0, 1, 0, 0, 0],
    [0, 1, 0, 1, 0, 1, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0],
    ]

# C3 = 48
notas = range(48,80,2)

bpm = 111
delay = 60.0 / bpm

try:
    n = 0
    while True:
        sleep(delay)
        
        for i in range(0,len(loop[n])):
            if loop[n][i]:
                if not loop[n-1][i]:
                    output.send(mido.Message('note_on', note=notas[i], velocity=64))
            else:
                if loop[n-1][i]:
                    output.send(mido.Message('note_off', note=notas[i], velocity=64))

        print loop[n]

        # increment or restart loop
        n+=1
        if n==len(loop):
            n = 0
            # output.send(mido.Message('program_change', program=random.randrange(70,80)))
        

except KeyboardInterrupt:
    for i in notas:
        output.send(mido.Message('note_on', note=notas[i], velocity=64))