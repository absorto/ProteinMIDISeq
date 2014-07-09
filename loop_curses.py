import curses
from time import sleep
import mido
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
print coding_dna.translate()


class loop:

    def __init__(self, stdscr, midi_output, loop):
        self.loop = loop
        self.stdscr = stdscr
        self.output = midi_output
        self.output.send(mido.Message('program_change', program=80))


    def curses_render(self, t):
        for x in range(0,len(self.loop)):
            for y in range(0,len(self.loop[0])):
                if x == t:
                    color = curses.A_REVERSE
                else:
                    color = curses.color_pair(1)

                if self.loop[x][y]:
                    self.stdscr.addstr(y,x,'#', color)
                else:
                    self.stdscr.addstr(y,x,' ', color)
        self.stdscr.refresh()


    def midi_messages(self, t):
        # C3 = 48
        notas = range(48,80,2)

        for i in range(0,len(self.loop[t])):
            if self.loop[t][i]:
                if not self.loop[t-1][i]:
                    self.output.send(mido.Message('note_on', note=notas[i], velocity=64))
            else:
                if self.loop[t-1][i]:
                    self.output.send(mido.Message('note_off', note=notas[i], velocity=64))


def main(stdscr):
    # initialize curses environment
    curses.curs_set(0)
    curses.init_pair(1, curses.COLOR_RED, curses.COLOR_WHITE)
    stdscr.nodelay(1)

    a_loop = [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
    ]
    # create a loop object
    l = loop(stdscr, 
             mido.open_output( u'TiMidity port 0'), 
             a_loop)

    # initialize main loop
    bpm = 111
    delay = 60.0 / bpm
    t = 0
    while True:
        l.curses_render(t)        
        l.midi_messages(t)

        # as time goes by
        sleep(delay)
        t+=1
        if t == len(l.loop):
            t = 0




curses.wrapper(main)



#         # get keyboard input, returns -1 if none available
#         c = stdscr.getch()
#         if c != -1:
#             if c == curses.KEY_DOWN:
#                 y+=1
            
#             if c == curses.KEY_UP:
#                 y-=1
            
#             if c == curses.KEY_LEFT:
#                 x-=1
            
#             if c == curses.KEY_RIGHT:
#                 x+=1

#             # print numeric value
#             stdscr.addstr(y,x,'#', curses.color_pair(1))

# #            stdscr.addstr(y,x,'#', curses.A_STANDOUT)
#             stdscr.refresh()
#             # return curser to start position
#             stdscr.move(0, 0)
