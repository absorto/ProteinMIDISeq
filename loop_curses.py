import curses

def main(stdscr):
    x = 5
    y = 5
    # undisplay cursor
    curses.curs_set(0)

    curses.init_pair(1, curses.COLOR_RED, curses.COLOR_WHITE)

    # do not wait for input when calling getch
    stdscr.nodelay(1)
    while True:
        # get keyboard input, returns -1 if none available
        c = stdscr.getch()
        if c != -1:
            if c == curses.KEY_DOWN:
                y+=1
            
            if c == curses.KEY_UP:
                y-=1
            
            if c == curses.KEY_LEFT:
                x-=1
            
            if c == curses.KEY_RIGHT:
                x+=1

            # print numeric value
            stdscr.addstr(y,x,curses.keyname(c), curses.color_pair(1))

            stdscr.addstr(0,0,curses.keyname(stdscr.inch(y,x)), curses.A_STANDOUT)

#            stdscr.addstr(y,x,'#', curses.A_STANDOUT)
            stdscr.refresh()
            # return curser to start position
            stdscr.move(0, 0)

if __name__ == '__main__':
    curses.wrapper(main)