# C3 = 48

# just a range of numbers with corresponding midi notes
scale_range = range(38,98,2)

# persian scale, using midi numbers
persian = [64+(1*0), 64+(1*1), 64+(1* 4), 64+(1*  5), 64+(1* 6), 64+(1* 8), 64+(1* 11), 64+(2*0), 64+(2*1), 64+(2* 4), 64+(2*  5), 64+(2* 6), 64+(2* 8), 64+(2* 11), 64+(3*0), 64+(3*1), 64+(3* 4), 64+(3*  5), 64+(3* 6), 64+(3* 8), 64+(3* 11)]


# c major scale using ABC notation
c_major_names = ['C0', 'D0', 'E0', 'F0', 'G0', 'A0', 'B0', 'C1', 'D1', 'E1', 'F1', 'G1', 'A1', 'B1', 'C2', 'D2', 'E2', 'F2', 'G2', 'A2', 'B2', 'C3', 'D3', 'E3', 'F3']


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


name2midi = NoteConverter()

c_major = [name2midi(n) for n in c_major_names]

