"""Bisection-based random access for sorted coverage files."""
import os


def _chrconv(chr_):
    return chr_ if "MT" in chr_ else "chr" + chr_.replace("chr", "")


def _index_ucsc(filename):
    chr_list = []
    with open(filename) as f:
        for line in f:
            chr_ = line.split()[0] if (line[0] == "c" or line[0] == "M") else line.split()[2]
            if chr_ not in chr_list:
                chr_list.append(chr_)
    return chr_list


class Bisect:
    """Bisection search over a position-sorted text file."""

    _DEFAULT_CHRS = (
        "chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
        "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chrM", "chrX", "chrY", "MT",
    )

    def __init__(self, fname, feeder=None):
        self.fname = fname
        self.size = os.path.getsize(fname)
        self.f = open(self.fname)
        self.chr_list = list(self._DEFAULT_CHRS)
        self.seek_chr = ""
        self.curr_chr = "chr1"
        self.pos = 0
        if feeder:
            self.chr_list = _index_ucsc(feeder)

    def create_index(self):
        new_chr_list = []
        self.f.seek(0)
        self.f.read(100)
        while True:
            chunk = self.f.read(5000)
            if len(chunk) <= 1:
                break
            lines = chunk.split("\n")
            for line in (lines[1], lines[-2]):
                try:
                    chr_ = _chrconv(line.split()[0])
                    if chr_ in self.chr_list and chr_ not in new_chr_list:
                        new_chr_list.append(chr_)
                except (IndexError, ValueError):
                    print("End_of_File. Continuing")
        self.chr_list = new_chr_list
        self.f.seek(0)

    def _which_half(self, chr_):
        try:
            if self.chr_list.index(chr_) > self.chr_list.index(self.seek_chr):
                return -1
            return 1
        except ValueError:
            print(f"Bad chromosome: {chr_}")
            return 1

    def find(self, chr_, start, save_pos=False):
        s = int(start)
        f = self.f
        f.seek(0)
        DIST = 1e0
        SMALL_JUMP = 1e3

        try:
            self.seek_chr = self.chr_list[self.chr_list.index(chr_)]
        except ValueError:
            print(f"{chr_} not in my list - skipped")
            return None

        if save_pos and self.seek_chr == self.curr_chr:
            pos = self.pos
            inc = SMALL_JUMP
        else:
            pos = self.size / 2
            inc = self.size / 4

        while True:
            pos = max(pos, 0)
            f.seek(int(pos))
            f.readline()
            _tmp_chr = ""
            while _tmp_chr not in self.chr_list:
                try:
                    _tmp_chr, _tmp_s, _tmp_e = f.readline().split()[:3]
                    _tmp_chr = _chrconv(_tmp_chr)
                except ValueError:
                    break

            self.curr_chr = _tmp_chr

            if _tmp_chr == self.seek_chr:
                if (int(_tmp_s) + DIST) < s:
                    while (int(_tmp_s) + DIST) < s and _tmp_chr == self.seek_chr:
                        pos += SMALL_JUMP
                        f.seek(int(pos))
                        f.readline()
                        try:
                            _tmp_chr, _tmp_s, _tmp_e = f.readline().split()[:3]
                            _tmp_chr = _chrconv(_tmp_chr)
                        except ValueError:
                            pos -= SMALL_JUMP
                            f.seek(int(pos))
                            f.readline()
                            print("Ouch!! It is the end of file")
                            self.pos = f.tell()
                            return self.f
                    pos -= SMALL_JUMP
                    pos = max(pos, 0)
                    f.seek(int(pos))
                    f.readline()
                    break

                if (int(_tmp_s) + DIST) >= s:
                    while int(_tmp_s) > s and _tmp_chr == self.seek_chr:
                        pos -= SMALL_JUMP
                        pos = max(pos, 0)
                        f.seek(int(pos))
                        if pos <= 0:
                            return self.f
                        f.readline()
                        _tmp_chr, _tmp_s, _tmp_e = f.readline().split()[:3]
                        _tmp_chr = _chrconv(_tmp_chr)
                    while _tmp_chr != self.seek_chr:
                        try:
                            _tmp_chr, _tmp_s, _tmp_e = f.readline().split()[:3]
                            _tmp_chr = _chrconv(_tmp_chr)
                        except ValueError:
                            break
                    break

            pos = pos + self._which_half(_tmp_chr) * inc
            inc = inc / 2
            if inc == 0:
                break

        self.pos = f.tell()
        return self.f

    def close(self):
        self.f.close()


class Regions(Bisect):
    """Random-access region extractor over a sorted coverage/BED file."""

    def __init__(self, fname):
        super().__init__(fname)
        self.create_index()
        self.mem_file = []

    def focus_single_prl(self, pos, step=1):
        """Yield all entries between (chr, start, end), thread-safe."""
        chr_, start, end = pos
        self.mem_file = []
        chr_ = _chrconv(chr_)
        f = self.find(chr_, int(start))
        count = 0
        value = None
        try:
            f.readline()
            while True:
                line = f.readline()
                fpos = f.tell()
                count += 1
                chr_r, e_start = line.split()[:2]
                chr_r = _chrconv(chr_r)
                if chr_ == chr_r:
                    if int(start) <= int(e_start) <= int(end):
                        if count % step == 0:
                            value = line.strip()
                            yield value
                            f.seek(fpos)
                    if int(e_start) > end:
                        return
                else:
                    if value:
                        return
        except Exception:
            pass

    def into(self, pos):
        chr_, start, end = pos
        start, end = int(start), int(end)
        off = (end - start) / 2
        for line in self.mem_file:
            chr_r, e_start, e_end = line.split()[:3]
            chr_r = _chrconv(chr_r)
            if chr_ == chr_r:
                if start >= int(e_start) - off and end <= int(e_end) + off:
                    return 1
                if int(e_start) > end:
                    self.mem_file = self.mem_file[self.mem_file.index(line) - 100:]
                    return 0
        return 0
