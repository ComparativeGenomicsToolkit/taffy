
class Wiggle:
    """
    Lazy implementation of a wiggle parser
    """
    def __init__(self, file):
        self.seq_intervals = {}
        with open(file, "r") as f:
            line = f.readline()
            while True:
                header = self.parse_header(line)
                seq_name = header["chrom"]
                if seq_name not in self.seq_intervals:
                    self.seq_intervals[seq_name] = {}
                values = self.seq_intervals[seq_name]
                span = int(header["span"]) if "span" in header else 1
                if header["fixed_step"]:
                    step = int(header["step"]) if "step" in header else 1
                    assert span <= step
                    i = int(header["start"])
                    for line in f:
                        tokens = line.split()
                        if len(tokens) > 1:
                            break
                        if len(tokens) == 1:
                            v = float(tokens[0])
                            for j in range(span):
                                values[i+j] = v
                            i += step
                else:
                    for line in f:
                        tokens = line.split()
                        if len(tokens) == 0:
                            continue
                        if len(tokens) > 2 or tokens[0] in ("fixedWith", "variableStep"):
                            break
                        if len(tokens) == 2:
                            i = int(tokens[0])
                            v = float(tokens[1])
                            for j in range(span):
                                values[i+j] = v

    def get(self, seq_name, index):
        return self.seq_intervals[seq_name][index]

    @staticmethod
    def parse_header(line):
        tags = {}
        tokens = line.split()
        step = tokens[0]
        if step not in ("fixedStep", "variableStep"):
            raise RuntimeError(f"Misformed wiggle header line: {line}")
        tags["fixed_step"] = step == "fixedStep"
        for tag in tokens[1:]:
            if len(tag.split("=")) != 2:
                raise RuntimeError(f"Misformed wiggle header line: {line}")
            key, value = tag.split("=")
            tags[key] = value
        return tags
