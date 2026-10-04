# ==========================================
# Task 1: GenomicFeature base class
# ==========================================
class GenomicFeature:
    def __init__(self, chromosome, start, end, strand):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("Start and end positions must be integers.")
        if start < 1 or end < 1:
            raise ValueError("Start and end positions must be 1-based positive integers.")
        if start > end:
            raise ValueError(f"Invalid range: start ({start}) must be <= end ({end}).")
        if strand not in ("+", "-"):
            raise ValueError(f"Invalid strand: '{strand}'. Must be '+' or '-'.")
            
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        
    def length(self):
        return self.end - self.start + 1
        
    def overlaps(self, other):
        if self.chromosome != other.chromosome:
            return False
        # Overlap happens if one feature starts before or when the other ends, 
        # and ends after or when the other starts.
        return self.start <= other.end and self.end >= other.start
        
    def describe(self):
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand})"


# ==========================================
# Task 2: Exon subclass
# ==========================================
class Exon(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, exon_number):
        super().__init__(chromosome, start, end, strand)
        self.exon_number = int(exon_number)
        
    def describe(self):
        return f"{super().describe()} exon #{self.exon_number}"


# ==========================================
# Main Execution Blocks
# ==========================================
def run_task1_and_2():
    print("--- Task 1 Tests ---")
    a = GenomicFeature("chr1", 1000, 5000, "+")
    b = GenomicFeature("chr1", 4800, 6000, "+")
    c = GenomicFeature("chr2", 1000, 5000, "+")

    print(a.describe())     # GenomicFeature chr1:1000-5000(+)
    print(a.length())       # 4001
    print(a.overlaps(b))    # True
    print(a.overlaps(c))    # False
    try:
        GenomicFeature("chr1", 5000, 1000, "+")  # should raise ValueError
    except ValueError as e:
        print(f"Caught expected ValueError: {e}")

    print("\n--- Task 2 Tests ---")
    features = [
        GenomicFeature("chr1", 1000, 5000, "+"),
        Exon("chr1", 1000, 1200, "+", 1),
        Exon("chr1", 3000, 3300, "+", 2),
    ]
    for feature in features:
        print(feature.describe())
    print()

 
if __name__ == "__main__":
    run_task1_and_2()