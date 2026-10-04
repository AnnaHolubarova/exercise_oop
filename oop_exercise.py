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
# Task 3: Gene and Variant subclasses
# ==========================================
class Gene(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, name):
        super().__init__(chromosome, start, end, strand)
        self.name = name
        self.exons = []
        
    def add_exon(self, exon):
        self.exons.append(exon)
        
    def total_exon_length(self):
        return sum(exon.length() for exon in self.exons)
        
    def describe(self):
        return f"{super().describe()} {self.name}, {len(self.exons)} exon(s)"


class Variant(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, ref_allele, alt_allele):
        super().__init__(chromosome, start, end, strand)
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        
    def variant_type(self):
        if len(self.ref_allele) == 1 and len(self.alt_allele) == 1:
            return "SNP"
        elif len(self.alt_allele) > len(self.ref_allele):
            return "insertion"
        elif len(self.alt_allele) < len(self.ref_allele):
            return "deletion"
        else:
            return "MNV"
            
    def describe(self):
        return f"{super().describe()} {self.ref_allele}>{self.alt_allele} ({self.variant_type()})"


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


def run_task3():
    print("--- Task 3: Processing oop_data.tsv ---")
    genes = {}
    report_features = [] # To hold Genes and Variants for polymorphic iteration

    try:
        with open("oop_data.tsv", "r") as f:
            next(f)  # Skip the header row
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
        
                if len(parts) < 7:
                    parts += [""] * (7 - len(parts))
                    
                ftype, chrom, start, end, strand, field_a, field_b = parts
                start, end = int(start), int(end)

                if ftype == "gene":
                    gene = Gene(chrom, start, end, strand, field_a)
                    genes[field_a] = gene
                    report_features.append(gene)
                elif ftype == "exon":
                    parent_name, exon_number = field_a, field_b
                    exon = Exon(chrom, start, end, strand, exon_number)
                    if parent_name in genes:
                        genes[parent_name].add_exon(exon)
                elif ftype == "variant":
                    ref, alt = field_a, field_b
                    variant = Variant(chrom, start, end, strand, ref, alt)
                    report_features.append(variant)

    except FileNotFoundError:
        print("Error: 'oop_data.tsv' not found. Ensure it is in the same directory.")
        return

    # Polymorphic Report
    print("Polymorphic report (genes and variants):")
    for feature in report_features:
        print(f"- {feature.describe()}")
        # Check specifically for gene to print total exon length next to it
        if isinstance(feature, Gene):
            print(f"Total exon length: {feature.total_exon_length()}")
            
    # Variant Location Analysis
    print("\nVariant context report:")
    variants = [f for f in report_features if isinstance(f, Variant)]
    
    for variant in variants:
        overlapping_genes = []
        for gene in genes.values():
            if variant.overlaps(gene):
                # Check if it also falls inside one of the gene's exons
                inside_exon = any(variant.overlaps(exon) for exon in gene.exons)
                overlapping_genes.append((gene.name, inside_exon))
                
        print(f"Variant: {variant.describe()}")
        if not overlapping_genes:
            print("  -> Context: intergenic")
        else:
            for gene_name, inside_exon in overlapping_genes:
                if inside_exon:
                    print(f"  -> Context: inside {gene_name} (in exon)")
                else:
                    print(f"  -> Context: inside {gene_name} (intron/UTR)")

 
if __name__ == "__main__":
    run_task1_and_2()
    run_task3()