"""Q1b helper: decode a SAM FLAG into its individual bits.

Usage:
    python q1_decode_flag.py <flag>
    samtools view sorted.bam | sed -n 10p | cut -f2 | xargs python q1_decode_flag.py
"""
import sys

FLAG_BITS = [
    (0x1, "PAIRED", "template has multiple segments (read is paired)"),
    (0x2, "PROPER_PAIR", "each segment properly aligned (proper pair)"),
    (0x4, "UNMAP", "this read is unmapped"),
    (0x8, "MUNMAP", "the mate (next segment) is unmapped"),
    (0x10, "REVERSE", "this read is mapped to the reverse strand"),
    (0x20, "MREVERSE", "the mate is mapped to the reverse strand"),
    (0x40, "READ1", "this is the first read in the pair"),
    (0x80, "READ2", "this is the second read in the pair"),
    (0x100, "SECONDARY", "secondary (non-primary) alignment"),
    (0x200, "QCFAIL", "read fails platform/vendor quality checks"),
    (0x400, "DUP", "PCR or optical duplicate"),
    (0x800, "SUPPLEMENTARY", "supplementary (chimeric) alignment"),
]


def decode(flag):
    print(f"FLAG = {flag} = 0b{flag:012b}")
    for bit, name, meaning in FLAG_BITS:
        state = "SET  " if flag & bit else "unset"
        print(f"  0x{bit:03x} ({bit:4d}) {state} {name:<13} {meaning}")
    if flag & 0x1:
        which = "first" if flag & 0x40 else "second" if flag & 0x80 else "unknown"
        print(f"=> {which} read in the pair; mate is "
              f"{'UNMAPPED' if flag & 0x8 else 'mapped'}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python q1_decode_flag.py <flag>")
    decode(int(sys.argv[1], 0))
