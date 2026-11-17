"""Read and validate paired FASTQ records."""
from telescope_common import require
VALID_QUALITIES = bytes(range(33,127))

def fastq_record(f, mate):
    lines = [f.readline() for _ in range(4)]
    if not any(lines):
        return None
    require(all(lines), 'Truncated FASTQ record')
    h,seq,plus,qual = [s.rstrip(b'\r\n') for s in lines]
    require(h.startswith(b'@') and plus.startswith(b'+'), 'Invalid FASTQ header/separator')
    require(len(seq)==len(qual)>0, 'FASTQ sequence/quality length mismatch')
    require(not qual.translate(None, VALID_QUALITIES), 'Invalid FASTQ Phred+33 quality character')
    parts=h[1:].split(); name=parts[0]
    if name.endswith((b'/1',b'/2')):
        require(name.endswith(f'/{mate}'.encode()), 'Wrong mate suffix')
        name=name[:-2]
    if len(parts)>1 and parts[1][:2] in (b'1:',b'2:'):
        require(parts[1].startswith(f'{mate}:'.encode()), 'Wrong Illumina mate field')
    return name, b''.join(lines)
